#!/soft/R/R-3.0.0/bin/Rscript

library(ggplot2)
suppressMessages(library(RDAVIDWebService)) #Avoid the annoying message

plotStuff <- function(tag, pngName, p){
  png(paste0(out, "iLoops/", tag, "/", pngName, ".png"), width=2500, height=2500)
  print(p)
  graphics.off()
}

args <- commandArgs(trailingOnly = TRUE)
load(paste0(args[1], "RWorkspaces/2_GetCandidates.RData"))
tag <- args[2]

InteraX <- read.table( paste0(out, "iLoops/InteraXChanges_", tag, ".tsv"), header=T, sep="\t", stringsAsFactors=F)

txExpressionDifference = log(interReplicate[["N"]]$Median_TPM + 0.0001)  - log(interReplicate[["T"]]$Median_TPM + 0.0001)
gnExpressionDifference = log(interReplicate[["N"]]$Median_tTPM + 0.0001) - log(interReplicate[["T"]]$Median_tTPM + 0.0001)
psiDiff = interReplicate[["N"]]$Median_PSI  - interReplicate[["T"]]$Median_PSI

diffs_df <- data.frame(
                        tx=interReplicate[["N"]]$Transcript, 
                        Entrez=interReplicate[["N"]]$Gene, 
                        nTxExp=log(interReplicate[["N"]]$Median_TPM + 0.0001), 
                        tTxExp=log(interReplicate[["T"]]$Median_TPM + 0.0001), 
                        nGnExp=log(interReplicate[["N"]]$Median_tTPM + 0.0001),
                        tGnExp=log(interReplicate[["T"]]$Median_tTPM + 0.0001),
                        nPSI=interReplicate[["N"]]$Median_PSI,
                        tPSI=interReplicate[["T"]]$Median_PSI
                      )

InteraX <- merge(InteraX, diffs_df, by.x = "Partner", by.y = "tx")
InteraX$dTxExpression <- InteraX$nTxExp - InteraX$tTxExp
InteraX$dGnExpression <- InteraX$nGnExp - InteraX$tGnExp
InteraX$dPSI <- InteraX$nPSI - InteraX$tPSI

write.table(InteraX, paste0(out, "iLoops/InteraXChanges_", tag, "_full1.tsv"), sep="\t", row.names=F, quote=F)
save(InteraX, file=paste0(out, "RWorkspaces/5_Interactions_", tag, "_1.RData"))

plotStuff(tag, paste0(tag, "_txExpression"), ggplot(InteraX, aes(x=dRC, y=dTxExpression, color=Annotation)) + geom_point(shape=1) + geom_text(data=InteraX[InteraX$Annotation == "Driver",], aes(label=Partner_gene),hjust=0, vjust=0) )
plotStuff(tag, paste0(tag, "_gnExpression"), ggplot(InteraX, aes(x=dRC, y=dGnExpression, color=Annotation)) + geom_point(shape=1) + geom_text(data=InteraX[InteraX$Annotation == "Driver",], aes(label=Partner_gene),hjust=0, vjust=0) )
plotStuff(tag, paste0(tag, "_psi"), ggplot(InteraX, aes(x=dRC, y=dPSI, color=Annotation)) + geom_point(shape=1) + geom_text(data=InteraX[InteraX$Annotation == "Driver",], aes(label=Partner_gene),hjust=0, vjust=0) )

allExpressedTranscripts <- as.character()
for (tx in intraReplicate){
  mask <- tx$TPM_N >= 0.1 | tx$TPM_T >= 0.1
  allExpressedTranscripts <- c(allExpressedTranscripts, tx$Gene[mask])
}

allExpressedTranscripts <- unique(allExpressedTranscripts)
expressedEntrezIds <- unlist(strsplit(as.character(allExpressedTranscripts[grepl("|", allExpressedTranscripts, fixed=TRUE)]), "\\|"))[c(F,T)]
tmp <- as.character(InteraX$Entrez)
tmp[!grepl("\\|", tmp)] <- paste0("xxx|", tmp[!grepl("\\|", tmp)])
entrezIds <- unlist(strsplit(as.character(tmp), "\\|"))[c(F,T)]

if (length(entrezIds) >= 3000){
  cat("More than 3000 genes. Abort. \n")
  quit(save = "no")
}

david <- DAVIDWebService$new(email="hector.climente@upf.edu")
bg <- addList(david, expressedEntrezIds,idType="ENTREZ_GENE_ID", listName="Expressed", listType="Background")
query <- addList(david, entrezIds, idType="ENTREZ_GENE_ID", listName=paste0(tag, "_Partners"), listType="Gene")

setAnnotationCategories(david, "GOTERM_BP_ALL")
termCluster <- getClusterReport(david, type="Term")

clusterComposition <- list()

for (clusterNum in seq(1, length(enrichment(termCluster)))){
  
  uniqGenes <- unique(unlist(strsplit(members(termCluster)[[clusterNum]]$Genes, ", ")))
  clusterComposition[[clusterNum]] <- entrezIds %in% uniqGenes
  if (min(members(termCluster)[[clusterNum]][,5]) > 0.1){
    cat("Couldn't process cluster", clusterNum, "\n")
    next
  }
  
  plotStuff(tag, paste0(tag, "_cluster", as.character(clusterNum) ), plot2D(termCluster, clusterNum) )
  davidGODag <- DAVIDGODag(members(termCluster)[[clusterNum]], "BP", pvalueCutoff=0.1, removeUnattached=TRUE)
  plotStuff(tag, paste0(tag, "_cluster", as.character(clusterNum), "_graph"), plotGOTermGraph(g=goDag(davidGODag), r=davidGODag, max.nchar=30, node.shape="ellipse"))
  
}

InteraX <- cbind(InteraX, do.call("cbind", clusterComposition))
write.table(InteraX, paste0(out, "iLoops/InteraXChanges_", tag, "_full2.tsv"), sep="\t", row.names=F, col.names=F, quote=F)
save(InteraX, termCluster, file=paste0(out, "RWorkspaces/5_Interactions_", tag, "_2.RData"))