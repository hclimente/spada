#!/soft/R/R-3.0.0/bin/Rscript

library(plyr)

args <- commandArgs(trailingOnly = TRUE)
minExpression <- as.numeric(args[1])
minCandidateExpression <- as.numeric(args[2])
load(paste0("Results/", args[3], "/RWorkspaces/1_ExploreData.RData"))

candidates <- list()
allGenes <- as.character()

for (replicate in seq(1,inputData[["Replicates"]])){
  
	cat("\t* Replicate", replicate)

  candidates[[replicate]] <- data.frame(Gene=as.character(), Genename=as.character(), Switch=as.numeric(), maxdPSI=as.character(), mindPSI=as.character())
  
  #Filter by deltaPSI and expression, based on the FPR
  #psiThreshold <- abs(intraReplicate[[replicate]]$deltaPSI) > 0.15
  psiThreshold <- abs(intraReplicate[[replicate]]$deltaPSI) > 4 * interReplicate$MAD
  psiThreshold[is.na(psiThreshold)] <- FALSE
  posPSI <- intraReplicate[[replicate]]$deltaPSI > 0 & psiThreshold
  negPSI <- intraReplicate[[replicate]]$deltaPSI < 0 & psiThreshold
  expressionThreshold <- intraReplicate[[replicate]]$la_tTPM > minCandidateExpression
  
  replicateCandidates <- vector('list', length(unique(intraReplicate[[replicate]]$Gene[psiThreshold & expressionThreshold])))

  for (aCandidate in unique(intraReplicate[[replicate]]$Gene[psiThreshold & expressionThreshold])){
    
    thisGeneData <- intraReplicate[[replicate]]$Gene == aCandidate

    if(length(intraReplicate[[replicate]]$deltaPSI[thisGeneData & posPSI]) == 0 || length(intraReplicate[[replicate]]$deltaPSI[thisGeneData & negPSI]) == 0){
      next
    }
    
    #Calculate max difference between the transcripts
    maxDeltaPsi <- max(intraReplicate[[replicate]]$deltaPSI[thisGeneData & posPSI], na.rm=T)
    minDeltaPsi <- min(intraReplicate[[replicate]]$deltaPSI[thisGeneData & negPSI], na.rm=T)
    maxSwitch <- maxDeltaPsi - minDeltaPsi
    
    #deltaPSI = PSI_ref - PSI_alt
    #	MaxDeltaPsiCond: predominant transcript in the normal sample
    #	MinDeltaPsiCond: predominant transcript in the tumor sample
    maxDeltaPsiCond <- intraReplicate[[replicate]]$deltaPSI == maxDeltaPsi
    minDeltaPsiCond <- intraReplicate[[replicate]]$deltaPSI == minDeltaPsi

    replicateCandidates[[aCandidate]] <- data.frame(Genename=intraReplicate[[replicate]]$Genename[thisGeneData & maxDeltaPsiCond], 
                                                    Gene=aCandidate,
                                                    maxdPSI=intraReplicate[[replicate]]$Transcript[thisGeneData & maxDeltaPsiCond], 
                                                    mindPSI=intraReplicate[[replicate]]$Transcript[thisGeneData & minDeltaPsiCond],
                                                    Switch=maxSwitch
                                                   )
  }

  candidates[[replicate]] <- do.call('rbind', replicateCandidates)
  
  #Expressed genes: transcript whose expression is above the threshold
  expressedGenes <-(log(as.numeric(intraReplicate[[replicate]]$TPM_N)) + log(as.numeric(intraReplicate[[replicate]]$TPM_T))) / 2 > minExpression
  allGenes <- unique(c(allGenes, intraReplicate[[replicate]]$Transcript[expressedGenes]))
  
  switchCut <- candidates[[replicate]]$Switch > 0.2
  
  candidates[[replicate]] <- candidates[[replicate]] [switchCut, c("Genename", "Gene", "maxdPSI","mindPSI", "Switch")]

  cat(":", nrow(candidates[[replicate]]), "candidates found\n")
  
}

candidateList <- do.call("rbind", candidates)
candidateList <- ddply(candidateList,.(Genename,Gene,maxdPSI,mindPSI), summarise, Replicated=length(Genename))
candidateList <- candidateList[with(candidateList, order(-Replicated)), ]

write.table(candidateList, file=paste0(out, "candidateList.tsv"), sep="\t", row.names=F, col.names=F, quote=F)
write(allGenes, paste0(out, "expressedGenes.lst"), sep="\n")

library(gplots)

top <- 30

topCandidates <- head(candidateList, n=top)
fig <- data.frame(matrix(nrow=top, ncol=inputData[["Replicates"]]))
rownames(fig) <- topCandidates$Genename
colnames(fig) <- seq(1,inputData[["Replicates"]])

for (replicate in seq(1,inputData[["Replicates"]])){
  for (gene in topCandidates$Genename){
    if (gene %in% candidates[[replicate]]$Genename) {
      fig[gene, replicate] <- candidates[[replicate]]$Switch[candidates[[replicate]]$Genename == gene]
    } else {
      fig[gene, replicate] <- NA
    }
  }
}

heatmap.2(as.matrix(fig), trace="none", scale="none", Rowv=NULL, Colv=NULL, 
          dendrogram="none", col=colorpanel(20, "white", "blue"), na.col="grey", 
          breaks=seq(0, max(fig, na.rm=T), length.out=21), main="PSI Switch"
         )

save(isoformExpression, intraReplicate, interReplicate, candidates, candidateList, inputData, wd, out, file=paste0(out, "RWorkspaces/2_GetCandidates.RData"))