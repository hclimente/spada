#!/soft/R/R-3.0.0/bin/Rscript

library(plyr)
suppressMessages( library(logging) )

args <- commandArgs(trailingOnly = TRUE)
minExpression <- as.numeric(args[1])
load(paste0(args[2], "RWorkspaces/1_ExploreData.RData"))

patientSet <- inputData$Replicates
if (length(inputData$unpairedReplicates) > 0) {
  patientSet <- inputData$unpairedReplicates
}

numOfReplicates <- length(patientSet)

logger <- getLogger(name="get_candidates", level=10) #Level debug

addHandler(writeToConsole, logger="get_candidates", level='INFO')
addHandler(writeToFile, logger="get_candidates", file=paste0(out, "rSmartAS.log"), level='DEBUG')

binomialTest <- function(x){

  successes <- as.numeric(x[4])
  a <- binom.test(successes,numOfReplicates,p = 0.05,alternative = "greater",conf.level = 0.95)
  return(unlist(a["p.value"]))

}

candidates <- list()
replicateExpressed  <- list()
candidateMask <- list()
candidateMask[["N"]] <- list()
candidateMask[["T"]] <- list()

loginfo("Searching isoform switches in %d patients.", numOfReplicates, logger="get_candidates")

candidatesPB <- txtProgressBar(min=1, max=numOfReplicates, initial = 1, style=3)
counter <- 1

for (replicate in patientSet){
  
  candidates[[replicate]] <- data.frame(Gene=as.character(), Switch=as.numeric(), maxdPSI=as.character(), mindPSI=as.character())
  
  #Filter by deltaPSI and expression, based on the FPR
  psiThreshold <- abs(intraReplicate[[replicate]]$deltaPSI) > interReplicate[["N"]]$FPR_5
  psiThreshold[is.na(psiThreshold)] <- FALSE
  
  norExpression <- log(intraReplicate[[replicate]]$TPM_N) > minExpression
  tumExpression <- log(intraReplicate[[replicate]]$TPM_T) > minExpression
  
  replicateCandidates <- vector('list', length(unique(intraReplicate[[replicate]]$Gene[psiThreshold]))) 
  
  for (aCandidate in unique(intraReplicate[[replicate]]$Gene[psiThreshold])){ 
    
    thisGeneData <- intraReplicate[[replicate]]$Gene == aCandidate
    
    norTranscripts <- intraReplicate[[replicate]]$deltaPSI > 0 & psiThreshold & norExpression & thisGeneData
    tumTranscripts <- intraReplicate[[replicate]]$deltaPSI < 0 & psiThreshold & tumExpression & thisGeneData
    
    if(length(intraReplicate[[replicate]]$deltaPSI[norTranscripts]) == 0 || length(intraReplicate[[replicate]]$deltaPSI[tumTranscripts]) == 0){
      next
    }
    
    #Calculate max difference between the transcripts and choose the candidates
    #     deltaPSI = PSI_ref - PSI_alt
    # MaxDeltaPsiCond: predominant transcript in the normal sample
    # MinDeltaPsiCond: predominant transcript in the tumor sample
    maxDeltaPsi <- max(intraReplicate[[replicate]]$deltaPSI[thisGeneData & norTranscripts], na.rm=T)
    norCandidate <- intraReplicate[[replicate]]$deltaPSI == maxDeltaPsi
    
    minDeltaPsi <- min(intraReplicate[[replicate]]$deltaPSI[thisGeneData & tumTranscripts], na.rm=T)
    tumCandidate <- intraReplicate[[replicate]]$deltaPSI == minDeltaPsi

    Switch <- maxDeltaPsi - minDeltaPsi
    
    if (Switch < 0.2){
      next
    }
    
    replicateCandidates[[aCandidate]] <- data.frame(Gene=aCandidate,
                                                    maxdPSI=intraReplicate[[replicate]]$Transcript[thisGeneData & norCandidate], 
                                                    mindPSI=intraReplicate[[replicate]]$Transcript[thisGeneData & tumCandidate],
                                                    Switch=Switch)
  }
  
  candidates[[replicate]] <- do.call('rbind', replicateCandidates)
  
  candidateMask[["N"]][[replicate]] <- interReplicate[["T"]]$Transcript %in% candidates[[replicate]]$maxdPSI
  candidateMask[["T"]][[replicate]] <- interReplicate[["T"]]$Transcript %in% candidates[[replicate]]$mindPSI
  
  #Expressed genes: transcripts whose expression are above the expression threshold in any of the conditions
  replicateExpressed[[replicate]] <- data.frame(Transcript=intraReplicate[[replicate]]$Transcript[norExpression | tumExpression],
                                                Gene=intraReplicate[[replicate]]$Gene[norExpression | tumExpression])
  logdebug("%d switches found at patient %s", nrow(candidates[[replicate]]), replicate, logger="get_candidates")
  setTxtProgressBar(candidatesPB, counter)
  counter <- counter + 1
  
}

close(candidatesPB)

allExpressedTranscripts <- unique(do.call('rbind', replicateExpressed))
allExpressedTranscripts <- allExpressedTranscripts[with(allExpressedTranscripts, order(Transcript)), ]

candidateMask[["N"]] <- do.call('cbind', candidateMask[["N"]])
candidateMask[["T"]] <- do.call('cbind', candidateMask[["T"]])

for (i in patientSet){candidates[[i]]$Origin = i}

candidateList <- do.call("rbind", candidates)
candidateList <- ddply(candidateList,.(Gene,maxdPSI,mindPSI), summarise, Replicated=length(Gene), Patients=paste(Origin, collapse = ",") )
candidateList <- candidateList[with(candidateList, order(-Replicated)), ]
candidateList$pval <- apply(candidateList,1,binomialTest)

write.table(candidateList, file=paste0(out, "candidateList.tsv"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(allExpressedTranscripts, paste0(out, "expressedGenes.lst"), sep="\t", row.names=F, col.names=F, quote=F)

cols <- c( paste0("TPM_",patientSet), paste0("PSI_",patientSet), paste0("tTPM_",patientSet), "FPR_4MAD" )
write.table(interReplicate[["N"]][,!(colnames(interReplicate$N) %in% cols)], file=paste0(out, "expression_normal.tsv"), sep="\t", row.names=F, quote=F)
write.table(interReplicate[["T"]][,!(colnames(interReplicate$T) %in% cols)], file=paste0(out, "expression_tumor.tsv"), sep="\t", row.names=F, quote=F)

save(isoformExpression, intraReplicate, interReplicate, candidates, candidateList, candidateMask, inputData, allExpressedTranscripts, wd, out, file=paste0(out, "RWorkspaces/2_GetCandidates.RData"))

#Plot heatmap
suppressMessages(library(gplots)) #Avoid the annoying message
library(RColorBrewer)

top <- length(candidateList$Gene[candidateList$Replicated >= patientSet * 0.2])

topCandidates <- ddply(candidateList,.(Gene), summarise, Replicated=sum(Replicated))
topCandidates <- topCandidates[with(topCandidates, order(-Replicated)), ]
topCandidates <- head(topCandidates, n=top)

fig <- data.frame(matrix(nrow=length(topCandidates$Gene), ncol=numOfReplicates))
rownames(fig) <- topCandidates$Gene
colnames(fig) <- patientSet

for (replicate in patientSet){
  for (gene in topCandidates$Gene){
    if (gene %in% candidates[[replicate]]$Gene) {
      fig[gene, replicate] <- head(candidates[[replicate]]$Switch[candidates[[replicate]]$Gene == gene],1)
    } else {
      fig[gene, replicate] <- 0
    }
  }
}

png(paste0(out, "DataExploration/topCandidateSwitch.png"), width=960, height=960)
myPalette <- colorRampPalette(c("white", "firebrick2"))(n = 14)
heatmap.2(as.matrix(fig), trace="none", scale="none", col=myPalette, na.col="grey", 
          breaks=seq(0, max(fig, na.rm=T), length.out=15), main="PSI Switch")
          
graphics.off()

save(isoformExpression, intraReplicate, interReplicate, candidates, candidateList, candidateMask, inputData, allExpressedTranscripts, wd, out, fig, file=paste0(out, "RWorkspaces/2_GetCandidates.RData"))