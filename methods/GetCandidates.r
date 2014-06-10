#!/soft/R/R-3.0.0/bin/Rscript

library(plyr)

args <- commandArgs(trailingOnly = TRUE)
minExpression <- as.numeric(args[1])
load(paste0("Results/", args[2], "/RWorkspaces/1_ExploreData.RData"))
unpairedReplicates <- as.numeric(args[3])
numOfReplicates <- inputData[["Replicates"]]

if (unpairedReplicates != 0) {  
  load(paste0("Results/", args[2], "/RWorkspaces/1_Tumor_Intereplicate.RData"))
  numOfReplicates <- unpairedReplicates
}

candidates <- list()
replicateExpressed  <- list()
candidateMask <- list()
candidateMask[["N"]] <- list()
candidateMask[["T"]] <- list()

for (replicate in seq(1,numOfReplicates)){
  
  cat("\t* Replicate", replicate)
  
  candidates[[replicate]] <- data.frame(Gene=as.character(), Genename=as.character(), Switch=as.numeric(), maxdPSI=as.character(), mindPSI=as.character())
  
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

    if (unpairedReplicates != 0) {
      norCandidate <- intraReplicate[[replicate]]$deltaPSI <= 0.05
      tumCandidate <- intraReplicate[[replicate]]$deltaPSI <= 0.05
    }
    
    Switch <- maxDeltaPsi - minDeltaPsi
    
    if (Switch < 0.2){
      next
    }
    
    replicateCandidates[[aCandidate]] <- data.frame(Genename=intraReplicate[[replicate]]$Genename[thisGeneData & norCandidate], 
                                                    Gene=aCandidate,
                                                    maxdPSI=intraReplicate[[replicate]]$Transcript[thisGeneData & norCandidate], 
                                                    mindPSI=intraReplicate[[replicate]]$Transcript[thisGeneData & tumCandidate],
                                                    Switch=Switch)
  }
  
  candidates[[replicate]] <- do.call('rbind', replicateCandidates)
  
  candidateMask[["N"]][[replicate]] <- interReplicate[["T"]]$Transcript %in% candidates[[replicate]]$maxdPSI
  candidateMask[["T"]][[replicate]] <- interReplicate[["T"]]$Transcript %in% candidates[[replicate]]$mindPSI
  
  #Expressed genes: transcripts whose expression are above the expression threshold in any of the conditions
  replicateExpressed[[replicate]] <- data.frame(Transcript=intraReplicate[[replicate]]$Transcript[norExpression | tumExpression],
                                                Genename=intraReplicate[[replicate]]$Genename[norExpression | tumExpression])
  
  cat(":", nrow(candidates[[replicate]]), "candidates found.\n")
  
}

allExpressedTranscripts <- unique(do.call('rbind', replicateExpressed))

candidateMask[["N"]] <- do.call('cbind', candidateMask[["N"]])
candidateMask[["T"]] <- do.call('cbind', candidateMask[["T"]])

candidateList <- do.call("rbind", candidates)
candidateList <- ddply(candidateList,.(Genename,Gene,maxdPSI,mindPSI), summarise, Replicated=length(Genename))
candidateList <- candidateList[with(candidateList, order(-Replicated)), ]

write.table(candidateList, file=paste0(out, "candidateList.tsv"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(allExpressedTranscripts, paste0(out, "expressedGenes.lst"), sep="\t", row.names=F, col.names=F, quote=F)

save(isoformExpression, intraReplicate, interReplicate, candidates, candidateList, candidateMask, inputData, allExpressedTranscripts, wd, out, file=paste0(out, "RWorkspaces/2_GetCandidates.RData"))

#Plot heatmap
suppressMessages(library(gplots)) #Avoid the annoying message
library(RColorBrewer)

top <- length(candidateList$Genename[candidateList$Replicated >= inputData[["Replicates"]] * 0.2])

topCandidates <- ddply(candidateList,.(Genename), summarise, Replicated=sum(Replicated))
topCandidates <- topCandidates[with(topCandidates, order(-Replicated)), ]
topCandidates <- head(topCandidates, n=top)

fig <- data.frame(matrix(nrow=length(topCandidates$Genename), ncol=inputData[["Replicates"]]))
rownames(fig) <- topCandidates$Genename
colnames(fig) <- seq(1,inputData[["Replicates"]])

for (replicate in seq(1,inputData[["Replicates"]])){
  for (gene in topCandidates$Genename){
    if (gene %in% candidates[[replicate]]$Genename) {
      fig[gene, replicate] <- head(candidates[[replicate]]$Switch[candidates[[replicate]]$Genename == gene],1)
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