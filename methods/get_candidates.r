#!/soft/R/R-3.0.0/bin/Rscript

library(plyr)
suppressMessages( library(logging) )

args <- commandArgs(trailingOnly = TRUE)
minExpression <- as.numeric(args[1])
load(paste0(args[2], "RWorkspaces/1_ExploreData.RData"))

patientSet <- c(inputData$Replicates,inputData$unpairedReplicates)
numOfReplicates <- length(patientSet)

logger <- getLogger(name="get_candidates", level=10) #Level debug

addHandler(writeToConsole, logger="get_candidates", level='INFO')
addHandler(writeToFile, logger="get_candidates", file=paste0(out, "rSmartAS.log"), level='DEBUG')

getRobustZscore <- function(x){

  s      <- as.numeric(x[2])
  median <- as.numeric(x[3])
  mad    <- as.numeric(x[4])
  mead   <- as.numeric(x[5])

  if(is.na(median)){
    return(NA)
  } else if (!is.na(mad)){
    if(mad != 0){
      z <- (s - median)/(1.486*mad)
    } else {
      z <- (s - median)/(1.253314*mead)
    } 
  } else if (!is.na(mead) & mead != 0) {
    z <- (s - median)/(1.253314*mead)
  } else {
    return(NA)
  }
  pvalue = pnorm(z)

  return(pvalue)
}

candidates <- list()
replicateExpressed  <- list()

loginfo("Searching isoform switches in %d patients.", numOfReplicates, logger="get_candidates")

candidatesPB <- txtProgressBar(min=1, max=numOfReplicates, initial = 1, style=3)
counter <- 1

for (replicate in patientSet){
  
  #Filter by deltaPSI and expression, based on the FPR
  data <- merge(intraReplicate[[replicate]][,c("Transcript", "deltaPSI")],interReplicate[["N"]][,c("Transcript","Median_dPSI","MAD_dPSI","MeAD_dPSI")])

  intraReplicate[[replicate]]$p <- apply(data,1,getRobustZscore)
  intraReplicate[[replicate]]$padj <- p.adjust(intraReplicate[[replicate]]$p, method="fdr")
  significant <- intraReplicate[[replicate]]$padj < 0.05 | intraReplicate[[replicate]]$padj > 0.95 
  significant[is.na(significant)] <- FALSE
  
  norExpression <- log(intraReplicate[[replicate]]$TPM_N) > minExpression
  tumExpression <- log(intraReplicate[[replicate]]$TPM_T) > minExpression
  
  replicateCandidates <- vector('list', length(unique(intraReplicate[[replicate]]$Gene[significant]))) 
  
  for (aCandidate in unique(intraReplicate[[replicate]]$Gene[significant])){ 
    
    thisGeneMask <- intraReplicate[[replicate]]$Gene == aCandidate
        
    norTranscriptMask <- intraReplicate[[replicate]]$deltaPSI < 0 & significant & norExpression & thisGeneMask
    tumTranscriptMask <- intraReplicate[[replicate]]$deltaPSI > 0 & significant & tumExpression & thisGeneMask
    
    if(sum(norTranscriptMask,na.rm=T) == 0 || sum(tumTranscriptMask,na.rm=T) == 0){
      next
    }
    
    norTranscripts <- intraReplicate[[replicate]][norTranscriptMask,]
    tumTranscripts <- intraReplicate[[replicate]][tumTranscriptMask,]
    
    #Calculate max difference between the transcripts and choose the candidates
    #     deltaPSI = PSI_alt - PSI_ref
    # MaxDeltaPsiCond: predominant transcript in the normal sample
    # MinDeltaPsiCond: predominant transcript in the tumor sample
    maxP <- norTranscripts$padj == max(norTranscripts$padj, na.rm=T)
    if (sum(maxP) > 1){
      maxTPM <- norTranscripts$TPM_N == max(norTranscripts$TPM_N, na.rm=T)
      norCandidate <- norTranscripts$Transcript[maxP & maxTPM]  
    } else {
      norCandidate <- norTranscripts$Transcript[maxP]
    }
    
    minP <- tumTranscripts$padj == max(tumTranscripts$padj, na.rm=T)
    if (sum(minP) > 1){
      minTPM <- tumTranscripts$TPM_T == max(tumTranscripts$TPM_T, na.rm=T)
      tumCandidate <- tumTranscripts$Transcript[minP & minTPM]
    } else {
      tumCandidate <- tumTranscripts$Transcript[minP]
    }
    
    if (sum(!is.na(norCandidate)) == 0 || sum(!is.na(tumCandidate)) == 0){
      next
    }
  
    Switch <- norTranscripts$deltaPSI[norTranscripts$Transcript==norCandidate] - tumTranscripts$deltaPSI[tumTranscripts$Transcript==tumCandidate]

    replicateCandidates[[aCandidate]] <- data.frame(Gene=aCandidate,maxdPSI=norCandidate, 
                                                    mindPSI=tumCandidate,Switch=Switch)
  }
  
  candidates[[replicate]] <- do.call('rbind', replicateCandidates)
  
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

for (i in patientSet){candidates[[i]]$Origin = i}

candidateList <- do.call("rbind", candidates)
candidateList <- ddply(candidateList,.(Gene,maxdPSI,mindPSI), summarise, Replicated=length(Gene), Patients=paste(Origin, collapse = ",") )
candidateList <- candidateList[with(candidateList, order(-Replicated)), ]

write.table(candidateList, file=paste0(out, "candidateList.tsv"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(allExpressedTranscripts, paste0(out, "expressedGenes.lst"), sep="\t", row.names=F, col.names=F, quote=F)

cols <- c("Gene","Transcript","Median_PSI","MAD_PSI","MeAD_PSI","Median_TPM","MAD_TPM","MeAD_TPM","Median_tTPM","MAD_tTPM","MeAD_tTPM")
write.table(interReplicate[["N"]][,cols], file=paste0(out, "expression_normal.tsv"), sep="\t", row.names=F, quote=F)
write.table(interReplicate[["T"]][,cols], file=paste0(out, "expression_tumor.tsv"), sep="\t", row.names=F, quote=F)

save(isoformExpression, intraReplicate, interReplicate, candidates, candidateList, inputData, allExpressedTranscripts, wd, out, file=paste0(out, "RWorkspaces/2_GetCandidates.RData"))

#Plot heatmap
suppressMessages(library(gplots)) #Avoid the annoying message
library(RColorBrewer) 

library(RColorBrewer) 
top <- length(candidateList$Gene[candidateList$Replicated >= length(patientSet) * 0.2])

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

save(isoformExpression, intraReplicate, interReplicate, candidates, candidateList, inputData, allExpressedTranscripts, wd, out, fig, file=paste0(out, "RWorkspaces/2_GetCandidates.RData"))