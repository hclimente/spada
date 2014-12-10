#!/soft/R/R-3.0.0/bin/Rscript

getRobustZscore <- function(x){

  tx     <- x[1]
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

  return(data.frame(Transcript=tx,p=pvalue))
}

args <- commandArgs(trailingOnly = TRUE)
patient <- args[2]
load(paste0(args[1],"RWorkspaces/0_InitialEnvironment.RData"))
load(paste0(args[1],"RWorkspaces/cohortInfo.RData"))
load(paste0(args[1],"RWorkspaces/",patient,".RData"))

data <- merge(patientInfo[,c("Transcript", "deltaPSI")],cohortInfo[,c("Transcript","Median_dPSI","MAD_dPSI","MeAD_dPSI")])

pvals <- apply(data,1,getRobustZscore)
pvals <- do.call(rbind,pvals)
patientInfo <- merge(patientInfo,pvals,by="Transcript")

patientInfo$padj_up <- p.adjust(1-patientInfo$p, method="fdr")
patientInfo$padj_dw <- p.adjust(patientInfo$p, method="fdr")
significant <- patientInfo$padj_up < 0.05 | patientInfo$padj_dw < 0.05 
significant[is.na(significant)] <- FALSE

norExpression <- log(patientInfo$TPM_N) > inputData$minExpression
tumExpression <- log(patientInfo$TPM_T) > inputData$minExpression

candidates <- vector('list', length(unique(patientInfo$Gene[significant]))) 

for (gene in unique(patientInfo$Gene[significant])){ 
  
  thisGeneMask <- patientInfo$Gene == gene
      
  norTranscriptMask <- patientInfo$deltaPSI < 0 & significant & norExpression & thisGeneMask
  tumTranscriptMask <- patientInfo$deltaPSI > 0 & significant & tumExpression & thisGeneMask
  
  if(sum(norTranscriptMask,na.rm=T) == 0 || sum(tumTranscriptMask,na.rm=T) == 0){
    next
  }
  
  norTranscripts <- patientInfo[norTranscriptMask,]
  tumTranscripts <- patientInfo[tumTranscriptMask,]
  
  #Calculate max difference between the transcripts and choose the candidates
  #     deltaPSI = PSI_alt - PSI_ref
  # MaxDeltaPsi: predominant transcript in the tumor sample
  # MinDeltaPsi: predominant transcript in the normal sample
  maxP <- norTranscripts$padj_dw == min(norTranscripts$padj_dw, na.rm=T)
  if (sum(maxP) > 1){
    highExpression <- norTranscripts$TPM_N == max(norTranscripts$TPM_N, na.rm=T)
    maxP <- maxP & highExpression
  } 
  
  minP <- tumTranscripts$padj_up == min(tumTranscripts$padj_up, na.rm=T)
  if (sum(minP) > 1){
    highExpression <- tumTranscripts$TPM_T == max(tumTranscripts$TPM_T, na.rm=T)
    minP <- minP & highExpression
  }

  norCandidate <- norTranscripts$Transcript[maxP]
  tumCandidate <- tumTranscripts$Transcript[minP]
  
  if (sum(!is.na(norCandidate)) != 1 || sum(!is.na(tumCandidate)) != 1){
    next
  }

  Switch <- norTranscripts$deltaPSI[norTranscripts$Transcript==norCandidate] - tumTranscripts$deltaPSI[tumTranscripts$Transcript==tumCandidate]

  candidates[[gene]] <- data.frame(Gene=gene,Normal_isoform=norCandidate, 
                                   Tumor_isoform=tumCandidate,Switch=Switch)
}

patientCandidates <- do.call('rbind', candidates)

#Expressed genes: transcripts whose expression are above the expression threshold in any of the conditions
patientExpressed <- data.frame(Transcript=patientInfo$Transcript[norExpression | tumExpression],
                        Gene=patientInfo$Gene[norExpression | tumExpression])

save(patientCandidates,patientExpressed,file=paste0(args[1],"/RWorkspaces/",patient,"_candidates_expressed.RData"))
save(patientInfo,file=paste0(args[1],"/RWorkspaces/",patient,"_more.RData"))