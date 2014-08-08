#!/soft/R/R-3.0.0/bin/Rscript

suppressMessages(library(logging))

args <- commandArgs(trailingOnly = TRUE)
load(paste0(args[1], "RWorkspaces/0_InitialEnvironment.RData"))
inputPath <- args[2]

logger <- getLogger(name="exploreData", level=10) #Level debug

addHandler(writeToConsole, logger="exploreData", level='INFO')
addHandler(writeToFile, logger="exploreData", file=paste0(out, "rSmartAS.log"), level='DEBUG')

intraReplicate <- list()
interReplicate <- list(N=NULL, T=NULL)
isoformExpression <- list()

printTPMHist <- function(x, xLab, pngName){
  png(paste0(out, "DataExploration/", pngName, ".png"), width=960, height=960)
  histogram <- hist(log10(x + 0.0001), 10000)
  histogram$counts <- log10(histogram$counts)
  plot(histogram$mids, histogram$counts, type="h", main=pngName, xlab=xLab, ylab="log10(Frequency)")
  graphics.off()
}

printLogFreqHist <- function(x, xLab, pngName){
  png(paste0(out, "DataExploration/", pngName,".png"), width=960, height=960)
  histogram <- hist(x, 10000)
  histogram$counts <- log10(histogram$counts)
  plot(histogram$mids, histogram$counts, type="h", main=tag, xlab=xLab, ylab="log10(Frequency)")
  graphics.off()
}

plotCorrelations <- function(x, y, lab, pngName){
  xLab=paste0(lab, " Replicate 1")
  yLab=paste0(lab, " Replicate 2")
  png(paste0(out, "DataExploration/", pngName, ".png"), width=960, height=960)
  plot(x, y, xlab=xLab, ylab=yLab)
  graphics.off()
  cor(x, y, use="complete.obs")
}

simplePlot <- function(x, y, title, xLab, yLab, pngName){
  png(pngName, width=960, height=960)
  plot(x, y, main=title, xlab=xLab, ylab=yLab)
  graphics.off()
}

calculateMAD <- function(x){

  result <- apply(x, 1, mad, na.rm=T)
  
  mad0_mask <- result == 0
  mad0_mask[is.na(mad0_mask)] <- TRUE
  mad0 <- x[mad0_mask, ]
  
  result[mad0_mask] <- apply(mad0, 1, function(x) mad(x,center = mean(as.numeric(x),na.rm=T), constant = 1.253314, na.rm = TRUE))
  
  return(result)

}

cat("Importing data from", inputData[["Replicates"]], "patients.\n")
explorationPB <- txtProgressBar(min=1, max=inputData[["Replicates"]], initial=1, style=3)
counter <- 1

for (replicate in seq(1, inputData[["Replicates"]])){

  for (sample in inputData[["Conditions"]]){

    tag <- paste0(replicate, sample)
    inputFile=paste0(inputPath, replicate, "_", sample, ".tsv")
    outputFile=paste0(out, replicate, "_", sample, ".tsv")
      
    #Read Sailfish table
    isoformExpression[[tag]] <- read.table(inputFile, header=F, sep="\t", stringsAsFactors=F)
    colnames(isoformExpression[[tag]]) <- c("Gene", "Transcript","TPM")
      
    #Calculate the PSI for each transcript and the total expression of the gene
    vtTPM <- aggregate(TPM ~ Gene, data=isoformExpression[[tag]], FUN = "sum")
    colnames(vtTPM) <- c("Gene", "tTPM")
    isoformExpression[[tag]] <- merge(isoformExpression[[tag]], vtTPM)
    isoformExpression[[tag]] <- transform(isoformExpression[[tag]], PSI = TPM / tTPM)
  
  }
   
  nTag <- paste0(replicate, "N")
  tTag <- paste0(replicate, "T")
  
  intraReplicate[[replicate]] <- merge(isoformExpression[[nTag]], isoformExpression[[tTag]], by=c("Gene", "Transcript"), suffixes=c("_N","_T"), all=T)
  intraReplicate[[replicate]]$deltaPSI <- intraReplicate[[replicate]]$PSI_N - intraReplicate[[replicate]]$PSI_T
  intraReplicate[[replicate]]$la_tTPM <- 0.5 * (log(intraReplicate[[replicate]]$tTPM_N) + log(intraReplicate[[replicate]]$tTPM_T))
    
  #Plot stuff
  printLogFreqHist(intraReplicate[[replicate]]$deltaPSI, "deltaPSI", paste0("deltaPSI_", replicate))
  printTPMHist(intraReplicate[[replicate]]$TPM_N, "log10(TPM_N+0.0001)", paste0("TPM_N_",replicate))
  printLogFreqHist(intraReplicate[[replicate]]$PSI_N, "PSI_N", paste0("PSI_N_",replicate))
  printTPMHist(intraReplicate[[replicate]]$TPM_T, "log10(TPM_T+0.0001)", paste0("TPM_T_",replicate))
  printLogFreqHist(intraReplicate[[replicate]]$PSI_N, "PSI_T", paste0("PSI_T_",replicate))

  setTxtProgressBar(explorationPB, counter)
  counter <- counter + 1

}

close(explorationPB)

tpmCols <- paste0("TPM_", seq(1,inputData[["Replicates"]]) )
ttpmCols <- paste0("tTPM_", seq(1,inputData[["Replicates"]]) )
psiCols <- paste0("PSI_", seq(1,inputData[["Replicates"]]) )

all <- c(rbind(tpmCols, ttpmCols, psiCols))

delete <- list()
delete[["N"]] <- c("TPM_T","tTPM_T", "PSI_T", "deltaPSI", "la_tTPM")
delete[["T"]] <- c("TPM_N","tTPM_N", "PSI_N", "deltaPSI", "la_tTPM")

cat("Summarizing data from", inputData[["Replicates"]], "patients.\n")
patientSumPB <- txtProgressBar(min=1, max=inputData[["Replicates"]] * 2, initial=1, style=3)
counter <- 1

for (sample in inputData[["Conditions"]]){

  for (replicate in seq(1,inputData[["Replicates"]])){

    replicate_c <- as.character(replicate)
        
    if(is.null(interReplicate[[sample]])){
      interReplicate[[sample]] <- intraReplicate[[replicate]]    
    } else {
      interReplicate[[sample]] <- merge(interReplicate[[sample]], intraReplicate[[replicate]], by=c("Gene", "Transcript"), suffixes=c("",paste0("_", replicate_c)), all=T)
    }

    interReplicate[[sample]] <- interReplicate[[sample]][,!(colnames(interReplicate[[sample]]) %in% delete[[sample]]), drop=FALSE]

    simplePlot(intraReplicate[[replicate]]$la_tTPM, intraReplicate[[replicate]]$deltaPSI, replicate, "0.5·(log(sum tTPM_N) + log(sum tTPM_T) )", 
               "deltaPSI", paste0(out, "DataExploration/latTPM_PSI_intrarreplicate",replicate,"_", sample,".png"))

    setTxtProgressBar(patientSumPB, counter)
    counter <- counter + 1

  }

  colnames(interReplicate[[sample]]) <- c("Gene", "Transcript", all)

  interReplicate[[sample]]$Median_PSI <- apply(interReplicate[[sample]][,psiCols], 1, median, na.rm=T)
  interReplicate[[sample]]$MAD_PSI <- calculateMAD(interReplicate[[sample]][,psiCols])
  interReplicate[[sample]]$Median_TPM <- apply(interReplicate[[sample]][,tpmCols], 1, median, na.rm=T)
  interReplicate[[sample]]$MAD_TPM <- calculateMAD(interReplicate[[sample]][,tpmCols])
  interReplicate[[sample]]$Median_tTPM <- apply(interReplicate[[sample]][,ttpmCols], 1, median, na.rm=T)
  interReplicate[[sample]]$MAD_tTPM <- calculateMAD(interReplicate[[sample]][,ttpmCols])

}

close(patientSumPB)

fpr1 <- as.numeric()
fpr5 <- as.numeric()
fprMAD <- as.numeric()

cat("Summarizing data from", length(interReplicate[["N"]]$Transcript), "transcripts.\n")
sampleSumPB <- txtProgressBar(min = 1, max = length(interReplicate[["N"]]$Transcript), initial = 1, style=3)
counter <- 1

for (tx in interReplicate[["N"]]$Transcript){
  thisTranscript <- interReplicate[["N"]][interReplicate[["N"]]$Transcript==tx,psiCols]
  this4MAD <- 4 * interReplicate[["N"]]$MAD_PSI[interReplicate[["N"]]$Transcript==tx]
  
  diffMatrix <- abs(outer(as.numeric(thisTranscript),as.numeric(thisTranscript),"-"))
  subtraction <- diffMatrix[lower.tri(diffMatrix, diag = FALSE)]
  fpr1 <- c(fpr1, as.numeric( quantile(subtraction, 0.99, na.rm=T) ) )
  fpr5 <- c(fpr5, as.numeric( quantile(subtraction, 0.95, na.rm=T) ) )
  if ( any( !is.na(subtraction) ) ) {
    fprMAD <- c(fprMAD, 1 - ecdf(subtraction)(this4MAD) )
  } else {
    fprMAD <- c(fprMAD, NA)
  }

  setTxtProgressBar(sampleSumPB, counter)
  counter <- counter + 1
}

close(sampleSumPB)

interReplicate[["N"]]$FPR_1 <- fpr1
interReplicate[["N"]]$FPR_5 <- fpr5
interReplicate[["N"]]$FPR_4MAD <- fprMAD

simplePlot(interReplicate[["N"]]$Median_PSI, interReplicate[["N"]]$MAD, "InterReplicate_N", "Median PSI", "MAD", paste0(out, "DataExploration/Interreplicate_mean_MAD_N.png"))
simplePlot(interReplicate[["N"]]$Median_PSI, interReplicate[["N"]]$FPR_1, "InterReplicate_N", "Median PSI", "FPR 1%", paste0(out, "DataExploration/Interreplicate_medianPSI_FRP1_N.png"))
simplePlot(interReplicate[["N"]]$Median_PSI, interReplicate[["N"]]$FPR_5, "InterReplicate_N", "Median PSI", "FPR 5%", paste0(out, "DataExploration/Interreplicate_medianPSI_FRP5_N.png"))

save(isoformExpression, intraReplicate, interReplicate, inputData, wd, out, file=paste0(out, "RWorkspaces/1_ExploreData.RData"))