#!/soft/R/R-3.0.0/bin/Rscript

suppressMessages(library(logging))

args <- commandArgs(trailingOnly = TRUE)
load(paste0(args[1], "RWorkspaces/0_InitialEnvironment.RData"))
inputPath <- args[2]

logger <- getLogger(name="exploreData", level=10) #Level debug

addHandler(writeToConsole, logger="explore_data", level='INFO')
addHandler(writeToFile, logger="explore_data", file=paste0(out, "rSmartAS.log"), level='DEBUG')

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

simplePlot <- function(x, y, title, xLab, yLab, pngName){
  png(pngName, width=960, height=960)
  plot(x, y, main=title, xlab=xLab, ylab=yLab)
  graphics.off()
}

getThresholds <- function(x){
  psiCols  <- paste0("PSI_",inputData$Replicates)
  
  diffMatrix <- abs(outer(as.numeric(x),as.numeric(x),"-"))
  subtraction <- diffMatrix[lower.tri(diffMatrix, diag = FALSE)]
  fpr1 <- as.numeric( quantile(subtraction, 0.99, na.rm=T) )
  fpr5 <- as.numeric( quantile(subtraction, 0.95, na.rm=T) )
  
  return(list(fpr1=fpr1,fpr5=fpr5))
}

calculateMAD <- function(x){

  result <- apply(x, 1, mad, na.rm=T)
  
  mad0_mask <- result == 0
  mad0_mask[is.na(mad0_mask)] <- TRUE
  mad0 <- x[mad0_mask, ]
  
  result[mad0_mask] <- apply(mad0, 1, function(x) mad(x,center = mean(as.numeric(x),na.rm=T), constant = 1.253314, na.rm = TRUE))
  
  return(result)

}

makeInterReplicateComparisons <- function(condition, samples){

  interRepComp <- NULL

  loginfo(paste("Summarizing data from",condition,as.character(length(samples)),"samples."),logger="explore_data")
  pairedPatientSumPB <- txtProgressBar(min=1, max=length(samples), initial=1, style=3)
  counter <- 1

  delete <- c(paste(c("TPM","tTPM", "PSI"),"T",sep="_"), "deltaPSI", "la_tTPM")

  for (replicate in samples){
    if(is.null(interRepComp)){
      interRepComp <- intraReplicate[[replicate]]    
    } else {
      interRepComp <- merge(interRepComp, intraReplicate[[replicate]],by=c("Gene", "Transcript"), 
                            all=TRUE,suffixes=c("",paste0("_", as.character(replicate))))
    }

    interRepComp <- interRepComp[,!(colnames(interRepComp) %in% delete), drop=FALSE]

    simplePlot(intraReplicate[[replicate]]$la_tTPM, intraReplicate[[replicate]]$deltaPSI, replicate, "0.5·(log(sum tTPM_N) + log(sum tTPM_T) )", 
               "deltaPSI", paste0(out, "DataExploration/latTPM_PSI_intrarreplicate",replicate,"_", condition,".png"))

    setTxtProgressBar(pairedPatientSumPB, counter)
    counter <- counter + 1

  }
  close(pairedPatientSumPB)
  return(interRepComp)
}

loginfo("Importing data from %d paired samples.",length(inputData$Replicates),logger="explore_data")
explorationPB <- txtProgressBar(min=1, max=length(inputData$Replicates), initial=1, style=3)
counter <- 1

for (replicate in inputData$Replicates){
  for (sample in inputData$Conditions){

    tag <- paste0(replicate, sample)
    inputFile=paste0(inputPath, replicate, "_", sample, ".tsv")
          
    #Read expression table
    isoformExpression[[tag]] <- read.table(inputFile, header=F, sep="\t", stringsAsFactors=F)
    colnames(isoformExpression[[tag]]) <- c("Gene", "Transcript","TPM")
      
    #Calculate the total expression of the genes and the PSI for each transcript
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
    
  #Plot deltaPSI distribution
  printLogFreqHist(intraReplicate[[replicate]]$deltaPSI, "deltaPSI", paste0("deltaPSI_", replicate))
  
  #Plot PSI distributions
  printLogFreqHist(intraReplicate[[replicate]]$PSI_N, "PSI_N", paste0("PSI_N_",replicate))
  printLogFreqHist(intraReplicate[[replicate]]$PSI_T, "PSI_T", paste0("PSI_T_",replicate))
  
  #Plot expression distribution
  printTPMHist(intraReplicate[[replicate]]$TPM_N, "log10(TPM_N+0.0001)", paste0("TPM_N_",replicate))
  printTPMHist(intraReplicate[[replicate]]$TPM_T, "log10(TPM_T+0.0001)", paste0("TPM_T_",replicate))
  
  setTxtProgressBar(explorationPB, counter)
  counter <- counter + 1

}

close(explorationPB)

tpmCols  <- paste0("TPM_",inputData$Replicates)
tTpmCols <- paste0("tTPM_",inputData$Replicates)
psiCols  <- paste0("PSI_",inputData$Replicates)

all <- c(rbind(tpmCols, tTpmCols, psiCols))

for (condition in inputData$Conditions){

  interReplicate[[condition]] <- makeInterReplicateComparisons(condition,inputData$Replicates)

  colnames(interReplicate[[condition]]) <- c("Gene", "Transcript", all)

  interReplicate[[condition]]$Median_PSI <- apply(interReplicate[[condition]][,psiCols], 1, median, na.rm=T)
  interReplicate[[condition]]$MAD_PSI <- calculateMAD(interReplicate[[condition]][,psiCols])
  interReplicate[[condition]]$Median_TPM <- apply(interReplicate[[condition]][,tpmCols], 1, median, na.rm=T)
  interReplicate[[condition]]$MAD_TPM <- calculateMAD(interReplicate[[condition]][,tpmCols])
  interReplicate[[condition]]$Median_tTPM <- apply(interReplicate[[condition]][,tTpmCols], 1, median, na.rm=T)
  interReplicate[[condition]]$MAD_tTPM <- calculateMAD(interReplicate[[condition]][,tTpmCols])

}

#Input unpaired samples, if any
if (length(inputData$unpairedReplicates) > 0){

  loginfo("Importing data from %d unpaired samples.",length(inputData$unpairedReplicates),logger="explore_data")
  unpairedPatientSumPB <- txtProgressBar(min=1, max=length(inputData$unpairedReplicates), initial=1, style=3)
  counter <- 1

  for (replicate in inputData$unpairedReplicates){
    inputFile=paste0(inputPath, replicate, "_T.tsv")
    
    intraReplicate[[replicate]] <- read.table(inputFile, header=F, sep="\t", stringsAsFactors=F)
    colnames(intraReplicate[[replicate]]) <- c("Gene","Transcript","TPM_T")
    vtTPM <- aggregate(TPM_T ~ Gene, data = intraReplicate[[replicate]], FUN = "sum")
    colnames(vtTPM) <- c("Gene", "tTPM_T")
    intraReplicate[[replicate]] <- merge(intraReplicate[[replicate]], vtTPM,by="Gene")
    intraReplicate[[replicate]] <- transform(intraReplicate[[replicate]], PSI_T = TPM_T / tTPM_T)

    intraReplicate[[replicate]] <- merge(intraReplicate[[replicate]], interReplicate[["N"]][,c("Gene","Transcript","Median_PSI","Median_TPM")], by=c("Gene","Transcript") )
    names(intraReplicate[[replicate]])[names(intraReplicate[[replicate]])=="Median_PSI"] <- "PSI_N"
    names(intraReplicate[[replicate]])[names(intraReplicate[[replicate]])=="Median_TPM"] <- "TPM_N"
    #Calculare a deltaPSI, the difference txtProgressBaretween the PSI in the tumor and the median PSI
    #for normal transcripts.
    intraReplicate[[replicate]]$deltaPSI <- intraReplicate[[replicate]]$PSI_N - intraReplicate[[replicate]]$PSI_T

    setTxtProgressBar(unpairedPatientSumPB, counter)
    counter <- counter + 1
  }

  close(unpairedPatientSumPB)

  interReplicate[["T"]] <- makeInterReplicateComparisons("T",inputData$unpairedReplicates)
}

loginfo("Summarizing data from %d transcripts.",length(interReplicate[["N"]]$Transcript),logger="explore_data")
fprThresholds <- apply(interReplicate[["N"]][,psiCols],1,getThresholds)
dfFprThresholds <- do.call(rbind, fprThresholds)

interReplicate[["N"]]$FPR_1 <- as.numeric(dfFprThresholds[,1])
interReplicate[["N"]]$FPR_5 <- as.numeric(dfFprThresholds[,2])

simplePlot(interReplicate[["N"]]$Median_PSI, interReplicate[["N"]]$MAD, "InterReplicate_N", "Median PSI", "MAD", paste0(out, "DataExploration/Interreplicate_mean_MAD_N.png"))
simplePlot(interReplicate[["N"]]$Median_PSI, interReplicate[["N"]]$FPR_1, "InterReplicate_N", "Median PSI", "FPR 1%", paste0(out, "DataExploration/Interreplicate_medianPSI_FRP1_N.png"))
simplePlot(interReplicate[["N"]]$Median_PSI, interReplicate[["N"]]$FPR_5, "InterReplicate_N", "Median PSI", "FPR 5%", paste0(out, "DataExploration/Interreplicate_medianPSI_FRP5_N.png"))

save(isoformExpression, intraReplicate, interReplicate, inputData, wd, out, file=paste0(out, "RWorkspaces/1_ExploreData.RData"))