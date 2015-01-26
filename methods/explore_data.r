#!/soft/R/R-3.0.0/bin/Rscript

simplePlot <- function(x, y, title, xLab, yLab, pngName){
  png(pngName, width=960, height=960)
  plot(x, y, main=title, xlab=xLab, ylab=yLab)
  graphics.off()
}

getPseudodeltaPSIs <- function(x){
  
  diffMatrix <- abs(outer(as.numeric(x),as.numeric(x),"-"))
  subtraction <- diffMatrix[lower.tri(diffMatrix, diag = FALSE)]
  medPSI <- median(subtraction)
  madPSI <- mad(subtraction,na.rm=T)
  meadPSI <- mad(subtraction,center=mean(as.numeric(subtraction),na.rm=T))

  return(data.frame(median=medPSI,mad=madPSI,mead=meadPSI))
}

makeInterReplicateComparisons <- function(condition,samples,intraReplicate){

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

suppressMessages(library(logging))

args <- commandArgs(trailingOnly = TRUE)
load(paste0(args[1], "RWorkspaces/0_InitialEnvironment.RData"))
inputPath <- args[2]

allPatients <- c(inputData$Replicates,inputData$unpairedReplicates)

logger <- getLogger(name="exploreData", level=10) #Level debug

addHandler(writeToConsole, logger="explore_data", level='INFO')
addHandler(writeToFile, logger="explore_data", file=paste0(out, "rSmartAS.log"), level='DEBUG')

intraReplicate <- list()
interReplicate <- list(N=NULL, T=NULL)
isoformExpression <- list()

#Input paired samples
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
  intraReplicate[[replicate]]$deltaPSI <- intraReplicate[[replicate]]$PSI_T - intraReplicate[[replicate]]$PSI_N
  intraReplicate[[replicate]]$la_tTPM <- 0.5 * (log(intraReplicate[[replicate]]$tTPM_N) + log(intraReplicate[[replicate]]$tTPM_T))
   
  setTxtProgressBar(explorationPB, counter)
  counter <- counter + 1

}

close(explorationPB)

interReplicate[["N"]] <- makeInterReplicateComparisons("N",inputData$Replicates,intraReplicate)

tpmCols  <- paste0("TPM_",inputData$Replicates)
tTpmCols <- paste0("tTPM_",inputData$Replicates)
psiCols  <- paste0("PSI_",inputData$Replicates)

all <- c(rbind(tpmCols, tTpmCols, psiCols))

colnames(interReplicate[["N"]]) <- c("Gene", "Transcript", all)

interReplicate[["N"]]$Median_PSI <- apply(interReplicate[["N"]][,psiCols], 1, median, na.rm=T)
interReplicate[["N"]]$MAD_PSI <- apply(interReplicate[["N"]][,psiCols], 1, mad, na.rm=T)
interReplicate[["N"]]$MeAD_PSI <- apply(interReplicate[["N"]][,psiCols], 1, function(x) mad(x,center=mean(as.numeric(x),na.rm=T), na.rm=TRUE))
interReplicate[["N"]]$Median_TPM <- apply(interReplicate[["N"]][,tpmCols], 1, median, na.rm=T)
interReplicate[["N"]]$MAD_TPM <- apply(interReplicate[["N"]][,tpmCols], 1, mad, na.rm=T)
interReplicate[["N"]]$MeAD_TPM <- apply(interReplicate[["N"]][,tpmCols], 1, function(x) mad(x,center=mean(as.numeric(x),na.rm=T), na.rm=TRUE))
interReplicate[["N"]]$Median_tTPM <- apply(interReplicate[["N"]][,tTpmCols], 1, median, na.rm=T)
interReplicate[["N"]]$MAD_tTPM <- apply(interReplicate[["N"]][,tTpmCols], 1, mad, na.rm=T)
interReplicate[["N"]]$MeAD_tTPM <- apply(interReplicate[["N"]][,tTpmCols], 1, function(x) mad(x,center=mean(as.numeric(x),na.rm=T), na.rm=TRUE))

dPSI <- apply(interReplicate[["N"]][,psiCols], 1, getPseudodeltaPSIs)
dPSIThresholds <- do.call(rbind, dPSI)

interReplicate[["N"]]$Median_dPSI <- dPSIThresholds$median
interReplicate[["N"]]$MAD_dPSI <- dPSIThresholds$mad
interReplicate[["N"]]$MeAD_dPSI <- dPSIThresholds$mead

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
    intraReplicate[[replicate]]$deltaPSI <- intraReplicate[[replicate]]$PSI_T - intraReplicate[[replicate]]$PSI_N

    #Estimate of la_TPM in N and la_TPM
    vtTPM <- aggregate(TPM_N ~ Gene, data=intraReplicate[[replicate]], FUN = "sum")
    colnames(vtTPM) <- c("Gene", "tTPM_N")
    intraReplicate[[replicate]] <- merge(intraReplicate[[replicate]],vtTPM,by="Gene")
    
    intraReplicate[[replicate]]$la_tTPM <- 0.5 * (log(intraReplicate[[replicate]]$tTPM_N) + log(intraReplicate[[replicate]]$tTPM_T))

    setTxtProgressBar(unpairedPatientSumPB, counter)
    counter <- counter + 1
  }

  close(unpairedPatientSumPB)

}

#Make interreplicate comparisons and get data
interReplicate[["T"]] <- makeInterReplicateComparisons("T",allPatients,intraReplicate)

tpmCols  <- paste0("TPM_",allPatients)
tTpmCols <- paste0("tTPM_",allPatients)
psiCols  <- paste0("PSI_",allPatients)

all <- c(rbind(tpmCols, tTpmCols, psiCols))

colnames(interReplicate[["T"]]) <- c("Gene", "Transcript", all)

interReplicate[["T"]]$Median_PSI <- apply(interReplicate[["T"]][,psiCols], 1, median, na.rm=T)
interReplicate[["T"]]$MAD_PSI <- apply(interReplicate[["T"]][,psiCols], 1, mad, na.rm=T)
interReplicate[["T"]]$MeAD_PSI <- apply(interReplicate[["T"]][,psiCols], 1, function(x) mad(x,center=mean(as.numeric(x),na.rm=T), na.rm=TRUE))
interReplicate[["T"]]$Median_TPM <- apply(interReplicate[["T"]][,tpmCols], 1, median, na.rm=T)
interReplicate[["T"]]$MAD_TPM <- apply(interReplicate[["T"]][,tpmCols], 1, mad, na.rm=T)
interReplicate[["T"]]$MeAD_TPM <- apply(interReplicate[["T"]][,tpmCols], 1, function(x) mad(x,center=mean(as.numeric(x),na.rm=T), na.rm=TRUE))
interReplicate[["T"]]$Median_tTPM <- apply(interReplicate[["T"]][,tTpmCols], 1, median, na.rm=T)
interReplicate[["T"]]$MAD_tTPM <- apply(interReplicate[["T"]][,tTpmCols], 1, mad, na.rm=T)
interReplicate[["T"]]$MeAD_tTPM <- apply(interReplicate[["T"]][,tTpmCols], 1, function(x) mad(x,center=mean(as.numeric(x),na.rm=T), na.rm=TRUE))

cols <- c("Gene","Transcript","Median_PSI","MAD_PSI","MeAD_PSI","Median_TPM","MAD_TPM","MeAD_TPM","Median_tTPM","MAD_tTPM","MeAD_tTPM")
write.table(interReplicate[["N"]][,cols], file=paste0(out, "expression_normal.tsv"), sep="\t", row.names=F, quote=F)
write.table(interReplicate[["T"]][,cols], file=paste0(out, "expression_tumor.tsv"), sep="\t", row.names=F, quote=F)

save(isoformExpression, intraReplicate, interReplicate, inputData, wd, out, file=paste0(out, "RWorkspaces/1_ExploreData.RData"))

cohortInfo <- interReplicate[["N"]]
save(cohortInfo,file=paste0(out,"RWorkspaces/cohortInfo.RData"))
for ( patient in names(intraReplicate)){
  patientInfo <- intraReplicate[[patient]]
  save(patientInfo,file=paste0(out,"RWorkspaces/",patient,".RData"))
}

### Plot data
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

for (replicate in inputData$Replicates){
    
  #Plot deltaPSI distribution
  printLogFreqHist(intraReplicate[[replicate]]$deltaPSI, "deltaPSI", paste0("deltaPSI_", replicate))
  
  #Plot PSI distributions
  printLogFreqHist(intraReplicate[[replicate]]$PSI_N, "PSI_N", paste0("PSI_N_",replicate))
  printLogFreqHist(intraReplicate[[replicate]]$PSI_T, "PSI_T", paste0("PSI_T_",replicate))
  
  #Plot expression distribution
  printTPMHist(intraReplicate[[replicate]]$TPM_N, "log10(TPM_N+0.0001)", paste0("TPM_N_",replicate))
  printTPMHist(intraReplicate[[replicate]]$TPM_T, "log10(TPM_T+0.0001)", paste0("TPM_T_",replicate))

}

simplePlot(interReplicate[["N"]]$Median_PSI, interReplicate[["N"]]$MAD, "InterReplicate_N", "Median PSI", "MAD", paste0(out, "DataExploration/Interreplicate_mean_MAD_N.png"))
simplePlot(interReplicate[["N"]]$Median_PSI, interReplicate[["N"]]$FPR_1, "InterReplicate_N", "Median PSI", "FPR 1%", paste0(out, "DataExploration/Interreplicate_medianPSI_FRP1_N.png"))
simplePlot(interReplicate[["N"]]$Median_PSI, interReplicate[["N"]]$FPR_5, "InterReplicate_N", "Median PSI", "FPR 5%", paste0(out, "DataExploration/Interreplicate_medianPSI_FRP5_N.png"))