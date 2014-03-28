#!/soft/R/R-3.0.0/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
load(paste0("Results/", args[1], "/RWorkspaces/0_InitialEnvironment.RData"))
inputPath <- args[2]

intraReplicate <- list()
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

for (replicate in seq(1, inputData[["Replicates"]])){
  cat("\t* Exploring replicate",replicate, "/", inputData[["Replicates"]], "\t")
  for (sample in inputData[["Conditions"]]){

    tag <- paste0(replicate, sample)
    cat(" ",sample)
    inputFile=paste0(inputPath, replicate, "_", sample, ".tsv")
    outputFile=paste0(out, replicate, "_", sample, ".tsv")
      
    #Read Sailfish table
    isoformExpression[[tag]] <- read.table(inputFile, header=F, sep="\t", stringsAsFactors=F)
    colnames(isoformExpression[[tag]]) <- c("Gene", "Transcript","Genename","TPM")
      
    #Calculate the PSI for each transcript and the total expression of the gene
    vtTPM <- aggregate(TPM ~ Gene, data = isoformExpression[[tag]], FUN = "sum")
    colnames(vtTPM) <- c("Gene", "tTPM")
    isoformExpression[[tag]] <- merge(isoformExpression[[tag]], vtTPM)
    isoformExpression[[tag]] <- transform(isoformExpression[[tag]], PSI = TPM / tTPM)
      
    #write.table(isoformExpression[[tag]], file=outputFile, sep="\t", row.names=F)
  }

  cat("\n")
   
  nTag <- paste0(replicate, "N")
  tTag <- paste0(replicate, "T")
  
  intraReplicate[[replicate]] <- merge(isoformExpression[[nTag]], isoformExpression[[tTag]], by=c("Gene", "Transcript", "Genename"), suffixes=c("_N","_T"), all=T)
  intraReplicate[[replicate]]$deltaPSI <- intraReplicate[[replicate]]$PSI_N - intraReplicate[[replicate]]$PSI_T
  intraReplicate[[replicate]]$la_tTPM <- 0.5 * (log(intraReplicate[[replicate]]$tTPM_N) + log(intraReplicate[[replicate]]$tTPM_T))
    
  #Plot stuff
  printLogFreqHist(intraReplicate[[replicate]]$deltaPSI, "deltaPSI", paste0("deltaPSI_", replicate))
  printTPMHist(intraReplicate[[replicate]]$TPM_N, "log10(TPM_N+0.0001)", paste0("TPM_N_",replicate))
  printLogFreqHist(intraReplicate[[replicate]]$PSI_N, "PSI_N", paste0("PSI_N_",replicate))
  printTPMHist(intraReplicate[[replicate]]$TPM_T, "log10(TPM_T+0.0001)", paste0("TPM_T_",replicate))
  printLogFreqHist(intraReplicate[[replicate]]$PSI_N, "PSI_T", paste0("PSI_T_",replicate))
    
  #write.table(intraReplicate[[replicate]], file=paste0(out, "IntraReplicate_",replicate,".tsv"), sep="\t", row.names=F)

}

psiCols <- as.character()
delete <- c("TPM_T","tTPM_T", "PSI_T", "deltaPSI", "la_tTPM", "TPM_N", "tTPM_N")

for (replicate in seq(1,inputData[["Replicates"]])){

  replicate_c <- as.character(replicate)
      
  if(!exists("interReplicate_N")){
    interReplicate_N <- intraReplicate[[replicate]]
    interReplicate_N <- interReplicate_N[,!(colnames(interReplicate_N) %in% delete), drop=FALSE]
    
    columns <- c("Gene", "Transcript", "Genename", paste0("PSI_", replicate_c))
    
  } else {
    interReplicate_N <- merge(interReplicate_N, intraReplicate[[replicate]], by=c("Gene", "Transcript", "Genename"), suffixes=c("",paste0("_", replicate_c)), all=T)
    interReplicate_N <- interReplicate_N[,!(colnames(interReplicate_N) %in% delete), drop=FALSE]
    
    columns <- c(columns, paste0("PSI_", replicate_c))
  }
  colnames(interReplicate_N) <- columns
  psiCols <- c(psiCols, paste0("PSI_", replicate))

  simplePlot(intraReplicate[[replicate]]$la_tTPM, intraReplicate[[replicate]]$deltaPSI, replicate, "0.5·(log(sum tTPM_N) + log(sum tTPM_T) )", 
             "deltaPSI", paste0(out, "DataExploration/latTPM_PSI_intrarreplicate",replicate,".png"))

}

interReplicate_N$MedianPSI <- apply(interReplicate_N[,psiCols], 1, median, na.rm=T)
interReplicate_N$MAD <- apply(interReplicate_N[,psiCols], 1, mad, na.rm=T)

q01 <- as.numeric()
q05 <- as.numeric()
q95 <- as.numeric()
q99 <- as.numeric()
fpr <- as.numeric()

for (tx in interReplicate_N$Transcript){
  thisTranscript <- interReplicate_N[interReplicate_N$Transcript==tx,]
  subtraction <- as.numeric()
  for (i in seq(1, inputData[["Replicates"]] - 1) ){
    s <- thisTranscript[,paste0("PSI_",i)] - as.numeric(thisTranscript[,psiCols])
    subtraction <- c(subtraction,s)
  }
  
  q01 <- c(q01, as.numeric( quantile(subtraction, 0.01, na.rm=T) ) )
  q05 <- c(q05, as.numeric( quantile(subtraction, 0.05, na.rm=T) ) )
  q95 <- c(q95, as.numeric( quantile(subtraction, 0.95, na.rm=T) ) )
  q99 <- c(q99, as.numeric( quantile(subtraction, 0.99, na.rm=T) ) )

  if ( any(!is.na(subtraction)) ) {
    fpr <- c(fpr, ecdf(subtraction)(4 * thisTranscript$MAD) )
  } else {
    fpr <- c(fpr, NA)
  }

}

interReplicate_N$FPR_01 <- q01
interReplicate_N$FPR_05 <- q05
interReplicate_N$FPR_95 <- q95
interReplicate_N$FPR_99 <- q99
interReplicate_N$FPR_4MAD <- fpr

simplePlot(interReplicate_N$MedianPSI, interReplicate_N$MAD, "InterReplicate_N", "Median PSI", "MAD", paste0(out, "DataExploration/Interreplicate_mean_MAD_N.png"))
simplePlot(interReplicate_N$MedianPSI, interReplicate_N$FPR_01, "InterReplicate_N", "Median PSI", "FPR 1%", paste0(out, "DataExploration/Interreplicate_medianPSI_FRP1_N.png"))
simplePlot(interReplicate_N$MedianPSI, interReplicate_N$FPR_05, "InterReplicate_N", "Median PSI", "FPR 5%", paste0(out, "DataExploration/Interreplicate_medianPSI_FRP5_N.png"))
simplePlot(interReplicate_N$MedianPSI, interReplicate_N$FPR_95, "InterReplicate_N", "Median PSI", "FPR 95%", paste0(out, "DataExploration/Interreplicate_medianPSI_FRP95_N.png"))
simplePlot(interReplicate_N$MedianPSI, interReplicate_N$FPR_99, "InterReplicate_N", "Median PSI", "FPR 99%", paste0(out, "DataExploration/Interreplicate_medianPSI_FRP99_N.png"))

save(isoformExpression, intraReplicate, interReplicate_N, inputData, wd, out, file=paste0(out, "RWorkspaces/1_ExploreData.RData"))