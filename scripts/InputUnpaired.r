#!/soft/R/R-3.0.0/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
load(paste0("Results/", args[1], "/RWorkspaces/1_ExploreData.RData"))
unpairedReplicates <- as.numeric(args[2])
inputPath <- args[3]

for (replicate in seq(1, unpairedReplicates)){
  cat("\t* Exploring tumor replicate",replicate, "/", inputData[["Replicates"]], "\t")

  inputFile=paste0(inputPath, replicate, "_T.tsv")
  outputFile=paste0(out, replicate, "_T.tsv")
    
  isoformExpression <- read.table(inputFile, header=F, sep="\t", stringsAsFactors=F)
  colnames(isoformExpression) <- c("Gene", "Transcript","Genename","TPM")
  vtTPM <- aggregate(TPM ~ Gene, data = isoformExpression, FUN = "sum")
  colnames(vtTPM) <- c("Gene", "tTPM")
  isoformExpression <- merge(isoformExpression, vtTPM)
  isoformExpression <- transform(isoformExpression, PSI = TPM / tTPM)
    
  if(!exists("interReplicate_T")){
    interReplicate_T <- isoformExpression
    interReplicate_T <- interReplicate_T[,!(colnames(interReplicate_T) %in% c("TPM", "tTPM")), drop=FALSE]
    
    columns <- c("Gene", "Transcript", "Genename", paste0("PSI_", replicate_c))
    
  } else {
    interReplicate_T <- merge(interReplicate_T, intraReplicate[[replicate]], by=c("Gene", "Transcript", "Genename"), suffixes=c("",paste0("_", replicate_c)), all=T)
    interReplicate_T <- interReplicate_T[,!(colnames(interReplicate_T) %in% delete), drop=FALSE]
    
    columns <- c(columns, paste0("PSI_", replicate_c))
  }
}

interReplicate_T$MedianPSI <- apply(interReplicate_T[,psiCols], 1, median)
interReplicate_T$MAD <- apply(interReplicate_T[,psiCols], 1, mad)

save(interReplicate_T, file=paste0(out, "RWorkspaces/1_Tumor_Intereplicate.RData"))