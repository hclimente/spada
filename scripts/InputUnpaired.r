#!/soft/R/R-3.0.0/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
load(paste0("Results/", args[1], "/RWorkspaces/1_ExploreData.RData"))
unpairedReplicates <- as.numeric(args[2])
inputPath <- args[3]

intraReplicate <- list()

for (replicate in seq(1, unpairedReplicates)){
  cat("\t* Exploring tumor replicate",replicate, "/", unpairedReplicates, "\n")

  inputFile=paste0(inputPath, replicate, "_T.tsv")
  outputFile=paste0(out, replicate, "_T.tsv")

  intraReplicate[[replicate]] <- read.table(inputFile, header=F, sep="\t", stringsAsFactors=F)
  colnames(intraReplicate[[replicate]]) <- c("Gene", "Transcript","Genename","TPM_T")
  vtTPM <- aggregate(TPM_T ~ Gene, data = intraReplicate[[replicate]], FUN = "sum")
  colnames(vtTPM) <- c("Gene", "tTPM_T")
  intraReplicate[[replicate]] <- merge(intraReplicate[[replicate]], vtTPM)
  intraReplicate[[replicate]] <- transform(intraReplicate[[replicate]], PSI_T = TPM_T / tTPM)
  intraReplicate[[replicate]]$TPM_N = 9999
  if all(intraReplicate[[replicate]]$Transcript == interReplicate_N$Transcript){
    intraReplicate[[replicate]]$deltaPSI <- 0.6745 * (intraReplicate[[replicate]]$PSI - interReplicate_N$MedianPSI) / interReplicate_N$MAD
  } else {
    cat("Error\n")
  }

}

save(intraReplicate, file=paste0(out, "RWorkspaces/1_Tumor_Intereplicate.RData"))