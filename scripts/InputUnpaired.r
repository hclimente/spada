#!/soft/R/R-3.0.0/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

load(paste0("Results/", args[1], "/RWorkspaces/1_ExploreData.RData"))
out <- paste0("Results/", args[1], "_unpaired/")
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
  intraReplicate[[replicate]] <- transform(intraReplicate[[replicate]], PSI_T = TPM_T / tTPM_T)
  intraReplicate[[replicate]]$TPM_N = 9999

  tmp <- merge(intraReplicate[[replicate]], interReplicate[["T"]], by=c("Gene", "Genename", "Transcript") )
  tmp <- tmp[,! names(tmp) %in% c(paste0("PSI_", seq(1,inputData[["Replicates"]]) ), paste0("TPM_", seq(1,inputData[["Replicates"]]) ), paste0("tTPM_", seq(1,inputData[["Replicates"]]) ) ) ]
  #rzs <- sweep(tmp$PSI_T, 1, tmp$Median_PSI, "-")
  #rzs <- sweep(rzs, 1, tmp$MAD_PSI, "/")
  #tmp$deltaPSI <- rzs
  intraReplicate[[replicate]]$deltaPSI <- (tmp$PSI_T - tmp$Median_PSI ) / tmp$MAD_PSI
}

save(intraReplicate, out, file=paste0(out, "RWorkspaces/1_Tumor_Intereplicate.RData"))