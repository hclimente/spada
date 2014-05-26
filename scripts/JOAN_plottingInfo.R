#!/soft/R/R-3.0.0/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
gene <- args[1]
nTx <- args[2]
tTx <- args[3]

load(paste0("~/SmartAS/Results/TCGA/luad_mE-1.0/RWorkspaces/5_Interactions_", gene, "_", nTx, "_", tTx, "_1.RData"))
load("~/SmartAS/Results/TCGA/luad_mE-1.0/RWorkspaces/2_GetCandidates.RData")

finalDf <- InteraX[,c("Partner","Partner_gene","RC_N","RC_T","dRC","Annotation")]
expression_N <- interReplicate[["N"]][,c("Transcript", "Gene", "Median_TPM")]
expression_T <- interReplicate[["T"]][,c("Transcript", "Median_TPM")]

finalDf <- merge(finalDf, expression_N, by.x="Partner", by.y="Transcript")
colnames(finalDf) <- c("Partner","Partner_gene","RC_N","RC_T","dRC","Annotation","Gene_id","pMedian_TPM_N")
finalDf <- merge(finalDf, expression_T, by.x="Partner", by.y="Transcript")
colnames(finalDf) <- c("Partner","Partner_gene","RC_N","RC_T","dRC","Annotation","Gene_id","pMedian_TPM_N","pMedian_TPM_T")
finalDf$cMedian_TPM_N=expression_N$Median_TPM[expression_N$Transcript==nTx]
finalDf$cMedian_TPM_T=expression_N$Median_TPM[expression_N$Transcript==tTx]

write.table(finalDf, file=paste0("~/Desktop/plottingInfo_", gene, ".tsv"), quote=F, sep="\t", row.names=F)