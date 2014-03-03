#!/soft/R/R-3.0.0/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
tag <- args[1]
load(args[2])

allGenes <- as.character()
for (replicate in seq(1,inputData[["Replicates"]])){
  expressedGenes <- log(as.numeric(intraReplicate[[replicate]]$TPM_N)) > -1
  allGenes <- unique(c(allGenes, intraReplicate[[replicate]]$Gene[expressedGenes]))
}
allEntrezIds <- as.character(lapply(strsplit(allGenes, "|", fixed=T), function(x)x[2]))
allEntrezIds <- allEntrezIds[!is.na(allEntrezIds)]

minReplicates <- candidateList$Replicated >= inputData[["Replicates"]]/5
candidateGenes <- as.character(candidateList$Gene[minReplicates])
candEntrezIds <- as.character(lapply(strsplit(candidateGenes, "|", fixed=T), function(x)x[2]))
candEntrezIds <- candEntrezIds[!is.na(candEntrezIds)]

write(allEntrezIds, paste0(tag, "_expressedGenes.lst"), sep="\n")
write(candEntrezIds, paste0(tag, "_candidateGenes.lst"), sep="\n")
