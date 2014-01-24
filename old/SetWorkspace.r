#!/soft/R/R-3.0.0/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]

inputData <- list()
inputData[["Conditions"]] <- c("10", "7")
inputData[["Compartments"]] <- c("C")
inputData[["Replicates"]] <- c("1", "2")
inputData[["K-mer"]] <- c("30")

save.image("SmartAS.RData")