#!/usr/bin/env Rscript

library(ggplot2)
library(edgeR)

tryCatch(source("~/smartas/pipeline/scripts/variablesAndFunctions.r"),error=function(e){})

args <- commandArgs(trailingOnly = TRUE)

switches.file <- args[1]
xpr.nt.file <- args[2]
xpr.t.file <- args[3]
cancer <- args[4]

# test
#switches.file <- "~/smartas/analyses/luad/candidateList_info.tsv"
#xpr.nt.file <- "/projects_rg/TCGA/pipeline/run11/luad_gene_read_paired-filtered.txt"
#xpr.t.file <- "/projects_rg/TCGA/pipeline/run11/luad_gene_read_tumor-filtered.txt"
#cancer <- "luad"
#####

# read files
switches <- read.table(switches.file,header=T)
switches$Symbol <- as.character(switches$Symbol)
switches$Normal_transcript <- as.character(switches$Normal_transcript)
switches$Tumor_transcript <- as.character(switches$Tumor_transcript)
switches$Patients_affected <- as.character(switches$Patients_affected)
xpr.nt <- read.table(xpr.nt.file, check.names=FALSE)
xpr.t <- read.table(xpr.t.file, check.names=FALSE)

xpr <- cbind(xpr.nt,xpr.t)
rm(xpr.nt,xpr.t)

# format files
y <- DGEList(counts=xpr)
y <- calcNormFactors(y)
xpr.norm <- cpm(y, normalized.lib.sizes=TRUE)

# get sample names
normal <- grep("^.{4}N$",colnames(xpr.norm), value=TRUE)
tumor <- grep("^.{4}T$",colnames(xpr.norm), value=TRUE)
tumor.paired <- gsub("N$","T",normal)
tumor.unpaired <- setdiff(tumor,tumor.paired)

# subset
## get patients
for (i in 1:nrow(switches)){
  x <- switches[i,]
  gene = paste(x[2],x[1],sep="|")
  niso = x[3]
  tiso = x[4]
  switch.patients = unlist(strsplit(as.character(x[18]),","))
  
  # get tpm
  switch.xpr.norm <- as.data.frame(xpr.norm[rownames(xpr.norm)==gene,])
  colnames(switch.xpr.norm) <- c("Expression")
  
  switch.xpr.norm$What <- "Normal"
  switch.xpr.norm$What[rownames(switch.xpr.norm) %in% tumor] <- "Tumor no switch"
  switch.xpr.norm$What[rownames(switch.xpr.norm) %in% switch.patients] <- "Tumor switch"
  
  # plot
  p <- ggplot(switch.xpr.norm,aes(x=What,y=log2(Expression))) + 
    geom_boxplot() +
    smartas_theme() +
    labs(x="Sample type", y="log2(gene cpm)")
  
  ggsave(paste(cancer,i,"reads",gene,niso,tiso,"png",sep="."),p,width = 12,height = 12)
}