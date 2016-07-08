#!/usr/bin/env Rscript

library(ggplot2)
library(grid)
library(gridExtra)

tryCatch(source("~/smartas/pipeline/scripts/variablesAndFunctions.r"),error=function(e){})

args <- commandArgs(trailingOnly = TRUE)

switches.file <- args[1]
psi.nt.file <- args[2]
psi.t.file <- args[3]
cancer <- args[4]

# test
# switches.file <- "~/smartas/analyses/brca/candidateList_info.tsv"
# psi.nt.file <- "/projects_rg/TCGA/pipeline/run11/brca_iso_psi_paired-filtered.txt"
# psi.t.file <- "/projects_rg/TCGA/pipeline/run11/brca_iso_psi_tumor-filtered.txt"
# cancer <- "brca"
# #####

# read files
switches <- read.table(switches.file,header=T)
switches$Symbol <- as.character(switches$Symbol)
switches$Normal_transcript <- as.character(switches$Normal_transcript)
switches$Tumor_transcript <- as.character(switches$Tumor_transcript)
switches$Patients_affected <- as.character(switches$Patients_affected)
psi.nt <- read.table(psi.nt.file, check.names=FALSE)
psi.t <- read.table(psi.t.file, check.names=FALSE)

psi <- cbind(psi.nt,psi.t)
rm(psi.nt,psi.t)

# get sample names
normal <- grep("^.{4}N$",colnames(psi), value=TRUE)
tumor <- grep("^.{4}T$",colnames(psi), value=TRUE)
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

  # get psi
  switch.psi <- psi[rownames(psi) %in% c(paste(gene,niso,sep=","),paste(gene,tiso,sep=",")),]
  switch.psi <- as.data.frame(t(switch.psi))
  
  switch.psi$What <- "Normal"
  switch.psi$What[rownames(switch.psi) %in% tumor] <- "Tumor no switch"
  switch.psi$What[rownames(switch.psi) %in% unlist(switch.patients)] <- "Tumor switch"
  
  switch.psi$Type <- "Paired"
  switch.psi$Type[rownames(switch.psi) %in% tumor.unpaired] <- "Unpaired"
  
  colnames(switch.psi)[colnames(switch.psi) == paste(gene,niso,sep=",")] <- "Normal"
  colnames(switch.psi)[colnames(switch.psi) == paste(gene,tiso,sep=",")] <- "Tumor"
  
  # plot
  g <- ggplot() + 
    geom_point(data=subset(switch.psi, What=="Tumor no switch"),aes(x=Normal,y=Tumor),color="#bdbdbd",size=2) +
    geom_point(data=subset(switch.psi, What!="Tumor no switch"),aes(x=Normal,y=Tumor,color=What,shape=Type),size=5) +
    smartas_theme() +
    labs(x="", y="") +
    scale_color_manual(values=c("Normal"="#3182bd", "Tumor switch"="#e6550d"))
  
  n <- ggplot(switch.psi,aes(x=Normal, fill=What)) + 
    geom_histogram(alpha=0.75) +
    labs(x="Normal isoform PSI",y="") +
    smartas_theme() +
    scale_fill_manual(values=c("Tumor no switch"="#7fc97f", "Normal"="#beaed4", "Tumor switch"="#fdc086"))
  
  t <- ggplot(switch.psi, aes(x=Tumor, fill=What)) + 
    geom_histogram(alpha=0.75) +
    smartas_theme() +
    labs(x="Tumor isoform PSI",y="") +
    coord_flip() +
    scale_fill_manual(values=c("Tumor no switch"="#7fc97f", "Normal"="#beaed4", "Tumor switch"="#fdc086"))
  
  x <- ggplotGrob(g + theme(legend.position="bottom") + 
                  labs(color="",shape="") +
                  guides(color=guide_legend(nrow=3,byrow=TRUE,override.aes = list(shape = 15, size = 8)),
                         shape=guide_legend(nrow=2,byrow=TRUE,override.aes = list(size = 8))))$grobs 
  
  legend <- x[[which(sapply(x, function(y) y$name) == "guide-box")]]
  
  p <- grid.arrange(t,g,legend,n,ncol=2,nrow=2, widths=c(1.5,5), heights=c(5,1.5),
                    top=textGrob(paste(gene,niso,tiso,sep=" "),gp=gpar(fontsize=20,font=3)))
  
  ggsave(paste(cancer,i,"psi",gene,niso,tiso,"png",sep="."),p,width = 12,height = 12)
}