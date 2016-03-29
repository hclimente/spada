#!/soft/R/R-3.2.3/bin/Rscript

library(ggplot2)
library(grid)
library(gridExtra)

tryCatch(source("~/smartas/pipeline/scripts/variablesAndFunctions.r"),error=function(e){})

args <- commandArgs(trailingOnly = TRUE)

switches.file <- args[1]
gene <- args[2]
niso <- args[3]
tiso <- args[4]
tpm.nt.file <- args[5]
tpm.t.file <- args[6]

# test
# switches.file <- "~/smartas/analyses/brca/candidateList.tsv"
# gene <- "CDK6|1021"
# niso <- "uc011khw.1"
# tiso <- "uc010lez.2"
# tpm.nt.file <- "/projects_rg/TCGA/pipeline/run11/brca_iso_tpm_paired-filtered.txt"
# tpm.t.file <- "/projects_rg/TCGA/pipeline/run11/brca_iso_tpm_tumor-filtered.txt"
#####

# read files
switches <- read.table(switches.file)
colnames(switches) <- c("Gene","Normal","Tumor","Patients")
tpm.nt <- read.table(tpm.nt.file, check.names=FALSE)
tpm.t <- read.table(tpm.t.file, check.names=FALSE)

tpm <- cbind(tpm.nt,tpm.t)
rm(tpm.nt,tpm.t)

# get sample names
normal <- grep("^.{4}N$",colnames(tpm), value=TRUE)
tumor <- grep("^.{4}T$",colnames(tpm), value=TRUE)
tumor.paired <- gsub("N$","T",normal)
tumor.unpaired <- setdiff(tumor,tumor.paired)

# subset
## get patients
x <- strsplit(as.character(switches$Patients),",")
p <- which(switches$Gene==gene & switches$Normal==niso & switches$Tumor==tiso)

switch.patients <- x[p]

# get tpm
switch.tpm <- tpm[rownames(tpm) %in% c(paste(gene,niso,sep=","),paste(gene,tiso,sep=",")),]
switch.tpm <- as.data.frame(t(switch.tpm))

switch.tpm$What <- "Normal"
switch.tpm$What[rownames(switch.tpm) %in% tumor] <- "Tumor no switch"
switch.tpm$What[rownames(switch.tpm) %in% unlist(switch.patients)] <- "Tumor switch"

colnames(switch.tpm)[colnames(switch.tpm) == paste(gene,niso,sep=",")] <- "Normal"
colnames(switch.tpm)[colnames(switch.tpm) == paste(gene,tiso,sep=",")] <- "Tumor"

# plot
g <- ggplot(switch.tpm,aes(x=Normal,y=Tumor,color=What)) + 
  geom_point() +
  smartas_theme() +
  labs(x="", y="") +
  scale_color_manual(values=c("Tumor no switch"="#7fc97f", "Normal"="#beaed4", "Tumor switch"="#fdc086"))

n <- ggplot(switch.tpm,aes(x=Normal, fill=What)) + 
  geom_histogram(alpha=0.75) +
  labs(x="Normal isoform TPM",y="") +
  smartas_theme() +
  scale_fill_manual(values=c("Tumor no switch"="#7fc97f", "Normal"="#beaed4", "Tumor switch"="#fdc086"))

t <- ggplot(switch.tpm, aes(x=Tumor, fill=What)) + 
  geom_histogram(alpha=0.75) +
  smartas_theme() +
  labs(x="Tumor isoform TPM",y="") +
  coord_flip() +
  scale_fill_manual(values=c("Tumor no switch"="#7fc97f", "Normal"="#beaed4", "Tumor switch"="#fdc086"))

x <- ggplotGrob(g + theme(legend.position="bottom") + 
                labs(color="") +
                  guides(color=guide_legend(nrow=3,byrow=TRUE,override.aes = list(shape = 15, size = 8)),
                         shape=guide_legend(nrow=2,byrow=TRUE,override.aes = list(size = 8))))$grobs 

legend <- x[[which(sapply(x, function(y) y$name) == "guide-box")]]

p <- grid.arrange(t,g,legend,n,ncol=2,nrow=2, widths=c(1.5,5), heights=c(5,1.5),
                  top=textGrob(paste(gene,niso,tiso,sep=" "),gp=gpar(fontsize=20,font=3)))

ggsave(paste("tpm",gene,niso,tiso,"png",sep="."),p,width = 12,height = 12)