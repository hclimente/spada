library(ggplot2)
library(grid)
library(gridExtra)

tryCatch(source("~/smartas/pipeline/scripts/variablesAndFunctions.r"),error=function(e){})

args <- commandArgs(trailingOnly = TRUE)

switches.file <- args[1]
gene <- args[2]
niso <- args[3]
tiso <- args[4]
psi.nt.file <- args[5]
psi.t.file <- args[6]

# test
switches.file <- "~/smartas/analyses/brca/candidateList.tsv"
gene <- "GFRA1|2674"
niso <- "uc001lci.2"
tiso <- "uc009xyr.2"
psi.nt.file <- "/projects_rg/TCGA/pipeline/run11/brca_iso_psi_paired-filtered.txt"
psi.t.file <- "/projects_rg/TCGA/pipeline/run11/brca_iso_psi_tumor-filtered.txt"
#####

# read files
switches <- read.table(switches.file)
colnames(switches) <- c("Gene","Normal","Tumor","Patients")
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
x <- strsplit(as.character(switches$Patients),",")
p <- which(switches$Gene==gene & switches$Normal==niso & switches$Tumor==tiso)

switch.patients <- x[p]

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
g <- ggplot(switch.psi,aes(x=Normal,y=Tumor,color=What,shape=Type)) + 
  geom_point() +
  smartas_theme() +
  labs(x="", y="") +
  scale_color_manual(values=c("Tumor no switch"="#7fc97f", "Normal"="#beaed4", "Tumor switch"="#fdc086"))

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

ggsave(paste("psi",gene,niso,tiso,"png",sep="."),p,width = 12,height = 12)