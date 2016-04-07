library(plyr)

tryCatch(source("~/smartas/pipeline/scripts/variablesAndFunctions.r"),error=function(e){})

drivers.file <- paste0("~/smartas/data/ucsc/intogen_cancer_drivers-2014.12b/Mutational_drivers_per_tumor_type.tsv")

alterations <- list()
for (cancer in cancerTypes){
  
  # read alternative splicing switches
  cds.recurrence.file <- paste0("~/smartas/analyses/",cancer,"/candidateList_recurrence.tsv")
  cds.me.file <- paste0("~/smartas/analyses/",cancer,"/candidateList_mutationME.tsv")
  cds.co.file <- paste0("~/smartas/analyses/",cancer,"/candidateList_mutationCoocurrence.tsv")
  cds.features.file <- paste0("~/smartas/analyses/",cancer,"/candidateList_mutatedFeatures.tsv")
  
  cds.recurrence <- read.delim(cds.recurrence.file)
  cds.recurrence <- cds.recurrence[cds.recurrence$padj.recurrence < 0.05,c("GeneId","Symbol","Normal_transcript","Tumor_transcript")]
  cds.features <- read.delim(cds.features.file)
  cds.features <- cds.features[as.logical(cds.features$AffectingMutatedFeature),c("GeneId","Symbol","Normal_transcript","Tumor_transcript")]
  
  cds <- unique(rbind(cds.recurrence,cds.me,cds.features))
  
  if (file.exists(cds.co.file)){
    cds.co <- read.delim(cds.co.file)
    cds.co <- cds.co[cds.co$p.o < 0.05,c("GeneId","Symbol","Normal_transcript","Tumor_transcript")]
    cds <- rbind(cds,cds.co)
  }
  
  if (file.exists(cds.me.file)){
    cds.me <- read.delim(cds.me.file)
    cds.me <- cds.me[as.logical(cds.me$candidate),c("GeneId","Symbol","Normal_transcript","Tumor_transcript")]
    cds <- rbind(cds,cds.me)
  }

  # correct tag according to Intogene nomenclature
  if (cancer=="coad"){
    tag="COREAD"
  } else if  (cancer=="kich") {
    next
  } else if (cancer=="kirp") {
    next
  } else if (cancer=="kirc") {
    tag <- "RCCC"
  } else if (cancer=="lihc") {
    tag <- "HC"
  } else {
    tag <- toupper(cancer)
  }
  
  # read switches
  switches.file <- paste0("~/smartas/analyses/",cancer,"/candidateList_info.tsv")
  switches <- read.delim(switches.file,header=TRUE)  
  
  ## consider only AS drivers
  switches <- merge(switches,cds)
  
  ## count number of AS drivers switched per patient
  switchesPerPatient <- table(unlist(strsplit(as.character(switches$Patients_affected),",")))
  
  # get sample names
  psi.nt.file <- paste0("/projects_rg/TCGA/pipeline/run11/",cancer,"_iso_psi_paired-filtered.txt")
  psi.t.file <- paste0("/projects_rg/TCGA/pipeline/run11/",cancer,"_iso_psi_tumor-filtered.txt")
  
  patients.nt <- readLines(psi.nt.file, n=1)
  patients.nt <- unlist(strsplit(patients.nt,"\t"))
  patients.t <- readLines(psi.t.file, n=1)
  patients.t <- unlist(strsplit(patients.t,"\t"))
  
  patients <- c(patients.nt,patients.t)
  
  normal <- grep("^.{4}N$",patients, value=TRUE)
  tumor.paired <- gsub("N$","T",normal)
  
  # read mutations  
  ## get cancer specific drivers
  drivers <- read.delim(drivers.file,skip=5,row.names=NULL,header=T)
  drivers <- drivers[drivers$Tumor_type==tag,]
  
  ## read mutations
  mutations.file <- paste0("~/smartas/data/ucsc/rawdata/",cancer,"_gene_mutation-functional-count_full.txt")
  mutations <- read.delim(mutations.file,header=F)
  colnames(mutations) <- c("chr","start","end","gene","wut1","wut2","wut3","wut4","wut5","alteration","wut6","wut7")
  mutations <- mutations[mutations$alteration!=".",]
  mutations$Symbol <- unlist(strsplit(as.character(mutations$gene),"|",fixed=T))[c(T,F)]
  mutations$Patient <- unlist(strsplit(as.character(mutations$alteration),";"))[c(T,F)]
  
  ## filter out mutations not in drivers
  mutations.drivers <- mutations[mutations$Symbol %in% drivers$geneHGNCsymbol,]
  
  ## consider only mutated genes, not several mutations in a gene
  mutations.drivers <- unique(mutations.drivers[,c("gene","Patient")])
  
  ## count number of mutations per patient
  mutationsPerPatient <- table(mutations.drivers$Patient)
  
  # get patients common between paired patients screened for switches, 
  # patients screened for mutations
  patients.mut <- unique(mutations$Patient)
  patients.common <- intersect(tumor.paired,patients.mut)
  
  ## filter out those patients without full information
  switchesPerPatient <- switchesPerPatient[patients.common]
  names(switchesPerPatient) <- patients.common
  switchesPerPatient[is.na(switchesPerPatient)] <- 0
  switchesPerPatient <- switchesPerPatient/nrow(switches)
  
  mutationsPerPatient <- mutationsPerPatient[patients.common]
  names(mutationsPerPatient) <- patients.common
  mutationsPerPatient[is.na(mutationsPerPatient)] <- 0
  mutationsPerPatient <- mutationsPerPatient/nrow(drivers)
  
  # merge lists
  alterationsPerPatient <- as.data.frame(t(rbind.fill.matrix(t(switchesPerPatient), t(mutationsPerPatient))))
  colnames(alterationsPerPatient) <- c("Switches","Mutations")
  alterationsPerPatient$Tumor <- cancer
  
  alterations[[cancer]] <- alterationsPerPatient
  
}

alterationsPerPatient <- do.call("rbind",alterations)
alterationsPerPatient.m <- melt(alterationsPerPatient)
alterationsPerPatient.z <- ddply(alterationsPerPatient,.(Tumor),summarise,
                                 Switches.z=(Switches-mean(Switches))/sd(Switches),
                                 Mutations.z=(Mutations-mean(Mutations))/sd(Mutations))

# plot raw frequency
p <- ggplot(alterationsPerPatient,aes(x=Switches*100,y=Mutations*100,color=Tumor)) + 
  geom_point() + 
  smartas_theme() +
  labs(x="% splicing drivers switched",y="% mutational drivers mutated") +
  scale_color_manual(values=colorPalette) +
  theme(legend.position="bottom")

out.scatterplot <- paste0("~/smartas/results/switches/figures/mutated_vs_switched_patients.png")
ggsave(out.scatterplot,p,width = 12,height = 12)

# plot distributions
q <- ggplot(alterationsPerPatient.m) + 
  geom_boxplot(aes(x=Tumor,y=value,fill=Tumor,alpha=variable)) + 
  smartas_theme() +
  labs(x="Tumor",y="% drivers altered/patient",alpha="Driver type") +
  scale_fill_manual(values=colorPalette) +
  scale_alpha_manual(values=c(1,0.6)) +
  theme(legend.position="bottom")

out.boxplot <- paste0("~/smartas/results/switches/figures/mutated_and_switched_patients_boxplot.png")
ggsave(out.boxplot,q,width = 12,height = 12)

# plot frequency corrected by background frequency (z-score)
# consider that background information might not need to be corrected
# e.g. if a cancer has less mutations might be because relies more on AS
z <- ggplot(alterationsPerPatient.z,aes(x=Switches.z,y=Mutations.z,color=Tumor)) + 
  geom_point() + 
  smartas_theme() +
  labs(x="Z(% splicing drivers switched)",y="Z(% mutational drivers mutated)") +
  scale_color_manual(values=colorPalette) +
  theme(legend.position="bottom")

out.z.scatterplot <- paste0("~/smartas/results/switches/figures/mutated_vs_switched_patients_z.png")
ggsave(out.z.scatterplot,z,width = 12,height = 12)