source("~/smartas/pipeline/scripts/variablesAndFunctions.r")
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
wdp <- args[2]

# read switches
candidateInfo <- read_tsv(paste0(root,"notebook/data/pancancer/candidateList_full.tsv")) %>%
  mutate(Coocurrent= ifelse(is.na(Coocurrent),FALSE,Coocurrent),
         Candidate = as.logical(Recurrent) | as.logical(AffectingMutatedFeature) | as.logical(Coocurrent)) %>%
  mutate(Candidate = revalue(as.character(Candidate),replace=c("TRUE"="Candidate","FALSE"="Non-candidate"))) %>%
  select(GeneId,Symbol,Normal_transcript,Tumor_transcript,Recurrent,AffectingMutatedFeature,Coocurrent,Candidate)

switches <- read_tsv(paste0(wdp,"candidateList_info.tumorSplit.tsv"))
switches.agg <- read_tsv(paste0(wdp,"candidateList_info.agg.tsv"))

# read drivers
drivers <- read_tsv(paste0(root,"data/ucsc/intogen_cancer_drivers-2014.12b/Mutational_drivers_per_tumor_type.tsv"),comment="#")

# read mutation data
allPatients <- strsplit(switches.agg$Patients_affected,",") %>% unlist %>% unique

mutations <- read_tsv(paste0(root,"notebook/data/mutations/wes_mutations.txt")) %>%
  select(Tumor,Patient,Symbol) %>%
  unique %>%
  # filter out patients without RNAseq data available
  filter(Patient %in% allPatients) %>%
  mutate(Alteration2="MUT")

# read expression
proteome <- read_tsv(paste0(wdp,"mutations/proteome_information.txt"))

# read interactions
ppi.file <- paste0(root,"notebook/data/structural_analysis/Switched_interactions_consensus.txt")

## get max number of columns (necessary for reading)
no_col <- max(count.fields(ppi.file,sep = "\t"))
no_col.ppi <- (no_col-6)/2
ppi.cols <- paste(c("Origin","Interaction"), floor(seq(1,no_col.ppi,0.5)), sep="_")

## read table
ppi <- read.table(ppi.file,header=F,fill=T,col.names=1:no_col) %>%
  set_colnames(c("GeneId","Symbol","Normal_transcript","Tumor_transcript","partnerId","partnerSymbol",ppi.cols)) %>%
  # all Origin columns contail "DDI_match", so we can disregard them
  select(-starts_with("Origin_")) %>%
  # convert from wide to long table format
  melt(id.vars = c("GeneId","Symbol","Normal_transcript","Tumor_transcript","partnerId","partnerSymbol"),
       value.name = "Interaction") %>%
  select(-variable) %>%
  # remove cases with no interaction described
  filter(Interaction!="") %>%
  # split interaction information
  separate(Interaction, into=c("What","Transcript","Pfams"), sep="-") %>%
  # remove pfams columns (account for different domains for the same interaction)
  select(-Pfams) %>%
  # remove several instances of the same isoform
  unique %>%
  # annotate with switch info
  merge(switches) %>%
  merge(candidateInfo) %>%
  # consider only the most abundant isoform as partner: one interaction per pair & only expressed genes
  merge(proteome,by.x=c("Tumor","Transcript"), by.y=c("Cancer","Transcript"), suffixes = c(".switch",".partner"))

save.image(file=paste0(wdp,"ppi/data.RData"))
as.data.frame(unique(switches$Symbol)) %>% 
  write_tsv(paste0(wdp,"ppi/genes"))