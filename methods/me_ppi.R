source("~/smartas/pipeline/scripts/variablesAndFunctions.r")
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
wdp <- args[2]

#############################
##        SWITCHES         ##
#############################

candidateInfo <- read_tsv(paste0(root,"notebook/data/pancancer/candidateList_full.tsv")) %>%
  mutate(Coocurrent= ifelse(is.na(Coocurrent),FALSE,Coocurrent),
         Candidate = as.logical(Recurrent) | as.logical(AffectingMutatedFeature) | as.logical(Coocurrent)) %>%
  mutate(Candidate = revalue(as.character(Candidate),replace=c("TRUE"="Candidate","FALSE"="Non-candidate"))) %>%
  select(GeneId,Symbol,Normal_transcript,Tumor_transcript,Recurrent,AffectingMutatedFeature,Coocurrent,Candidate)

switches <- read_tsv(paste0(wdp,"candidateList_info.tumorSplit.tsv"))
switches.agg <- read_tsv(paste0(wdp,"candidateList_info.agg.tsv"))

#############################
##         DRIVERS         ##
#############################

drivers <- read_tsv(paste0(root,"data/ucsc/intogen_cancer_drivers-2014.12b/Mutational_drivers_per_tumor_type.tsv"),comment="#") %>%
  mutate(Tumor_type = ifelse(Tumor_type=="COREAD", "coad", Tumor_type),
         Tumor_type = ifelse(Tumor_type=="HC", "lihc", Tumor_type),
         Tumor_type = ifelse(Tumor_type=="RCCC", "kirc", Tumor_type),
         Tumor_type = tolower(Tumor_type) ) %>%
  set_colnames(c("Symbol","Tumor"))

#############################
##        MUTATIONS        ##
#############################

allPatients <- strsplit(switches.agg$Patients_affected,",") %>% unlist %>% unique

mutations <- read_tsv(paste0(root,"notebook/data/mutations/wes_mutations.txt")) %>%
  select(Tumor,Patient,Symbol) %>%
  unique %>%
  # filter out patients without RNAseq data available
  filter(Patient %in% allPatients) %>%
  mutate(Alteration2="MUT")

#############################
##       EXPRESSION        ##
#############################

proteome <- read_tsv(paste0(wdp,"mutations/proteome_information.txt"))

#############################
##      INTERACTIONS       ##
#############################

ppi.file <- paste0(root,"notebook/data/structural_analysis/Switched_interactions_consensus.txt")

# get max number of columns (necessary for reading)
no_col <- max(count.fields(ppi.file,sep = "\t"))
no_col.ppi <- (no_col-6)/2
ppi.cols <- paste(c("Origin","Interaction"), floor(seq(1,no_col.ppi,0.5)), sep="_")

# read table
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

#############################
##   CREATE WIDE MATRIX    ##
#############################

# select the appropriate columns
affectedSwitches <- switches %>%
  filter(IsFunctional==1) %>%
  select(Tumor,Symbol,Patients_affected,PatientNumber)

# create long format matrix for switches
## split switch patients
nocols <- max(affectedSwitches$PatientNumber)
namescols <- paste("Patient", 1:nocols, sep="_")
suppressWarnings( affectedSwitches.long <- affectedSwitches %>%
                    separate(Patients_affected, namescols, sep=",") )

affectedSwitches.long <- affectedSwitches.long %>%
  select(Symbol,Tumor,starts_with("Patient_")) %>%
  melt(id.vars = c("Symbol","Tumor")) %>%
  select(-variable) %>%
  set_colnames(c("Symbol","Tumor","Patient")) %>%
  filter(!is.na(Patient)) %>%
  # add "SPLICING" as tag
  mutate(Alteration="SPLICING")

# create long matrix for mutations
affectedMutations.long <- mutations %>%
  merge(drivers)

# merge both matrixes
affected.long <- merge(affectedSwitches.long, affectedMutations.long,all=T) %>%
  mutate(Alteration=ifelse(is.na(Alteration),"",Alteration),
         Alteration2=ifelse(is.na(Alteration2),"",Alteration2)) %>%
  mutate(Alteration=paste(Alteration,Alteration2,sep=";"))

# convert to wide format
affected.wide <- affected.long %>%
  select(Symbol,Patient,Alteration) %>%
  spread(Patient, Alteration) %>%
  set_rownames(.$Symbol) %>%
  select(-Symbol) %>%
  as.matrix

## substitute NA's
affected.wide[is.na(affected.wide)] <- ""

write_tsv(affected.wide, paste0(wdp,"ppi/mutation_matrix.tsv"))

#############################
##    CALCULATE OVERLAP    ##
##     WITH MUTATIONS      ##
#############################

mutation.overlap <- list()

for (gene in unique(drivers$Symbol)){
  # add switches that affect the target
  studiedGenes <- ppi %>%
    # add interactors
    filter(partnerSymbol==gene & What!="Kept") %>%
    .$Symbol.switch %>%
    as.character %>%
    c(gene) %>%
    unique %>%
    # pick only those with a mutation or switch
    intersect(rownames(affected.wide))
  
  if(length(studiedGenes) < 2)
    next
  
  this.affected.wide <- affected.wide[studiedGenes,]
  
  # remove those genes with alterations that affect less than 5% of the patients
  ## remove empty columns
  this.affected.wide <- this.affected.wide[,colSums(this.affected.wide != "") > 0]
  minPatients <- floor(0.05 * ncol(this.affected.wide))
  
  if(sum(rowSums(this.affected.wide != "") > minPatients) < 2)
    next
  
  ## remove new empty columns
  this.affected.wide <- this.affected.wide[,colSums(this.affected.wide != "") > 0]
  
  nonMutated <- affected.wide[, colnames(this.affected.wide)]
  nonMutated[!grepl("MUT",nonMutated)] <- ""
  mutatedPatients <- sum(colSums(nonMutated != "") > 0)/ncol(nonMutated)
  
  mutation.overlap[[gene]] <- c(length(studiedGenes), ncol(this.affected.wide), mutatedPatients)
}

mutation.overlap <- do.call("rbind",mutation.overlap) %>%
  as.data.frame %>%
  set_colnames(c("num.genes","num.patients","prop.mutated")) %>%
  mutate(., Driver=rownames(.))

write_tsv(mutation.overlap, paste0(wdp,"ppi/me_ppi.tsv"))