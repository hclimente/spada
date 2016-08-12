#!/usr/bin/env Rscript

source("~/smartas/pipeline/scripts/variablesAndFunctions.r")
wd <- "~/smartas/analyses/pancancer/"
nb <- "~/smartas/notebook/data/"
library(tidyr)

#############################
##      STROMA/IMMUNE      ##
#############################

# switches with high correlation with:
## stromal content
stromaCorr <- lapply(cancerTypes, function(tumor){
  filename <- paste0(nb,"eporta/estimate/", tumor, "_ouptut_r_stroma.txt")
  read.delim(filename) %>%
    mutate(., Switch=row.names(.), Tumor=tumor, SwitchCorrelation="Stroma") %>%
    filter(FDR_linear < 0.05 | FDR_wilcox < 0.05) %>%
    separate(Switch,into=c("Normal_transcript","Tumor_transcript"), sep="-")
}) %>% do.call("rbind", .)

## immune content
immuneCorr <- lapply(cancerTypes, function(tumor){
  filename <- paste0(nb,"eporta/estimate/", tumor, "_ouptut_r_immune.txt")
  read.delim(filename) %>%
    mutate(., Switch=row.names(.), Tumor=tumor, SwitchCorrelation="Immune") %>%
    filter(FDR_linear < 0.05 | FDR_wilcox < 0.05) %>%
    separate(Switch,into=c("Normal_transcript","Tumor_transcript"), sep="-")
}) %>% do.call("rbind", .)

contaminantSwitches <- rbind(stromaCorr[,c("Normal_transcript","Tumor_transcript","SwitchCorrelation")],
                             immuneCorr[,c("Normal_transcript","Tumor_transcript","SwitchCorrelation")])

# cell-lineage specific genes
cellSpecificGenes <- paste0(nb,"eporta/estimate/geneset_ESTIMATE.txt") %>%
  read_tsv(col_names = FALSE) %>% 
  t %>%
  as.data.frame %>%
  set_colnames(c("StromalSignature","ImmuneSignature")) %>%
  .[3:nrow(.),]

stromaGenes <- data.frame(Symbol=cellSpecificGenes$StromalSignature, Genetype="Stroma")
immuneGenes <- data.frame(Symbol=cellSpecificGenes$ImmuneSignature, Genetype="Immune")

contaminantGenes <- rbind(stromaGenes,immuneGenes)

#############################
##        AGGREGATED       ##
#############################

# read switches
switches <- read_tsv(paste0(wd,"candidateList_info.agg.tsv"))
candidates <- read_tsv("~/smartas/notebook/data/switches/driverEvidence.tsv") %>%
  group_by(Normal_transcript,Tumor_transcript) %>%
  summarise(Recurrence = max(Recurrence), PPI = min(PPI),
            Affects_mutated_feature = max(Affects_mutated_feature),
            Pannegative = max(Pannegative))

wes <- read_tsv(paste0(wd,"mutations/gene_functional_mutations_all_switches.txt")) %>%
	set_colnames(c("GeneId","Symbol","Normal_transcript","Tumor_transcript","MS.pam","M.pam","S.pam","N.pam","p.pam.me"))
wgs <- read_tsv(paste0(wd,"mutations/gene_wgs_mutations_all_switches.txt")) %>%
	set_colnames(c("GeneId","Symbol","Normal_transcript","Tumor_transcript","MS.mut","M.mut","S.mut","N.mut","p.mut.o"))

# add candidate, wes and wgs information
switches <- merge(switches,candidates,all.x=T) %>%
	merge(wes) %>%
	merge(wgs,all.x=T) %>%
  # add candidate information
  mutate(Candidate = as.numeric((Recurrence+Affects_mutated_feature+PPI+Pannegative) > 0),
         Candidate = ifelse(IsFunctional == 1, Candidate, 0))

switches %>%
  left_join(contaminantSwitches) %>%
  left_join(contaminantGenes) %>%
  select(Normal_transcript, Tumor_transcript, Genetype,SwitchCorrelation) %>%
  mutate(Genetype = ifelse(is.na(Genetype), "Tumor", as.character(Genetype)),
         SwitchCorrelation = ifelse(is.na(SwitchCorrelation), "Tumor", SwitchCorrelation),
         # cases where the switch and the gene analysis do not agree are "Uncertain"
         Origin = ifelse(Genetype==SwitchCorrelation, SwitchCorrelation, "Uncertain"),
         Origin = ifelse(Origin=="Uncertain" & SwitchCorrelation=="Tumor", Genetype, Origin),
         Origin = ifelse(Origin=="Uncertain" & Genetype=="Tumor", SwitchCorrelation, Origin)) %>%
  group_by(Normal_transcript,Tumor_transcript) %>%
  summarise(Origin = paste0(sort(unique(Origin)), collapse="&")) %>%
  mutate(Origin = ifelse(grepl("&Uncertain",Origin), "Uncertain", Origin)) %>%
  merge(switches) %>%
  # rearrange columns
  select(GeneId,Symbol,Normal_transcript,Tumor_transcript,Normal_protein:IsFunctional,Origin,
         Driver:Pannegative, Candidate, MS.pam:p.mut.o) %>%
  write_tsv(paste0(wd,"candidateList_full.tsv"))

#############################
##       TUMOR-SPLIT       ##
#############################

switches.split <- read_tsv(paste0(wd,"candidateList_info.tumorSplit.tsv"))
candidates <- read_tsv("~/smartas/notebook/data/switches/driverEvidence.tsv")
wes <- read_tsv("~/smartas/notebook/data/mutations/gene_functional_mutations_all_switches.txt") %>%
  set_colnames(c("Tumor","GeneId","Symbol","Normal_transcript","Tumor_transcript","MS.pam","M.pam","S.pam","N.pam","p.pam.me","OR.pam","Switched","PAM"))
wgs <- read_tsv("~/smartas/notebook/data/mutations/gene_wgs_mutations_all_switches.txt") %>%
  set_colnames(c("Tumor","GeneId","Symbol","Normal_transcript","Tumor_transcript","MS.mut","M.mut","S.mut","N.mut","p.mut.o"))

candidateInfo <- switches %>%
  select(Normal_transcript,Tumor_transcript,Candidate,Reliable)

switches.split <- merge(switches.split,candidates,all.x=T) %>%
  merge(wes) %>%
  merge(wgs,all.x=T) %>%
  # add candidate information
  merge(candidateInfo)

switches.split %>%
  left_join(contaminantSwitches) %>%
  left_join(contaminantGenes) %>%
  select(Normal_transcript, Tumor_transcript, Genetype,SwitchCorrelation) %>%
  mutate(Genetype = ifelse(is.na(Genetype), "Tumor", as.character(Genetype)),
         SwitchCorrelation = ifelse(is.na(SwitchCorrelation), "Tumor", SwitchCorrelation),
         # cases where the switch and the gene analysis do not agree are "Uncertain"
         Origin = ifelse(Genetype==SwitchCorrelation, SwitchCorrelation, "Uncertain"),
         Origin = ifelse(Origin=="Uncertain" & SwitchCorrelation=="Tumor", Genetype, Origin),
         Origin = ifelse(Origin=="Uncertain" & Genetype=="Tumor", SwitchCorrelation, Origin)) %>%
  group_by(Normal_transcript,Tumor_transcript) %>%
  summarise(Origin = paste0(sort(unique(Origin)), collapse="&")) %>%
  mutate(Origin = ifelse(grepl("&Uncertain",Origin), "Uncertain", Origin)) %>%
  merge(switches.split) %>%
  # rearrange columns
  select(Tumor,GeneId,Symbol,Normal_transcript,Tumor_transcript,Normal_protein:IsFunctional,
         Reliable,Origin,Driver:Pannegative, Candidate, MS.pam:p.mut.o) %>%
  write_tsv(paste0(wd,"candidateList_full.tumorSplit.tsv")) 