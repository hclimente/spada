#!/usr/bin/env Rscript

source("~/smartas/pipeline/scripts/variablesAndFunctions.r")
wd <- "~/smartas/analyses/pancancer/"
nb <- "~/smartas/notebook/data/"
library(tidyr)

# read stroma/immune information
origin <- read_tsv(paste0(wd,"candidateList.origin.tsv"))
recurrence <- read_tsv(paste0(wd,"candidateList_recurrence.tsv")) %>%
  mutate(Spurious = (padj.recurrence < 0.05) & (what == "less")) %>%
  select(Normal_transcript,Tumor_transcript,Spurious)

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
  merge(origin) %>%
  merge(recurrence) %>%
  mutate(Reliable = ifelse(Spurious, 0, Reliable)) %>%
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
  merge(origin) %>%
  merge(recurrence) %>%
  mutate(EnoughRecurrence = ifelse(Spurious, 0, 1)) %>%
  # rearrange columns
  select(Tumor,GeneId,Symbol,Normal_transcript,Tumor_transcript,Normal_protein:IsFunctional,
         EnoughRecurrence,Origin,Driver:Pannegative, Candidate, MS.pam:p.mut.o) %>%
  write_tsv(paste0(wd,"candidateList_full.tumorSplit.tsv"))
