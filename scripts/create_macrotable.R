#!/usr/bin/env Rscript

source("~/smartas/pipeline/scripts/variablesAndFunctions.r")
wd <- "~/smartas/analyses/pancancer/"

#############################
##        AGGREGATED       ##
#############################

switches <- read_tsv(paste0(wd,"candidateList_info.agg.tsv"))
candidates <- read_tsv("~/smartas/notebook/data/switches/driverEvidence.tsv") %>%
                         group_by(Normal_transcript,Tumor_transcript) %>%
                         summarise(Recurrence = max(Recurrence), PPI = max(PPI),
                                   Affects_mutated_feature = max(Affects_mutated_feature),
                                   Pannegative = max(Pannegative))
wes <- read_tsv(paste0(wd,"mutations/gene_functional_mutations_all_switches.txt")) %>%
	set_colnames(c("GeneId","Symbol","Normal_transcript","Tumor_transcript","MS.pam","M.pam","S.pam","N.pam","p.pam.me"))
wgs <- read_tsv(paste0(wd,"mutations/gene_wgs_mutations_all_switches.txt")) %>%
	set_colnames(c("GeneId","Symbol","Normal_transcript","Tumor_transcript","MS.mut","M.mut","S.mut","N.mut","p.mut.o"))

merge(switches,candidates,all.x=T) %>%
	merge(wes) %>%
	merge(wgs,all.x=T) %>%
  write_tsv(paste0(wd,"candidateList_full.tsv"))

#############################
##       TUMOR-SPLIT       ##
#############################

switches <- read_tsv(paste0(wd,"candidateList_info.tumorSplit.tsv"))
candidates <- read_tsv("~/smartas/notebook/data/switches/driverEvidence.tsv")
wes <- read_tsv("~/smartas/notebook/data/mutations/gene_functional_mutations_all_switches.txt") %>%
  set_colnames(c("Tumor","GeneId","Symbol","Normal_transcript","Tumor_transcript","MS.pam","M.pam","S.pam","N.pam","p.pam.me","OR.pam","Switched","PAM"))
wgs <- read_tsv("~/smartas/notebook/data/mutations/gene_wgs_mutations_all_switches.txt") %>%
  set_colnames(c("Tumor","GeneId","Symbol","Normal_transcript","Tumor_transcript","MS.mut","M.mut","S.mut","N.mut","p.mut.o"))

merge(switches,candidates,all.x=T) %>%
  merge(wes) %>%
  merge(wgs,all.x=T) %>%
  write_tsv(paste0(wd,"candidateList_full.tumorSplit.tsv"))