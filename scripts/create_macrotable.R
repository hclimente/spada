#!/usr/bin/env Rscript

source("~/smartas/pipeline/scripts/variablesAndFunctions.r")

wd <- "~/smartas/analyses/pancancer/"

switches <- read_tsv(paste0(wd,"candidateList_info.agg.tsv"))
candidates <- read_tsv("~/smartas/notebook/data/switches/driverEvidence.tsv") %>%
                         group_by(Normal_transcript,Tumor_transcript) %>%
                         summarise(Recurrence = max(Recurrence), DriverME = max(DriverME),
                                   Affects_mutated_feature = max(Affects_mutated_feature),
                                   PPI = max(PPI), Pannegative = max(Pannegative)) %>%
                         as.data.frame
wes <- read_tsv(paste0(wd,"mutations/gene_functional_mutations_all_switches.txt")) %>%
	set_colnames(c("GeneId","Symbol","Normal_transcript","Tumor_transcript","MS.pam","M.pam","S.pam","N.pam","p.pam.me"))
wgs <- read_tsv(paste0(wd,"mutations/gene_wgs_mutations_all_switches.txt")) %>%
	set_colnames(c("GeneId","Symbol","Normal_transcript","Tumor_transcript","MS.mut","M.mut","S.mut","N.mut","p.mut.o"))

x <- merge(switches,candidates,all.x=T) %>%
	merge(wes) %>%
	merge(wgs,all.x=T)

write_tsv(x,paste0(wd,"candidateList_full.tsv"))