#!/usr/bin/env Rscript

source("~/smartas/pipeline/scripts/variablesAndFunctions.r")
library(tidyr)
wd <- "~/smartas/analyses/pancancer/"
nb <- "~/smartas/notebook/data/"

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

read_tsv(paste0(wd,"candidateList_info.agg.tsv")) %>%
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
  select(Normal_transcript,Tumor_transcript,Origin) %>%
  write_tsv(paste0(wd,"candidateList.origin.tsv"))
