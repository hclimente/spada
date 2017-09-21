library(reshape2)
library(plyr)
library(dplyr)
library(gtools)
library(readr)
library(magrittr)

# calculate entropy among all genes with the domain
shannon.entropy <- function(p){
  if (invalid(p) || min(p) < 0 || sum(p) <= 0)
    return(NaN)
  p.norm <- p[p>0]/sum(p)
  H <- -sum(log2(p.norm)*p.norm)
  if (H == 0){
    0
  } else {
    maxH <- log2(length(p))
    H/maxH
  }
}

get.entropy <- function(m){
  totalMuts <- sum(m)
  p <- m/totalMuts
  shannon.entropy(p)
}

my.binom.test <- function(x,testNumber){ 
  p <- binom.test(x[1],testNumber,x[2],"greater")
  p$p.value
}

args <- commandArgs(trailingOnly = TRUE)
out <- args[1]

proteome.muts.file <- paste0(out,"mutations/proteome_mutations.txt")
proteome.fts.file <- paste0(out,"mutations/proteome_features.txt")
switch.prosite.file <- paste0(out,"structural_analysis/prosite_analysis.tsv")
switch.pfam.file <- paste0(out,"structural_analysis/interpro_analysis.tsv")
out.file <- paste0(out,"candidateList_mutatedFeatures.tsv")

switch.file <- paste0(out,"candidateList_info.tsv")
if(!file.exists(switch.file)){
  switch.file <- paste0(out,"candidateList_info.agg.tsv")
}

switches <- read_tsv(switch.file)

# read mutations
proteome.muts <- read_tsv(proteome.muts.file)
allMuts <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Frame_Shift_Del_out","Frame_Shift_Ins_out","Nonsense_Mutation_out")
inFeatureMuts <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation")

proteome.muts.wide <- dcast(data=proteome.muts,formula=Feature+Analysis+Cancer+Transcript~Type,fun.aggregate=length,value.var="Patient")

absent.mutTypes <- setdiff(c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"),proteome.muts$Type)
if (length(absent.mutTypes)!=0)
  proteome.muts.wide[,absent.mutTypes] <- 0
proteome.muts.wide$TotalMutations <- rowSums(proteome.muts.wide[,inFeatureMuts])

proteome.muts.agg <- ddply(proteome.muts.wide,.(Feature,Analysis),summarize, 
                           Frame_Shift_Del = sum(Frame_Shift_Del),
                           Frame_Shift_Ins = sum(Frame_Shift_Ins),
                           In_Frame_Del = sum(In_Frame_Del),
                           In_Frame_Ins = sum(In_Frame_Ins),
                           Missense_Mutation = sum(Missense_Mutation),
                           Nonsense_Mutation = sum(Nonsense_Mutation),
                           Nonstop_Mutation = sum(Nonstop_Mutation),
                           Frame_Shift_Del_out = sum(Frame_Shift_Del_out),
                           Frame_Shift_Ins_out = sum(Frame_Shift_Ins_out),
                           Nonsense_Mutation_out = sum(Nonsense_Mutation_out),
                           H_m = get.entropy(TotalMutations),
                           TotalMutations = sum(TotalMutations))

# read features
proteome.fts <- read_tsv(proteome.fts.file)

# read switches
switch.prosite <- read_tsv(switch.prosite.file) %>%
  filter(What!="Nothing")
switch.pfam <- read_tsv(switch.pfam.file) %>%
  filter(What!="Nothing")

switch.fts <- rbind(switch.pfam,switch.prosite)

switchesFull <- merge(switches,subset(switch.fts, select=c("Gene","Symbol","NormalTranscript","TumorTranscript","What","Feature")),by.x=c("GeneId","Symbol","Normal_transcript","Tumor_transcript"),by.y=c("Gene","Symbol","NormalTranscript","TumorTranscript"))

# aggregate by transcript and cancer to get the proteome size and expected frequencies
tests <- list()
for (f in c("prosite","Pfam")){
  ft.agg <- proteome.fts %>%
    filter(Analysis==f) %>%
    group_by(Feature) %>%
    summarize(TotalLength = sum(FeatureLength))
  
  # corrected: only mutations in domains are counted; only domains should be counted
  L <- sum(ft.agg$TotalLength)
  ft.agg$ExpectedMutFrequency <- ft.agg$TotalLength/L
  
  ## Test features enriched in mutations
  # aggregate mutations
  ft.mut.enrichment <- proteome.muts.agg %>%
    filter(Analysis==f) %>%
    merge(ft.agg)
  ft.allMutations <- sum(ft.mut.enrichment$TotalMutations)
  
  # calculate statistics
  ft.mut.enrichment <- ft.mut.enrichment %>%
    mutate(., fc_m = TotalMutations/ft.allMutations/ExpectedMutFrequency,
           p_m = apply(.[,c("TotalMutations","ExpectedMutFrequency")],1, my.binom.test, ft.allMutations)) %>%
    mutate(adjp_m = p.adjust(p_m))
  
  ft.mut.enrichment %>%
    select(Feature,TotalMutations,TotalLength,ExpectedMutFrequency,fc_m,p_m,adjp_m,
           H_m,Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,
           Missense_Mutation,Nonsense_Mutation,Nonstop_Mutation,Frame_Shift_Del_out,
           Frame_Shift_Ins_out,Nonsense_Mutation_out) %>%
  write_tsv(paste0(out,"mutations/",f,".mutation.enrichment.txt"))
  
  tests[[f]] <- ft.mut.enrichment
}

enrichment <- do.call("rbind",tests)

df <- switchesFull %>%
  filter(Feature %in% enrichment$Feature[enrichment$adjp_m<0.05 & enrichment$H_m > 0]) %>%
  select(c(GeneId,Symbol,Normal_transcript,Tumor_transcript)) %>%
  unique

switches$AffectingMutatedFeature <- 0
switches$AffectingMutatedFeature[switches$Normal_transcript %in% df$Normal_transcript & switches$Tumor_transcript %in% df$Tumor_transcript] <- 1

switches %>%
  select(GeneId,Symbol,Normal_transcript,Tumor_transcript,AffectingMutatedFeature) %>%
  write_tsv(out.file)