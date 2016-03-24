library(reshape2)
library(plyr)
library(gtools)

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
#out <- args[1]
out <- "~/smartas/analyses/coad/"

proteome.muts.file <- paste0(out,"mutations/proteome_mutations.txt")
proteome.fts.file <- paste0(out,"mutations/proteome_features.txt")
switch.prosite.file <- paste0(out,"structural_analysis/prosite_analysis.tsv")
switch.pfam.file <- paste0(out,"structural_analysis/interpro_analysis.tsv")
switch.file <- paste0(out,"candidateList_info.tsv")
out.file <- paste0(out,"candidateList_mutatedFeatures.tsv")

switches <- read.delim(switch.file,row.names=NULL)

# read mutations
proteome.muts <- read.delim(proteome.muts.file,row.names=NULL)
allMuts <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Frame_Shift_Del_out","Frame_Shift_Ins_out","Nonsense_Mutation_out")
inFeatureMuts <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation")

proteome.muts.wide <- dcast(data=proteome.muts,formula=Feature+Analysis+Cancer+Transcript~Type)

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
proteome.fts <- read.delim(proteome.fts.file,row.names=NULL)

# read switches
switch.prosite <- read.delim(switch.prosite.file,row.names=NULL)
switch.prosite <- subset(switch.prosite, What!="Nothing")
switch.pfam <- read.delim(switch.pfam.file,row.names=NULL)
switch.pfam <- subset(switch.pfam, What!="Nothing")
switch.fts <- rbind(switch.pfam,switch.prosite)

switchesFull <- merge(switches,subset(switch.fts, select=c("Gene","Symbol","NormalTranscript","TumorTranscript","What","Feature")),by.x=c("GeneId","Symbol","Normal_transcript","Tumor_transcript"),by.y=c("Gene","Symbol","NormalTranscript","TumorTranscript"))

# aggregate by transcript and cancer to get the proteome size and expected frequencies
tests <- list()
for (f in c("prosite","Pfam")){
  ft.agg <- ddply(subset(proteome.fts, Analysis==f),.(Feature),summarize, TotalLength = sum(FeatureLength) )
  
  # corrected: only mutations in domains are counted; only domains should be counted
  L <- sum(ft.agg$TotalLength)
  ft.agg$ExpectedMutFrequency <- ft.agg$TotalLength/L
  
  ## Test features enriched in mutations
  # aggregate mutations
  ft.mut.enrichment <- subset(proteome.muts.agg, Analysis==f)
  ft.allMutations <- sum(ft.mut.enrichment$TotalMutations)
  
  # calculate statistics
  ft.mut.enrichment <- merge(ft.mut.enrichment,ft.agg)
  ft.mut.enrichment$fc_m <- ft.mut.enrichment$TotalMutations/ft.allMutations/ft.mut.enrichment$ExpectedMutFrequency
  ft.mut.enrichment$p_m <- apply(ft.mut.enrichment[,c("TotalMutations","ExpectedMutFrequency")],1, my.binom.test, ft.allMutations)
  ft.mut.enrichment$adjp_m <- p.adjust(ft.mut.enrichment$p_m)
  
  write.table(ft.mut.enrichment[,c("Feature","TotalMutations","TotalLength","ExpectedMutFrequency","fc_m","p_m","adjp_m","H_m","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Frame_Shift_Del_out","Frame_Shift_Ins_out","Nonsense_Mutation_out")],paste0(out,"mutations/",f,".mutation.enrichment.txt"),sep="\t",row.names = FALSE,quote=F)
  
  tests[[f]] <- ft.mut.enrichment
}

enrichment <- do.call("rbind",tests)

df <- unique(switchesFull[switchesFull$Feature %in% enrichment$Feature[enrichment$adjp_m<0.05 & enrichment$H_m > 0],c("GeneId","Symbol","Normal_transcript","Tumor_transcript")])
switches$AffectingMutatedFeature <- 0
switches$AffectingMutatedFeature[switches$Normal_transcript %in% df$Normal_transcript & switches$Tumor_transcript %in% df$Tumor_transcript] <- 1

write.table(switches[,c("GeneId","Symbol","Normal_transcript","Tumor_transcript","AffectingMutatedFeature")],out.file,sep="\t",row.names = FALSE,quote=F)