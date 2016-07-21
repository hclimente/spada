#!/usr/bin/env Rscript

library(find.me)
library(tidyr)
source("~/smartas/pipeline/scripts/variablesAndFunctions.r")

args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
wdp <- args[2]
gene <- args[3]

# root <- "~/smartas/"
# wdp <- paste0(root,"analyses/pancancer/")
# gene <- "TAF9"

load(file=paste0(wdp,"ppi/data.RData"))

# add switches that affect the target
affectedSwitches <- ppi %>%
  filter(partnerSymbol==gene & What!="Kept") %>%
  select(Tumor,Symbol.switch,Patients_affected,PatientNumber) %>%
  set_colnames(c("Tumor","Symbol","Patients_affected","PatientNumber"))

# if present, add switches in the target gene
affectedSwitches <- switches %>%
  filter(Symbol==gene & IsFunctional==1) %>%
  select(Tumor,Symbol,Patients_affected,PatientNumber) %>%
  rbind(affectedSwitches)

# if we do not have at least two genes, skip
stopifnot(length(unique(affectedSwitches$Symbol)) >= 2)

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
  mutate(Alteration="SPLICING")

affectedMutations.long <- mutations %>%
  filter(Symbol %in% affectedSwitches.long$Symbol & Symbol %in% drivers$geneHGNCsymbol)

affected.long <- merge(affectedSwitches.long, affectedMutations.long,all=T) %>%
  mutate(Alteration=ifelse(is.na(Alteration),"",Alteration),
         Alteration2=ifelse(is.na(Alteration2),"",Alteration2)) %>%
  mutate(Alteration=paste(Alteration,Alteration2,sep=";"))

affected.wide <- affected.long %>%
  select(Symbol,Patient,Alteration) %>%
  spread(Patient, Alteration)

affected.wide[is.na(affected.wide)] <- ""

rownames(affected.wide) <- affected.wide$Symbol
affected.wide <- as.matrix(affected.wide)
affected.wide <- affected.wide[, !colnames(affected.wide) %in% c("Tumor","Symbol")]
affected.wide <- affected.wide[,colSums(affected.wide=="") < nrow(affected.wide)]

# remove those genes with alterations that affect less than 5% of the patients
minPatients <- floor(0.05 * ncol(affected.wide))
affected.wide <- affected.wide[rowSums(affected.wide != "") > minPatients,]
## remove new empty columns
affected.wide <- affected.wide[,colSums(affected.wide != "") > 0]

mutmat <- getSortedMatrix(affected.wide)$mutmat
# add empty matrix by the side
mutmat <- cbind(mutmat,matrix(0,nrow=nrow(mutmat), ncol=nPatients$total - nrow(mutmat)))
sampling.pval <- me.test.permutateSamples(mutmat)
fisher.pval <- me.test.fisher(mutmat)

if (min(sampling.pval,fisher.pval) < 0.05){
  
  # get patients for tumor bar at top
  patients <- affected.long %>% 
    filter(Patient %in% colnames(affected.wide)) %>%
    select(Tumor,Patient)
  
  ngenes <- nrow(affected.wide)
  
  plot.colors <- c(colorPalette, "amp" = "firebrick", "del" = "blue", "up" = NA, "down" = NA,
                   "splicing" = "forestgreen", "germline" = "purple", "somatic" = "#36454F")
  
  p <- oncoprint(affected.wide) + 
    geom_tile(data=patients, aes(x=Patient,y=ngenes+1,fill=Tumor), height=0.3) +
    scale_fill_manual(values = plot.colors) +
    smartas_theme() +
    theme(axis.text.x=element_blank())
  
  if (sampling.pval < 0.05){
    sampling.pval <- format(sampling.pval,scientific=TRUE)
    q <- p + labs(title=paste(gene,"p <= ",sampling.pval), x="")
    ggsave(paste0(wdp,"ppi/resampling.",gene,".png"),q,  width=10, height=10)
  }
  
  if (fisher.pval < 0.05){
    fisher.pval <- format(fisher.pval,scientific=TRUE)
    q <- p + labs(title=paste(gene,"p <= ",fisher.pval), x="")
    ggsave(paste0(wdp,"ppi/fisher.",gene,".png"), q, width = 10, height = 10)
  }
}

result <- data.frame(Gene=gene,fisher.p=fisher.pval,sampling.p=sampling.pval,
                     n.patients=sum(colSums(mutmat) > 0), n.genes=nrow(mutmat), 
                     n.drivers=sum(rownames(mutmat) %in% drivers$geneHGNCsymbol))
write_tsv(result,paste0(wdp,"ppi/",gene,".txt"))