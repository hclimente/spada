library(plyr)
library(tidyr)

source("~/smartas/pipeline/scripts/variablesAndFunctions.r")
out.co <- "~/smartas/analyses/pancancer/candidateList_mutationCoocurrence.tsv"

coocurrence <- list()
for (cancer in setdiff(cancerTypes,c("coad","kirp","lihc"))){
  wgs.in <- paste0("~/smartas/analyses/",cancer,"/mutations/wgs_mutations.txt")
  switches.in <- paste0("~/smartas/analyses/",cancer,"/candidateList_info.tsv")
  
  # read WGS mutations
  wgs <- read_tsv(wgs.in) %>%
    set_colnames(c("Tumor","GeneId","Symbol","Patient","Position","Reference","Variant")) %>%
    mutate(Patient = paste0(Patient,"T"))
  
  wgs.patients <- wgs %>%
    select(GeneId,Symbol,Patient) %>%
    unique %>%
    set_colnames(c("GeneId","Symbol","MutPatient")) %>%
    mutate(Patient = MutPatient)
  
  # get patients affected by any switch
  ## read switches
  switches.info <- read_tsv(switches.in)
  nocols <- max(unlist(lapply(strsplit(switches.info$Patients_affected,","),length)))
  
  ## convert to same format as mutations
  switch.patients <- separate(switches.info, Patients_affected, paste("Patient", 1:nocols, sep="_"), sep=",") %>%
    melt(id.vars = c("GeneId","Symbol","Normal_transcript","Tumor_transcript"), measure.vars = paste("Patient", 1:nocols, sep="_")) %>%
    filter(!is.na(value)) %>%
    select(-variable) %>%
    set_colnames(c("GeneId","Symbol","Normal_transcript","Tumor_transcript","SwitchedPatient")) %>%
    mutate(Patient = SwitchedPatient)
  
  # get overlap between studied mechanisms
  wgs.patients <- subset(wgs.patients,MutPatient %in% switch.patients$SwitchedPatient)
  switch.patients <- subset(switch.patients,SwitchedPatient %in% wgs.patients$MutPatient)
  
  # calculate unbalance
  ## get contingency table
  patients <- merge(wgs.patients,switch.patients,all=TRUE)
  co.tmp <- patients %>%
    group_by(GeneId,Symbol,Normal_transcript,Tumor_transcript) %>%
    summarise(MS=sum(!is.na(SwitchedPatient) & !is.na(MutPatient)),
              M=sum(is.na(SwitchedPatient) & !is.na(MutPatient)),
              S=sum(!is.na(SwitchedPatient) & is.na(MutPatient)))
  
  ### get M cases (where no pair of transcripts can be inferred)
  co.s <- subset(co.tmp, !is.na(Normal_transcript) & !is.na(Tumor_transcript))
  co.m <- subset(co.tmp, is.na(Normal_transcript) & is.na(Tumor_transcript),select=c("GeneId","M"))
  
  co <- merge(co.m,co.s,by=c("GeneId")) %>%
    mutate(M = M.x)
  
  discount <- ddply(patients,.(GeneId,Symbol),summarise,
                    anySwitch=sum(!is.na(SwitchedPatient)))
  
  totalPatients <- length(unique(patients$Patient))
  co <- merge(co,discount,all.x=TRUE) %>%
    mutate(N = totalPatients - (MS+M+S) - (anySwitch - (MS+S))) %>%
    select(GeneId,Symbol,Normal_transcript,Tumor_transcript,MS,M,S,N)
  
  ## fisher test
  f <- apply(co[,c("MS","M","S","N")],1, function(x){
    y <- x + 0.5  
    f <- fisher.test(x=matrix(x,nrow=2,ncol=2),alternative="greater")
    or <- y[1]*y[4]/(y[2]*y[3])
    
    c(p.o=f$p.value,or=f$estimate,my.or=or)})
  
  f.df <- as.data.frame(t(f))
  colnames(f.df) <- c("p.o","OR","eOR")
  
  co <- cbind(co,f.df)
  co$Tumor <- cancer
  coocurrence[[cancer]] <- co
}

coocurrence <- do.call("rbind",coocurrence)

co.agg <- ddply(coocurrence,.(GeneId,Symbol,Normal_transcript,Tumor_transcript),summarise,
                MS=sum(MS),M=sum(M),S=sum(S),N=sum(N))

f <- apply(co.agg[,c("MS","M","S","N")],1, function(x){
  y <- x + 0.5  
  f <- fisher.test(x=matrix(x,nrow=2,ncol=2),alternative="greater")
  or <- y[1]*y[4]/(y[2]*y[3])
  
  c(p.o=f$p.value,or=f$estimate,my.or=or)}) %>%
  as.data.frame %>%
  t %>%
  set_colnames(c("p.o","OR","eOR"))

co.agg <- cbind(co.agg,f)

write_tsv(co.agg, out.co)
