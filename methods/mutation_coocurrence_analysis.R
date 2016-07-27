source("~/smartas/pipeline/scripts/variablesAndFunctions.r")

args <- commandArgs(trailingOnly = TRUE)
switches.in <- args[1]
wgs.in <- args[2]
out.candidates <- args[3]

#wgs.in <- "/home/hector/smartas/analyses/luad/mutations/wgs_mutations.txt"
#switches.in <- "/home/hector/smartas/analyses/luad/candidateList_info.tsv"
#out.candidates <- "/home/hector/smartas/analyses/luad/candidateList_mutationCoocurrence.tsv"

# read WGS mutations
wgs <- read_tsv(wgs.in) %>%
	set_colnames(c("Tumor","GeneId","Symbol","Patient","Position","Reference","Variant")) %>%
	mutate(Patient = paste0(Patient,"T"))

wgs.patients <- unique(wgs[,c("GeneId","Symbol","Patient")]) %>%
	set_colnames(c("GeneId","Symbol","MutPatient")) %>%
	mutate(Patient = MutPatient)

# get patients affected by any switch
## read switches
switches.info <- read_tsv(switches.in)

## convert to same format as mutations
nocol <- switches.info$Patients_affected %>%
  strsplit(",") %>%
  lapply(length) %>%
  unlist %>%
  max

switch.patients <- switches.info %>%
  select(GeneId,Symbol,Normal_transcript,Tumor_transcript,Patients_affected) %>%
  separate(Patients_affected, into=paste0("Patient_",1:nocol)) %>%
  melt(measure.vars = paste0("Patient_",1:nocol)) %>%
  select(-variable) %>%
  set_colnames(c("GeneId","Symbol","Normal_transcript","Tumor_transcript","SwitchedPatient")) %>%
  mutate(Patient = SwitchedPatient) %>%
  filter(!is.na(SwitchedPatient))

# get overlap between studied mechanisms
wgs.patients <- subset(wgs.patients, MutPatient %in% switch.patients$SwitchedPatient)
switch.patients <- subset(switch.patients, SwitchedPatient %in% wgs.patients$MutPatient)

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

discount <- patients %>%
	group_by(GeneId,Symbol) %>%
	summarise(anySwitch = sum(!is.na(SwitchedPatient)))

totalPatients <- length(unique(patients$Patient))
co <- merge(co,discount,all.x=TRUE) %>%
	mutate(N = totalPatients - (MS+M+S) - (anySwitch - (MS + S)))

co <- co %>%
  select(GeneId,Symbol,Normal_transcript,Tumor_transcript,MS,M,S,N)

## fisher test
co <- apply(co[,c("MS","M","S","N")],1, function(x){
  y <- x + 0.5  
  f <- fisher.test(x=matrix(x,nrow=2,ncol=2),alternative="greater")
  or <- y[1]*y[4]/(y[2]*y[3])
  
  c(p.o=f$p.value,or=f$estimate,my.or=or)}) %>%
	t %>%
	as.data.frame %>%
	set_colnames(c("p.o","OR","eOR")) %>%
  cbind(co,.)

# complete information and save
co.info <- switches.info %>%
	select(GeneId,Symbol,Normal_transcript,Tumor_transcript,Annotation,DriverAnnotation,Driver,Druggable) %>%
	merge(co,all.x=TRUE) %>%
	arrange(p.o) %>%
  mutate(Tumor=unique(wgs$Tumor)) %>%
	select(Tumor,GeneId,Symbol,Normal_transcript,Tumor_transcript,MS,M,S,N,p.o)

write.table(co.info, out.candidates, sep="\t", row.names=F, quote=F)