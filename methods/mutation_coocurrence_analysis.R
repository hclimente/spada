library(plyr)

args <- commandArgs(trailingOnly = TRUE)
switches.in <- args[1]
wgs.in <- args[2]
out.candidates <- args[3]

#wgs.in <- "/home/hector/smartas/analyses/luad/mutations/wgs_mutations.txt"
#switches.in <- "/home/hector/smartas/analyses/luad/candidateList_info.tsv"
#out.candidates <- "/home/hector/smartas/analyses/luad/candidateList_mutationCoocurrence.tsv"

# read WGS mutations
wgs <- read.delim(wgs.in)
colnames(wgs) <- c("Tumor","GeneId","Symbol","Patient","Position","Reference","Variant")
wgs$Patient <- paste0(wgs$Patient,"T")

wgs.patients <- unique(wgs[,c("GeneId","Symbol","Patient")])
colnames(wgs.patients) <- c("GeneId","Symbol","MutPatient")
wgs.patients$Patient <- wgs.patients$MutPatient

# get patients affected by any switch
## read switches
switches.info <- read.delim(switches.in)

## convert to same format as mutations
switch.patients <- ddply(switches.info,.(GeneId,Symbol,Normal_transcript,Tumor_transcript),summarise,PatientAffected=unlist(strsplit(paste(Patients_affected,collapse=","),",")))
colnames(switch.patients) <- c("GeneId","Symbol","Normal_transcript","Tumor_transcript","SwitchedPatient")
switch.patients$Patient <- switch.patients$SwitchedPatient

# get overlap between studied mechanisms
wgs.patients <- subset(wgs.patients,MutPatient %in% switch.patients$SwitchedPatient)
switch.patients <- subset(switch.patients,SwitchedPatient %in% wgs.patients$MutPatient)

# calculate unbalance
## get contingency table
patients <- merge(wgs.patients,switch.patients,all=TRUE)
co.tmp <- ddply(patients,.(GeneId,Symbol,Normal_transcript,Tumor_transcript),summarise,
            MS=sum(!is.na(SwitchedPatient) & !is.na(MutPatient)),
            M=sum(is.na(SwitchedPatient) & !is.na(MutPatient)),
            S=sum(!is.na(SwitchedPatient) & is.na(MutPatient)))

### get M cases (where no pair of transcripts can be inferred)
co.s <- subset(co.tmp, !is.na(Normal_transcript) & !is.na(Tumor_transcript))
co.m <- subset(co.tmp, is.na(Normal_transcript) & is.na(Tumor_transcript),select=c("GeneId","M"))

co <- merge(co.m,co.s,by=c("GeneId"))
co$M <- co$M.x

discount <- ddply(patients,.(GeneId,Symbol),summarise,
                  anySwitch=sum(!is.na(SwitchedPatient)))

totalPatients <- length(unique(patients$Patient))
co <- merge(co,discount,all.x=TRUE)
co$N <- totalPatients - (co$MS+co$M+co$S) - (co$anySwitch - (co$MS+co$S))

co <- co[,c("GeneId","Symbol","Normal_transcript","Tumor_transcript","MS","M","S","N")]

## fisher test
f <- apply(co[,c("MS","M","S","N")],1, function(x){
  y <- x + 0.5  
  f <- fisher.test(x=matrix(x,nrow=2,ncol=2),alternative="greater")
  or <- y[1]*y[4]/(y[2]*y[3])
  
  c(p.o=f$p.value,or=f$estimate,my.or=or)})

f.df <- as.data.frame(t(f))
colnames(f.df) <- c("p.o","OR","eOR")

co <- cbind(co,f.df)

# complete information and save
co.info <- merge(co,switches.info[,c("GeneId","Symbol","Normal_transcript","Tumor_transcript","Annotation","DriverAnnotation","Driver","Druggable")],all.y=TRUE)
co.info <- co.info[order(co.info$p.o),]
co.info <- co.info[,c("GeneId","Symbol","Normal_transcript","Tumor_transcript","MS","M","S","N","p.o")]

write.table(co.info, out.candidates, sep="\t", row.names=F, quote=F)