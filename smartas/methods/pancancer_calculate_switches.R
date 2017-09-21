tryCatch(source("~/smartas/pipeline/scripts/variablesAndFunctions.r"),error=function(e){})

args <- commandArgs(trailingOnly = TRUE)
switches.filename <- args[1]

all.switches <- list()
for (cancer in cancerTypes){
  switches.file <- paste0("~/smartas/analyses/",cancer,"/",switches.filename,".tsv")
  all.switches[[cancer]] <- read_tsv(switches.file) %>%
    mutate(PatientNumber = unlist(lapply(strsplit(Patients_affected,",",fixed=T),length)),
           Percentage = PatientNumber/nPatients[[cancer]], Tumor = cancer)
}

all.switches <- do.call("rbind",all.switches)

all.switches.file <- paste0("~/smartas/analyses/pancancer/",switches.filename,".tumorSplit.tsv")
write.table(all.switches, all.switches.file, sep="\t", row.names=F, quote=F)

# a switch is regarded as Reliable pancancer if it's NotNoise and Model in at least a tumor
switches.agg <- all.switches %>%
    group_by(GeneId,Symbol,Normal_transcript,Tumor_transcript,Normal_protein,
             Tumor_protein,DriverAnnotation,Driver,Druggable,IsFunctional,
             CDS_Normal,CDS_Tumor,CDS_change,UTR_change) %>%
    summarise(Reliable=max((NotNoise+IsModel)==2), 
              Tumors=paste(Tumor,collapse = ","),
              Patients_affected=paste(Patients_affected,collapse = ","), 
              PatientNumber=sum(PatientNumber), 
              Percentage=sum(PatientNumber)/nPatients[["total"]]) %>%
    arrange(desc(Percentage)) %>%
    select(GeneId,Symbol,Normal_transcript,Tumor_transcript,
           Normal_protein,Tumor_protein,DriverAnnotation,
           Reliable,IsFunctional,Driver,Druggable,CDS_Normal,
           CDS_Tumor,CDS_change,UTR_change,Tumors,PatientNumber,
           Percentage,Patients_affected)

agg.switches.file <- paste0("~/smartas/analyses/pancancer/",switches.filename,".agg.tsv")
write.table(switches.agg,agg.switches.file,quote=F,row.names=F, sep="\t")
