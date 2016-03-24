args <- commandArgs(trailingOnly = TRUE)
switches.file <- args[1]
out.file <- args[2]

# read switches
switches <- read.delim(switches.file)

## calculate patient number
switches$Patient_number <- unlist(lapply(strsplit(as.character(switches$Patients_affected),","),length))

# calculate expected frequency of a switch
patientNumber <- length(unique(unlist(strsplit(as.character(switches$Patients_affected),","))))
geneNumber <- length(unique(switches$GeneId))
p <- sum(switches$Patient_number)/(geneNumber*patientNumber)

# binomial test
tests <- lapply(switches$Patient_number, binom.test, patientNumber, p, "greater")
switches$p.recurrence <- unlist(lapply(tests,function(x){x$p.value}))
switches$padj.recurrence <- p.adjust(switches$p.recurrence)

# save results
write.table(switches[,c("GeneId","Symbol","Normal_transcript","Tumor_transcript","p.recurrence","padj.recurrence")], out.file, sep="\t", row.names=F, quote=F)