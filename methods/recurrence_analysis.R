source("~/smartas/pipeline/scripts/variablesAndFunctions.r")

args <- commandArgs(trailingOnly = TRUE)
switches.file <- args[1]
switches.origin <- args[2]
out.file <- args[3]

# read switches
origin <- read_tsv(switches.origin)

switches <- read_tsv(switches.file) %>%
  merge(origin) %>%
  filter(Origin == "Tumor") %>%
  ## calculate patient number
  mutate(Patient_number = Patients_affected %>% strsplit(",") %>% lapply(length) %>% unlist)

# calculate expected frequency of a switch
patientNumber <- switches$Patients_affected %>% strsplit(",") %>% unlist %>% unique %>% length
geneNumber <- switches$GeneId %>% unique %>% length
f.exp <- sum(switches$Patient_number)/(geneNumber*patientNumber)

# binomial test
tests <- lapply(switches$Patient_number, binom.test, patientNumber, f.exp)
switches %>%
  mutate(p.recurrence = unlist(lapply(tests,function(x){x$p.value})),
         padj.recurrence = p.adjust(p.recurrence),
         what = ifelse(switches$Patient_number/patientNumber > f.exp, "greater", "less")) %>%
  # save results
  select(GeneId,Symbol,Normal_transcript,Tumor_transcript,p.recurrence,padj.recurrence,what) %>%
  write_tsv(out.file)
