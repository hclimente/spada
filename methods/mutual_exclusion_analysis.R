source("~/smartas/pipeline/scripts/variablesAndFunctions.r")

args <- commandArgs(trailingOnly = TRUE)
wes.in <- args[1]
wgs.in <- args[2]
out.file <- args[3]

wes <- read_tsv(wes.in) %>%
  select(Tumor,GeneId,Symbol,Normal_transcript,Tumor_transcript,MS,M,N,p) %>%
  mutate( p.me = p)
wgs <- read_tsv(wgs.in)

meScore <- merge(wgs, wes, all.x=TRUE,
                 by=c("Tumor","GeneId","Symbol","Normal_transcript","Tumor_transcript")) %>%
  mutate(score = -log10(p.me) + log10(p.o)) %>%
  arrange(score)

top <- quantile(meScore$score,0.99,na.rm=T)

meScore <- mutate(meScore, candidate = ifelse( meScore$score > top, 1, 0))

meScore %>%
  select(Tumor,GeneId,Symbol,Normal_transcript,Tumor_transcript,candidate) %>%
  write_tsv(out.file)