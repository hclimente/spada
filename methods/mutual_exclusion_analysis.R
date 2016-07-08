args <- commandArgs(trailingOnly = TRUE)
wes.in <- args[1]
wgs.in <- args[2]
out.file <- args[3]

wes <- read.delim(wes.in)
wes <- wes[,c("Tumor","GeneId","Symbol","Normal_transcript","Tumor_transcript","MS","M","N","p")]
wes$p.me <- wes$p
wgs <- read.delim(wgs.in)

meScore <- merge(wgs,wes,all.x=TRUE,by=c("GeneId","Symbol","Normal_transcript","Tumor_transcript"))

meScore$score <- -log10(meScore$p.me) + log10(meScore$p.o)
meScore <- meScore[order(-meScore$score),]

top <- quantile(meScore$score,0.99,na.rm=T)
meScore$candidate <- 0
meScore$candidate[meScore$score > top] <- 1

write.table(meScore[,c("GeneId","Symbol","Normal_transcript","Tumor_transcript","candidate")], out.file, sep="\t", row.names=F, quote=F)