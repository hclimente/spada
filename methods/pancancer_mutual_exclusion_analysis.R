source("~/smartas/pipeline/scripts/variablesAndFunctions.r")

library(plyr)

out.file <- "~/smartas/analyses/pancancer/candidateList_mutationME.tsv"
wgs.agg.out <- "~/smartas/analyses/pancancer/mutations/gene_wgs_mutations_all_switches.txt"
wes.agg.out <- "~/smartas/analyses/pancancer/mutations/gene_functional_mutations_all_switches.txt"

wgs <- list()
for (cancer in setdiff(cancerTypes,c("coad","kirp","lihc"))){
  
  wgs.cancer.file <- paste0("~/smartas/analyses/",cancer,"/mutations/gene_wgs_mutations_all_switches.txt")
  wgs.cancer <- read.delim(wgs.cancer.file)
  wgs.cancer$Tumor <- cancer
  wgs.cancer$GeneId <- as.character(wgs.cancer$GeneId)
  wgs.cancer$Symbol <- as.character(wgs.cancer$Symbol)
  wgs.cancer$Normal_transcript <- as.character(wgs.cancer$Normal_transcript)
  wgs.cancer$Tumor_transcript <- as.character(wgs.cancer$Tumor_transcript)
  
  wgs[[cancer]] <- wgs.cancer
  
}

wgs <- do.call("rbind",wgs)
wgs <- wgs[!as.logical(rowSums(is.na(wgs[,c("MS","M","S","N")]))),]

wgs.agg <- ddply(wgs,.(GeneId,Symbol,Normal_transcript,Tumor_transcript),summarise,
                 MS=sum(MS),M=sum(M),S=sum(S),N=sum(N))
wgs.agg$p.o <- apply(wgs.agg[,c("MS","M","S","N")],1,function(x){
  fisher.test(matrix(x,2,2),alternative="greater")$p.value})
write.table(wgs.agg, wgs.agg.out, sep="\t", row.names=F, quote=F)

wes <- list()
for (cancer in cancerTypes){
  
  wes.cancer.file <- paste0("~/smartas/analyses/",cancer,"/mutations/gene_functional_mutations_all_switches.txt")
  wes.cancer <- read.delim(wes.cancer.file)
  wes.cancer$Tumor <- as.character(wes.cancer$Tumor)
  wes.cancer$GeneId <- as.character(wes.cancer$GeneId)
  wes.cancer$Symbol <- as.character(wes.cancer$Symbol)
  wes.cancer$Normal_transcript <- as.character(wes.cancer$Normal_transcript)
  wes.cancer$Tumor_transcript <- as.character(wes.cancer$Tumor_transcript)
  
  wes[[cancer]] <- wes.cancer 
}

wes <- do.call("rbind",wes)

totals <- ddply(subset(wes, Normal_transcript=="None"),.(GeneId,Symbol),summarise,
                M.all=sum(M),All=sum(N,M))
wes.agg <- ddply(wes,.(GeneId,Symbol,Normal_transcript,Tumor_transcript),
                 summarise,MS=sum(MS),S=sum(S),Switched=paste(Switched,collapse=","),
                 Mutated=paste(Mutated,collapse=","))

wes.agg <- merge(wes.agg,totals)
wes.agg$M <- wes.agg$M.all - wes.agg$MS
wes.agg$N <- wes.agg$All - (wes.agg$MS+wes.agg$S+wes.agg$M)

# remove those genes with dummy switch (None,None) that have other swith described
wes.agg <- wes.agg[wes.agg$Tumor_transcript!="None" | (wes.agg$Tumor_transcript=="None" & !(wes.agg$GeneId %in% wes.agg$GeneId[wes.agg$Tumor_transcript!="None"])),]

wes.agg$p.me <- apply(wes.agg[,c("MS","M","S","N")],1,function(x){
  fisher.test(matrix(x,2,2),alternative="less")$p.value})

wes.agg <- wes.agg[,c("GeneId","Symbol","Normal_transcript","Tumor_transcript","MS","M","S","N","p.me")]
write.table(wes.agg, wes.agg.out, sep="\t", row.names=F, quote=F)

meScore <- merge(wes.agg,wgs.agg,by=c("GeneId","Symbol","Normal_transcript","Tumor_transcript"))
meScore$score <- -log10(meScore$p.me) + log10(meScore$p.o)

meScore <- meScore[order(-meScore$score),]

top <- quantile(meScore$score,0.99,na.rm=T)
meScore$candidate <- 0
meScore$candidate[meScore$score > top] <- 1

ggplot(wes.agg,aes(x=log2(M),y=log2(S),color=(M+S)/(M+S+MS))) + geom_point() + geom_text(aes(label=Symbol))

write.table(meScore[,c("GeneId","Symbol","Normal_transcript","Tumor_transcript","candidate")], out.file, sep="\t", row.names=F, quote=F)