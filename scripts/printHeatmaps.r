#!/soft/R/R-3.0.0/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
tag <- args[1]
load(args[2])

library(gplots)
library(RColorBrewer)
library(plyr)

top <- 40

topCandidates <- ddply(candidateList,.(Genename), summarise, Replicated=sum(Replicated))
topCandidates <- topCandidates[with(topCandidates, order(-Replicated)), ]
topCandidates <- head(topCandidates, n=top)

fig <- data.frame(matrix(nrow=top, ncol=inputData[["Replicates"]]))
rownames(fig) <- topCandidates$Genename
colnames(fig) <- seq(1,inputData[["Replicates"]])

for (replicate in seq(1,inputData[["Replicates"]])){
  for (gene in topCandidates$Genename){
    if (gene %in% candidates[[replicate]]$Genename) {
      fig[gene, replicate] <- head(candidates[[replicate]]$Switch[candidates[[replicate]]$Genename == gene],1)
    } else {
      fig[gene, replicate] <- NA
    }
  }
}

png(paste0(tag, "_heatmap.png"), width=960, height=960)
myPalette <- colorRampPalette(c("white", "firebrick2"))(n = 14)
heatmap.2(as.matrix(fig), trace="none", scale="none", col=myPalette, na.col="grey", Rowv=NULL, Colv=NULL, 
          dendrogram="none", breaks=seq(0, max(fig, na.rm=T), length.out=15), main="PSI Switch")
graphics.off()

save(fig, tag, file=paste0(tag, ".RData"))