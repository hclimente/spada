#gene <- "IGF2BP3"
#tx1 <- "uc003swf.2"
#tx2 <- "uc003swg.2"

gene <- "ARL1"
tx1 <- "uc001tib.2"
tx2 <- "uc001tic.2"

InteraX <- read.delim(paste0("~/SmartAS/Results/TCGA/luad_mE-1.0/iLoops/InteraXChanges_", gene, "_", tx1, "_", tx2, ".tsv"))
kk <- vector("list", length(unique(InteraX$Partner_gene)))

mask1 <- abs(InteraX$dRC) != 9999

for (aGene in unique(InteraX$Partner_gene)){
  mask2 <- InteraX$Partner_gene == aGene
  score <- 0
  if (any(mask1 & mask2)){
    score <- score + max(abs(InteraX$dRC[mask1 & mask2]))/50  
  }
  else{
    score <- score + 0.5
  }
  
  if (any(InteraX$Annotation[mask2] == "Driver")){
    score <- score + 0.5
  }
  kk[[aGene]] <- c(gene, aGene, score)
}

perrote <- do.call("rbind", kk)
write.table(perrote, file=paste0("~/Desktop/Baldo_", gene, ".tsv"), sep="\t", row.names=F, col.names=F, quote=F)