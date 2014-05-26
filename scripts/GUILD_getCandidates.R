library(plyr)
Baldo_rawTable <- read.delim("~/Desktop/Baldo_rawTable.tsv")
Baldo_rawTable <- ddply(Baldo_rawTable,.(Gene), summarise, Replicated=sum(Replicates))
Baldo_rawTable <- Baldo_rawTable[with(Baldo_rawTable, order(-Replicated)), ]
Baldo_rawTable$Score <- 0.2 + (Baldo_rawTable$Replicated - 6)/57 * 0.3

View(Baldo_rawTable)
write.table(Baldo_rawTable[,c(1,3)], file="~/Desktop/Baldo_candidates.tsv", sep="\t", row.names=F, col.names=F, quote=F)