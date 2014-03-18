#!/soft/R/R-3.0.0/bin/Rscript

kk <- vector('list', length(unique(intraReplicate[[8]]$Gene)))
counter <- 1
for(aRow in unique(intraReplicate[[8]]$Gene)){
	tTPM <- intraReplicate[[8]]$la_tTPM[intraReplicate[[8]]$Gene == aRow][1]
	MAD <- interReplicate_N$MAD[interReplicate_N$Gene == aRow][1]
	kk[[counter]] <- data.frame(MAD=MAD,tTPM=tTPM)
	counter <- counter + 1
}
df <- do.call("rbind", kk)
View(df)
qplot(df$MAD, df$tTPM)