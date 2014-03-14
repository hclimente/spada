load("~/SmartAS/Results/TCGA/luad_mE0.0_mCE-1.0/RWorkspaces/2_GetCandidates.RData")

library(ggplot2)

genesOfInterest <- c("NUMB", "QKI", as.vector(head(candidateList$Genename, n=30)) )
thisKansur <- "luad"

for (gene in genesOfInterest){
  
  cat(gene,"\n")
  
	mindPSI <- as.character(candidateList$mindPSI[ candidateList$Genename == gene ][1])
	maxdPSI <- as.character(candidateList$maxdPSI[ candidateList$Genename == gene ][1])
  
	ni_dist_N <- as.numeric()
	ni_dist_T <- as.numeric()
	tumCand <- as.numeric()
	norCand <- as.numeric()
  
	for (rep in seq(1, inputData[["Replicates"]])){
	  if(gene %in% candidates[[rep]]$Genename){
	    tumCand <- c(tumCand, intraReplicate[[rep]]$PSI_T[intraReplicate[[rep]]$Transcript == mindPSI])
	    norCand <- c(norCand, intraReplicate[[rep]]$PSI_T[intraReplicate[[rep]]$Transcript == maxdPSI])
	  }
	  
	  ni_dist_N <- c(ni_dist_N, intraReplicate[[rep]]$PSI_N[intraReplicate[[rep]]$Transcript == maxdPSI])
	  ni_dist_T <- c(ni_dist_T, intraReplicate[[rep]]$PSI_N[intraReplicate[[rep]]$Transcript == mindPSI])
	}

	normalDist <- data.frame(PSI=ni_dist_N, Condition="Normal")
  normalDist <- rbind(normalDist, data.frame(PSI=norCand, Condition="Tumor"))
  
  tumorDist <- data.frame(PSI=ni_dist_T, Condition="Normal")
  tumorDist <- rbind(tumorDist, data.frame(PSI=tumCand, Condition="Tumor"))

	cat("Normal isoform p-values\n")
	print(ecdf(normalDist$PSI[normalDist$Condition=="Normal"])(normalDist$PSI[normalDist$Condition=="Tumor"]))

	cat("Tumor isoform p-values\n")
  print(ecdf(tumorDist$PSI[tumorDist$Condition=="Normal"])(tumorDist$PSI[tumorDist$Condition=="Tumor"]))

	png(paste0("~/", thisKansur, "_", gene, "_isoT(",mindPSI ,").png"), width=1080, height=1080)
  print(ggplot(tumorDist, aes(x=PSI, fill=Condition)) + geom_density(alpha=.3) + scale_fill_manual( values = c("steelblue3","firebrick2")) + scale_x_continuous(limits = c(0, 1)) + theme_minimal(base_size=40) + guides(fill=FALSE))
  #print(ggplot(NULL) + geom_histogram(data=normalDist, aes(x=PSI_T, fill = ..count..), binwidth = 0.01) + scale_fill_gradient("Count", low = "black", high = "red") + geom_point(data=ti_candidates, aes(x=PSI, y=Y, size=5)))
	dev.off()

	png(paste0("~/", thisKansur, "_", gene, "_isoN(",maxdPSI ,").png"), width=1080, height=1080)
  print(ggplot(normalDist, aes(x=PSI, fill=Condition)) + geom_density(alpha=.3) + scale_fill_manual( values = c("steelblue3","firebrick2")) + scale_x_continuous(limits = c(0, 1)) + theme_minimal(base_size=40) + guides(fill=FALSE))
  #print(ggplot(NULL) + geom_histogram(data=normalDist, aes(x=PSI_N, fill = ..count..), binwidth = 0.01) + scale_fill_gradient("Count", low = "black", high = "red") + geom_point(data=ni_candidates, aes(x=PSI, y=Y, size=5)))
	dev.off()
}

ggplot(normalDist, aes(x=PSI_N, fill=condition)) + geom_histogram(binwidth=.01, alpha=.5, position="identity")
ggplot(tumorDist, aes(x=PSI_T, fill=condition)) + geom_histogram(binwidth=.01, alpha=.5, position="identity")
ggplot(tumorDist, aes(x=PSI_T, fill=condition)) + geom_density(alpha=.3)

ggplot(normalDist, aes(x=PSI_N, fill=condition)) + geom_density(alpha=.3) + scale_fill_manual( values = c("steelblue3","firebrick2"))