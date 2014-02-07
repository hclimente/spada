#!/soft/R/R-3.0.0/bin/Rscript

load("SmartAS.RData")
setwd(wd)

calculateEntropy <- function(x){
	H <- 0
	for (p in x){
		if(p == 0)
			next
		H <- H + p * log(p, base=10)
	}
	
	maxH <- log(length(x[x!=0]), base=10)
	if( maxH == 0 & length(x) > 1 ) {
		return(0)
	} else if(maxH == 0){
		return(9999)
	} else {
		return(-H/maxH)
	}
}

candidates <- list()
allGenes <- as.character()

args <- commandArgs(trailingOnly = TRUE)
minExpression <- as.numeric(args[1])
minCandidateExpression <- as.numeric(args[2])
minPSI <- as.numeric(args[3])

for (replicate in seq(1,inputData[["Replicates"]])){

	candidates[[replicate]] <- data.frame(Gene=as.character(), Genename=as.character(), Entropy_Ref=as.numeric(), Entropy_Alt=as.numeric(), 
										  Switch=as.numeric(), maxdPSI=as.character(), mindPSI=as.character())

	#Filter by deltaPSI and expression, based on the FPR
	psiThreshold <- abs(intraReplicate[[replicate]]$deltaPSI) > minPSI
	expressionThreshold <- intraReplicate[[replicate]]$la_tTPM > minCandidateExpression
			
	for (aCandidate in unique(intraReplicate[[replicate]]$Gene[psiThreshold & expressionThreshold])){

		thisGeneData <- intraReplicate[[replicate]]$Gene == aCandidate
		posPSI <- intraReplicate[[replicate]]$deltaPSI > 0
		negPSI <- intraReplicate[[replicate]]$deltaPSI < 0

		#Calculate max difference between the transcripts
		maxDeltaPsi <- max(intraReplicate[[replicate]]$deltaPSI[thisGeneData & posPSI], na.rm=T)
		minDeltaPsi <- min(intraReplicate[[replicate]]$deltaPSI[thisGeneData & negPSI], na.rm=T)
		maxSwitch <- maxDeltaPsi - minDeltaPsi

		#deltaPSI = PSI_ref - PSI_alt
		#	MaxDeltaPsiCond: predominant transcript in the normal
		#	MinDeltaPsiCond: predominant transcript in the tumor
		maxDeltaPsiCond <- intraReplicate[[replicate]]$deltaPSI == maxDeltaPsi
		minDeltaPsiCond <- intraReplicate[[replicate]]$deltaPSI == minDeltaPsi

		candidates[[replicate]] <- rbind( candidates[[replicate]], data.frame(Gene=aCandidate, Entropy_Ref=calculateEntropy(intraReplicate[[replicate]]$PSI_ref[thisGeneData]), 
										  Entropy_Alt=calculateEntropy(intraReplicate[[replicate]]$PSI_alt[thisGeneData]), Switch=maxSwitch, 
										  maxdPSI=intraReplicate[[replicate]]$Transcript[thisGeneData & maxDeltaPsiCond], 
										  Genename=intraReplicate[[replicate]]$Genename[thisGeneData & maxDeltaPsiCond], 
										  mindPSI=intraReplicate[[replicate]]$Transcript[thisGeneData & minDeltaPsiCond])
										)
	}

	#Expressed genes: transcript whose expression is above the threshold
	expressedGenes <-(log(intraReplicate[[replicate]]$TPM_ref) + log(intraReplicate[[replicate]]$TPM_alt)) / 2 > minExpression
	allGenes <- unique(c(allGenes, intraReplicate[[replicate]]$Transcript[expressedGenes]))

	entropyCutRef <- candidates[[replicate]]$Entropy_Ref < median(candidates[[replicate]]$Entropy_Ref)  
 	entropyCutAlt <- candidates[[replicate]]$Entropy_Alt < median(candidates[[replicate]]$Entropy_Alt)
	switchCut <- candidates[[replicate]]$Switch > median(candidates[[replicate]]$Switch)

	candidates[[replicate]] <- candidates[[replicate]] [entropyCutRef & entropyCutAlt & switchCut, c("Genename", "Gene", "maxdPSI","mindPSI")])

}

candidateList <- data.frame(maxdPSI=as.character(), mindPSI=as.character())

for (aCondition in candidates){
	candidateList <- rbind(candidateList, aCondition[entropyCutRef & entropyCutAlt & switchCut, c("Genename", "Gene", "maxdPSI","mindPSI")])
}

candidateList <- unique(candidateList)

write.table(candidateList, file=paste0(wd, "/Results/candidateList.tsv"), sep="\t", row.names=F, col.names=F, quote=F)
write(allGenes, paste0(wd, "/Results/expressedGenes.lst"), sep="\n")

save(isoformExpression, intraReplicate, interReplicate, candidates, candidateList, inputData, wd, file="SmartAS.RData")
