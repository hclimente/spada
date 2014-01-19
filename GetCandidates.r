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

for (kmer in inputData[["K-mer"]]){
	for (replicate in inputData[["Replicates"]]){

		tag <- paste0(inputData[["Compartments"]][1], replicate, "_", kmer)
		candidates[[tag]] <- data.frame(Gene=as.character(), Entropy_Ref=as.numeric(), Entropy_Alt=as.numeric(), Switch=as.numeric(), maxdPSI=as.character(), mindPSI=as.character())

		#Filter by deltaPSI and expression, based on the FPR
		psiThreshold <- abs(intraReplicate[[tag]]$deltaPSI) > minPSI
		expressionThreshold <- intraReplicate[[tag]]$la_tTPM > minCandidateExpression
	
		for (aCandidate in unique(intraReplicate[[tag]]$Gene[psiThreshold & expressionThreshold])){

			thisGeneData <- intraReplicate[[tag]]$Gene == aCandidate
			posPSI <- intraReplicate[[tag]]$deltaPSI > 0
			negPSI <- intraReplicate[[tag]]$deltaPSI < 0

			#Calculate max difference between the transcripts
			maxDeltaPsi <- max(intraReplicate[[tag]]$deltaPSI[thisGeneData & posPSI], na.rm=T)
			minDeltaPsi <- min(intraReplicate[[tag]]$deltaPSI[thisGeneData & negPSI], na.rm=T)
			maxSwitch <- maxDeltaPsi - minDeltaPsi

			maxDeltaPsiCond <- intraReplicate[[tag]]$deltaPSI == maxDeltaPsi
			minDeltaPsiCond <- intraReplicate[[tag]]$deltaPSI == minDeltaPsi

			candidates[[tag]] <- rbind( candidates[[tag]], data.frame(Gene=aCandidate, Entropy_Ref=calculateEntropy(intraReplicate[[tag]]$PSI_ref[thisGeneData]), 
										Entropy_Alt=calculateEntropy(intraReplicate[[tag]]$PSI_alt[thisGeneData]), Switch=maxSwitch, 
										maxdPSI=intraReplicate[[tag]]$Transcript[thisGeneData & maxDeltaPsiCond], 
										mindPSI=intraReplicate[[tag]]$Transcript[thisGeneData & minDeltaPsiCond])	)
		}

		expressedGenes <- intraReplicate[[tag]]$la_tTPM >= minExpression
		allGenes <- unique(c(allGenes, intraReplicate[[tag]]$Transcript[expressedGenes]))

	}
}

candidateList <- data.frame(maxdPSI=as.character(), mindPSI=as.character())

for (aCondition in candidates){

  entropyCutRef <- aCondition$Entropy_Ref < median(aCondition$Entropy_Ref)  
  entropyCutAlt <- aCondition$Entropy_Alt < median(aCondition$Entropy_Alt)
  switchCut <- aCondition$Switch > median(aCondition$Switch)
  
  candidateList <- rbind(candidateList, aCondition[entropyCutRef & entropyCutAlt & switchCut, c("maxdPSI","mindPSI", "Gene")])
  
}

write.table(candidateList, file=paste0(wd, "/Results/candidateList.lst"), sep="\t", row.names=F, col.names=F, quote=F)
write(allGenes, paste0(wd, "/Results/expressedGenes.lst"), sep="\n")

save(isoformExpression, intraReplicate, interReplicate, candidates, candidateList, inputData, wd, file="SmartAS.RData")
