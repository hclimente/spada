#!/usr/bin/env Rscript

# get the empirical distribution of the differences between a vector
#   x numerical vector
#   minV minimum value to consider that measure
# Returns:
#   * function returning NA when there are less than 10 valid values
#     (ie non NA and higher than threshold)
#   * function returning 1 to any value higher than 0 when all the
#     differences are 0.
#   * return ecdf otherwise, counting that the minimum p-value will be 1/n+1
getEmpiricalDistribution <- function(x,minV){

  v <- as.numeric(x)

  # discard cases with less than 10 valid cases
  if ( sum(!is.na(v)) < 10 | sum(v[!is.na(v)] >= minV) < 10){
    return(function(x){return(rep(NA,length(x)))})
  } else {
    diffMatrix <- abs(outer(v,v,"-"))
    subtraction <- diffMatrix[lower.tri(diffMatrix, diag = FALSE)]

    if (all(subtraction[!is.na(subtraction)]==0)){
      return(function(x){
        p <- rep(NA,length(x))
        p[x==0] <- 0
        p[x>0] <- 1
        return(p)})
    } else {
      return(ecdf(subtraction))
    }
  }
}

readTxExpression <- function(nTxFile, tTxFile, tx2geneFile) {

  tpmRead <- function(f) {
    read_tsv(f, col_types = cols(transcript = "c", .default = "d")) %>%
	gather(key = sample, value = expr, -1)
  }

  tx2gene <- read_csv(tx2geneFile, col_types = 'cc')

  left_join(tpmRead(tTxFile) %>% rename(tpmTumor = expr),
            tpmRead(nTxFile) %>% rename(tpmNormal = expr),
		    by = c("transcript", "sample")) %>%
    left_join(tx2gene, by = "transcript") %>%
    select(gene, transcript, sample, tpmTumor, tpmNormal) %>%
    arrange(sample, gene, transcript)

}

calculateDPSI <- function(expression) {

  dpsi <- expression %>%
    group_by(gene, sample) %>%
    mutate(psiNormal = tpmNormal/sum(tpmNormal),
           psiTumor = tpmTumor/sum(tpmTumor)) %>%
    ungroup

  normalMedianPsi <- dpsi %>%
    group_by(transcript) %>%
    summarize(psiNormalMedian = median(psiNormal, na.rm = T))

  left_join(dpsi, normalMedianPsi, by = "transcript") %>%
    mutate(psiNormal = ifelse(is.na(psiNormal), psiNormalMedian, psiNormal),
           dPSI = psiNormal - psiTumor) %>%
    select(-psiNormalMedian)

}

scoreDPSI <- function(dpsi) {

  pdPSI <- by(dpsi, dpsi$transcript, function(x) { getEmpiricalDistribution(x$dPSI, 0) })

  dpsi %>%
    group_by(transcript) %>%
    mutate(pdPSI = pdPSI[[unique(transcript)]](dPSI)) %>%
    ungroup

}

calculateDE <- function(normalGeneExpressionFile, tumorGeneExpressionFile) {

	# read expression files
	normalGeneExpression <- read_tsv(normalGeneExpressionFile,
									 col_types = cols(gene = 'c', .default = 'i'))
	normalSamples <- colnames(normalGeneExpression)
 	tumorGeneExpression <- read_tsv(tumorGeneExpressionFile,
	  							    col_types = cols(gene = 'c', .default = 'i'))
	tumorSamples <- colnames(tumorGeneExpression)

	# generate expression matrix
	geneExpression <- left_join(tumorGeneExpression, normalGeneExpression,
	  						  by = 'gene')
	genes <- geneExpression$gene
	geneExpression <- select(geneExpression, -gene) %>% as.matrix
	rownames(geneExpression) <- genes

	# normalize expression and calculate log2
	y <- DGEList(counts = geneExpression)
	y <- calcNormFactors(y)

	logNormGeneExpression <- cpm(y, normalized.lib.sizes=TRUE) %>%
		log2
	logNormGeneExpression[logNormGeneExpression == -Inf] <- NA

	# calculate differential expression
	tibble(gene = rownames(logNormGeneExpression),
		   pDE = apply(logNormGeneExpression, 1, function(x) {
			   wilcox.test(x[normalSamples], x[tumorSamples])$p.value
		   }))

}

calculateSwitches <- function(txInfo) {

  txInfo %>%
    # remove differentially expressed genes
    filter(pDE > 0.01) %>%
    # remove low deltaPSI
    filter(pdPSI < 0.01 & dPSI > 0.05) %>%
    # remove lowly expressed transcript
    filter((dPSI > 0 & tpmTumor > 0.1) | (dPSI < 0 & tpmNormal > 0.1)) %>%
    group_by(gene, sample) %>%
    # remove unelligible genes
    filter(n() > 1 & any(dPSI > 0) & any(dPSI < 0)) %>%
    summarise(normal = transcript[which.min(dPSI)],
              tumor = transcript[which.max(dPSI)]) %>%
    ungroup %>%
    group_by(gene, normal, tumor) %>%
    summarise(samples = paste(samples, collapse = ','))

}

library(tidyverse)
library(edgeR)

# cd /Users/hclimente/projects/spada/spada/tests/data/calculate_switches
# /Users/hclimente/projects/spada/scripts/calculate_switches.R tx_normal tx_tumor gn_normal gn_tumor tx2gene switches.tsv
args <- commandArgs(trailingOnly = TRUE)
nTxFile <- args[1]
tTxFile <- args[2]
nGnFile <- args[3]
tGnFile <- args[4]
tx2geneFile <- args[5]
outfile <- args[6]

print("Differential transcript-expression")
txDPSI <- readTxExpression(nTxFile, tTxFile, tx2geneFile) %>%
	calculateDPSI %>%
	scoreDPSI

print("Calculate differential gene-expression")
gnDiffExpr <- calculateDE(nGnFile, tGnFile)

print("Compute switches")
left_join(txDPSI, gnDiffExpr, by = "gene") %>%
	calculateSwitches %>%
	write_tsv(outfile)
