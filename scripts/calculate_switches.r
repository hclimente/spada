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
  if ( sum(!is.na(v))< 10 | sum(v[!is.na(v)]>=minV)<10){
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

readExpression <- function(nTxFile, tTxFile, tx2geneFile) {
  
  tpmRead <- function(file) {
    read_tsv(file) %>% gather(key = sample, value = expr, -1)
  }
  
  tx2gene <- read_tsv(tx2geneFile)
  
  left_join(tpmRead(tTxFile) %>% rename(tpmTumor = expr),
            tpmRead(nTxFile) %>% rename(tpmNormal = expr)) %>% 
    left_join(tx2gene) %>% 
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
  
  left_join(dpsi, normalMedianPsi) %>%
    mutate(psiNormal = ifelse(is.na(psiNormal), psiNormalMedian, psiNormal),
           deltaPsi = psiNormal - psiTumor) %>%
    select(-psiNormalMedian)
  
}

scoreDPSI <- function(dpsi) {
  
  pDPSI <- by(dpsi, dpsi$transcript, function(x) { getEmpiricalDistribution(x$deltaPsi, 0) })
  
  dpsi %>%
    group_by(transcript) %>%
    mutate(pDPSI = pDPSI[[unique(transcript)]](deltaPsi)) %>%
    ungroup
  
}

calculateDE <- function(normalGeneExpressionFile, tumorGeneExpressionFile) {
  
  tumorGeneExpression <- read.table(tumorGeneExpressionFile, check.names=FALSE)
  tumorSamples <- colnames(tumorGeneExpression)
  normalGeneExpression <- read.table(normalGeneExpressionFile, check.names=FALSE)
  normalSamples <- colnames(normalGeneExpression)
  
  geneExpression <- cbind(tumorGeneExpression, normalGeneExpression)
  y <- DGEList(counts = geneExpression)
  y <- calcNormFactors(y)
  normGeneExpression <- cpm(y, normalized.lib.sizes=TRUE)
  
  logNormGeneExpression <- log2(normGeneExpression)
  logNormGeneExpression[logNormGeneExpression == -Inf] <- NA
  
  tibble(gene = rownames(logNormGeneExpression),
         pDE = apply(logNormGeneExpression, 1, function(x) {
           wilcox.test(x[normalSamples], x[tumorSamples])$p.value
         }))
  
}

calculateSwitches <- function(txInfo) {
  
  txInfo %>%
    # remove differentially expressed genes
    filter(dGE > 0.01)
    # remove low deltaPSI
    filter(pDPSI < 0.01 & deltaPsi > 0.05) %>%
    # remove lowly expressed transcript
    filter((deltaPsi > 0 & tpmTumor > 0.1) | (deltaPsi < 0 & tpmNormal > 0.1)) %>%
    group_by(gene, sample) %>%
    # remove unelligible genes
    filter(n() > 1 & any(deltaPsi > 0) & any(deltaPsi < 0)) %>%
    summarise(normal = transcript[which.min(deltaPsi)],
              tumor = transcript[which.max(deltaPsi)]) %>%
    ungroup %>%
    group_by(gene, normal, tumor) %>%
    summarise(samples = paste(samples, collapse = ','))
  
}

library(tidyverse)
library(edgeR)

setwd("~/projects/spada/spada/tests/data")
nTxFile <- "expression"
tTxFile <- "expression"
tx2geneFile <- "tx2gene"
outfile <- args[4]

# args <- commandArgs(trailingOnly = TRUE)
# nTxFile <- args[1]
# tTxFile <- args[2]
# tx2geneFile <- args[3]
# outfile <- args[4]

txExpression <- readExpression(nTxFile, tTxFile, tx2geneFile)
gnDiffExpr <- calculateDE(nGnFile, tGnFile)
txDPSI <- calculateDPSI(txExpression)
txDPSI <- scoreDPSI(txDPSI)

txInfo <- left_join(txDPSI, gnDiffExpr)

calculateSwitches(txInfo)