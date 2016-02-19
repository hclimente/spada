#!/soft/R/R-3.0.0/bin/Rscript

# get the empirical distribution of the differences between a vector
#   x numerical vector
#   minV minimum value to consider that measure
# Returns:
#   * function returning NA when there are less than 10 valid values 
#     (ie non NA and higher than threshold)
#   * function returning 1 to any value higher than 0 when all the 
#     differences are 0.
#   * return ecdf otherwise
getEmpiricalDistribution <- function(x,minV){
  
  v <- as.numeric(x)
  
  # discard cases with less than 10 valid cases
  if ( sum(!is.na(v))< 10 | sum(v[!is.na(v)]>=minV)<10){
    return(function(x){return(NA)})
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

library(plyr)

args <- commandArgs(trailingOnly = TRUE)
outfile <- args[1]
tumor <- args[2]

# Prepare the data
## read psi
psi.nt <- read.table(paste0("/projects_rg/TCGA/pipeline/run11/",tumor,"_iso_psi_paired-filtered.txt"), check.names=FALSE)
psi.t <- read.table(paste0("/projects_rg/TCGA/pipeline/run11/",tumor,"_iso_psi_tumor-filtered.txt"), check.names=FALSE)
psi <- cbind(psi.nt,psi.t)

rm(psi.nt,psi.t)

## get genes and transcripts
genes <- unlist(strsplit(row.names(psi),","))[c(T,F)]
transcripts <- unlist(strsplit(row.names(psi),","))[c(F,T)]
psi$Gene <- genes

## get patients
pats.n <- grep("^.{4}N$",colnames(psi), value=TRUE)
pats.t <- grep("^.{4}T$",colnames(psi), value=TRUE)
pats.t.nt <- gsub("N$","T",pats.n)
pats.t.t <- setdiff(pats.t,pats.t.nt)

## read isoform expression
xpr.nt <- read.table("/projects_rg/TCGA/pipeline/run11/luad_iso_tpm_paired-filtered.txt", check.names=FALSE)
xpr.t <- read.table("/projects_rg/TCGA/pipeline/run11/luad_iso_tpm_tumor-filtered.txt", check.names=FALSE)
xpr <- cbind(xpr.nt,xpr.t)

rm(xpr.nt,xpr.t)

## read gene expression
xpr.gene.nt <- read.table("/projects_rg/TCGA/pipeline/run11/luad_gene_tpm_paired-filtered.txt", check.names=FALSE)
xpr.gene.t <- read.table("/projects_rg/TCGA/pipeline/run11/luad_gene_tpm_tumor-filtered.txt", check.names=FALSE)
xpr.gene <- cbind(xpr.gene.nt,xpr.gene.t)

rm(xpr.gene.nt,xpr.gene.t)

# Calculate switches
## Gene expression
xpr.gene.diff <- data.frame(matrix(nrow=nrow(xpr.gene),ncol=length(pats.t)))
colnames(xpr.gene.diff) <- pats.t
xpr.gene.diff[,pats.t.nt] <- abs(xpr.gene[,pats.t.nt] - xpr.gene[,pats.n])
xpr.gene.diff[,pats.t.t] <- abs(xpr.gene[,pats.t.t] - apply(xpr.gene[,pats.n],1,median))

### get p
xpr.gene.ecdf <- apply(xpr.gene[,pats.n],1,getEmpiricalDistribution,0.1)
xpr.gene.diff.p <- mapply(do.call, xpr.gene.ecdf, apply(xpr.gene.diff,1,list))
xpr.gene.diff.p <- lapply(xpr.gene.diff.p, function(x) 1 - x)
xpr.gene.diff.padj <- lapply(xpr.gene.diff.p,p.adjust)
xpr.gene.diff.padj <- do.call("rbind",xpr.gene.diff.padj)

### prepare df
xpr.gene.diff.padj <- as.data.frame(xpr.gene.diff.padj)
colnames(xpr.gene.diff.padj) <- pats.t
xpr.gene.diff.padj$Gene <- rownames(xpr.gene.diff.padj)

rm(xpr.gene.diff,xpr.gene.ecdf,xpr.gene.diff.p)

## PSI
### calculate differences
psi.diff <- data.frame(matrix(nrow=nrow(psi),ncol=length(pats.t)))
colnames(psi.diff) <- pats.t
psi.diff[,pats.t.nt] <- abs(psi[,pats.t.nt] - psi[,pats.n])
psi.diff[,pats.t.t] <- abs(psi[,pats.t.t] - apply(psi[,pats.n],1,median))

### get p
psi.ecdf <- apply(psi[,pats.n],1,getEmpiricalDistribution,0)
psi.diff.p <- mapply(do.call, psi.ecdf, apply(psi.diff,1,list))
psi.diff.p <- lapply(psi.diff.p, function(x) 1 - x)
psi.diff.padj <- lapply(psi.diff.p,p.adjust)
psi.diff.padj <- do.call("rbind",psi.diff.padj)

### prepare df
psi.diff.padj <- as.data.frame(psi.diff.padj)
colnames(psi.diff.padj) <- pats.t
psi.diff.padj$Gene <- genes
psi.diff.padj$Transcript <- transcripts

rm(psi.diff,psi.ecdf,psi.diff.p)

## Sign
psi.sign <- data.frame(matrix(nrow=nrow(psi),ncol=length(pats.t)))
colnames(psi.sign) <- pats.t
psi.sign[,pats.t.nt][psi[,pats.t.nt] > psi[,pats.n]] <- "Tumor"
psi.sign[,pats.t.t][psi[,pats.t.t] > apply(psi[,pats.n],1,median)] <- "Tumor"
psi.sign[,pats.t.nt][psi[,pats.t.nt] < psi[,pats.n]] <- "Normal"
psi.sign[,pats.t.t][psi[,pats.t.t] < apply(psi[,pats.n],1,median)] <- "Normal"

## Order
### Normal
psi.order.n <- data.frame(medianPsi=apply(psi[,pats.n],1,median))
psi.order.n$Gene <- unlist(strsplit(rownames(psi),","))[c(T,F)]
psi.order.n <- ddply(psi.order.n,.(Gene), summarise, orderNormal=order(-medianPsi))
rownames(psi.order.n) <- rownames(psi)

### Tumor
x <- by(psi[,pats.t], psi$Gene, function(y){
  apply(y[,colnames(y)!="Gene"],2,function(z){
    v <- rep(F,length(z))
    v[which.max(z)] <- T
    v
  })
})
psi.order.t <- do.call("rbind",x)
colnames(psi.order.t) <- pats.t

# Define switches
## Conditions
### transcripts are among the top expressed in the prefered condition
topTumor <- psi.sign=="Tumor" & psi.order.t
topNormal <- psi.sign=="Normal" & psi.order.n$orderNormal==1
### significant deltaPSI
bigChange <- psi.diff.padj[,pats.t] < 0.05
### non-significant change in gene expression
x <- merge(psi.diff.padj,xpr.gene.diff.padj,by="Gene",suffix=c(".psi.diff",".xpr.diff"))
noExpressionChange <- x[,paste0(pats.t,".xpr.diff")] > 0.05
### transcripts are expressed
expressedTumor <- xpr[,pats.t] > 0.1
medianN <- apply(xpr[,pats.n],1,median)
expressedNormal <- cbind(xpr[,pats.n], replicate(length(pats.t.t),medianN)) > 0.1

## Filter transcripts
elegibleTxs <- psi.sign
rownames(elegibleTxs) <- transcripts
elegibleTxs[!((topNormal & expressedNormal | topTumor & expressedTumor) & bigChange & noExpressionChange)] <- NA
elegibleTxs$Transcript <- transcripts
elegibleTxs$Gene <- genes

### filter out transcripts that do not pass any of the conditions
elegibleTxs <- elegibleTxs[rowSums(!is.na(elegibleTxs[,!colnames(elegibleTxs) %in% c("Gene","Transcript")]))>0,]

## Calculate switches
switches <- by(elegibleTxs,elegibleTxs$Gene, function(x){
  z <- apply(x[,! colnames(x) %in% c("Gene","Transcript")],2,function(y){
    z <- y[!is.na(y)]
    if("Normal" %in% z & "Tumor" %in% z)
      return(list("Normal"=names(z[z=="Normal"]), "Tumor"=names(z[z=="Tumor"])))
    else
      return(list("Normal"=NA, "Tumor"=NA))
  })
  
  z <- data.frame(do.call("rbind",z))
  z$Patient <- rownames(z)
  validCases <- rowSums(is.na(z[,c("Normal","Tumor")])) < 2
  z <- z[validCases,]
  if (nrow(z)){
    z$Gene <- unique(x$Gene)
    z$Normal <- as.character(z$Normal)
    z$Tumor <- as.character(z$Tumor)
    z 
  } else
    NA
  })

switches.formatted <- list()
for (g in names(switches)){
  if (class(switches[[g]])=="data.frame"){
    switches[[g]]$Gene=g
    switches.formatted[[g]] <- switches[[g]]
  }
}

switches.formatted <- do.call("rbind",switches.formatted)
switches.formatted <- switches.formatted[,c("Gene","Normal","Tumor","Patient")]

write.table(switches.formatted, file=outfile, sep="\t", row.names=F, col.names=F, quote=F)

stromal.correlation <- read.table(paste0("/projects_rg/TCGA/pipeline/run11/",tumor,"_gene_gsea_full.txt"), check.names=FALSE, header=TRUE)
