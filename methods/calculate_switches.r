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
    
    if (all(subtraction==0)){
      return(function(x){
        if(x==0){
          return(0)
        } else if (x>0) {
          return(1)
        } else {
          return(NA)
        }
      })
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

## change psi colnames
colnames(psi) <- paste0(colnames(psi),".psi")

## read gene xpr
xpr.gene.nt <- read.table("/projects_rg/TCGA/pipeline/run11/luad_gene_tpm_paired-filtered.txt", check.names=FALSE)
xpr.gene.t <- read.table("/projects_rg/TCGA/pipeline/run11/luad_gene_tpm_tumor-filtered.txt", check.names=FALSE)
xpr.gene <- cbind(xpr.gene.nt,xpr.gene.t)
colnames(xpr.gene) <- paste0(colnames(xpr.gene),".xpr")

rm(xpr.gene.nt,xpr.gene.t)

# Calculate switches
## Gene expression
xpr.gene.diff <- data.frame(matrix(nrow=nrow(xpr.gene),ncol=length(pats.t)))
colnames(xpr.gene.diff) <- pats.t
xpr.gene.diff[,pats.t.nt] <- abs(xpr.gene[,paste0(pats.t.nt,".xpr")] - xpr.gene[,paste0(pats.n,".xpr")])
xpr.gene.diff[,pats.t.t] <- abs(xpr.gene[,paste0(pats.t.t,".xpr")] - apply(xpr.gene[,paste0(pats.n,".xpr")],1,median))

### get p
xpr.gene.ecdf <- apply(xpr.gene[,paste0(pats.n,".xpr")],1,getEmpiricalDistribution,0.1)
xpr.gene.diff.p <- mapply(do.call, xpr.gene.ecdf, apply(xpr.gene.diff,1,list))
xpr.gene.diff.p <- lapply(xpr.gene.diff.p, function(x) 1 - x)
xpr.gene.diff.padj <- lapply(xpr.gene.diff.p,p.adjust)
xpr.gene.diff.padj <- do.call("rbind",xpr.gene.diff.padj)

### prepare df
xpr.gene.diff.padj <- as.data.frame(xpr.gene.diff.padj)
colnames(xpr.gene.diff.padj) <- paste0(pats.t,".xpr.diff")
xpr.gene.diff.padj$Gene <- rownames(xpr.gene.diff.padj)

rm(xpr.gene.diff,xpr.gene.ecdf,xpr.gene.diff.p)

## PSI
### calculate differences
psi.diff <- data.frame(matrix(nrow=nrow(psi),ncol=length(pats.t)))
colnames(psi.diff) <- pats.t
psi.diff[,pats.t.nt] <- abs(psi[,paste0(pats.t.nt,".psi")] - psi[,paste0(pats.n,".psi")])
psi.diff[,pats.t.t] <- abs(psi[,paste0(pats.t.t,".psi")] - apply(psi[,paste0(pats.n,".psi")],1,median))

### get p
psi.ecdf <- apply(psi[,paste0(pats.n,".psi")],1,getEmpiricalDistribution,0)
psi.diff.p <- mapply(do.call, psi.ecdf, apply(psi.diff,1,list))
psi.diff.p <- lapply(psi.diff.p, function(x) 1 - x)
psi.diff.padj <- lapply(psi.diff.p,p.adjust)
psi.diff.padj <- do.call("rbind",psi.diff.padj)

### prepare df
psi.diff.padj <- as.data.frame(psi.diff.padj)
colnames(psi.diff.padj) <- paste0(pats.t,".psi.diff")
psi.diff.padj$Gene <- genes
psi.diff.padj$Transcript <- transcripts

rm(psi.diff,psi.ecdf,psi.diff.p)

## Sign
psi.sign <- data.frame(matrix(nrow=nrow(psi),ncol=length(pats.t)))
colnames(psi.sign) <- paste0(pats.t,".sign")
psi.sign[,paste0(pats.t.nt,".sign")][psi[,paste0(pats.t.nt,".psi")] > psi[,paste0(pats.n,".psi")]] <- "Tumor"
psi.sign[,paste0(pats.t.t,".sign")][psi[,paste0(pats.t.t,".psi")] > apply(psi[,paste0(pats.n,".psi")],1,median)] <- "Tumor"
psi.sign[,paste0(pats.t.nt,".sign")][psi[,paste0(pats.t.nt,".psi")] < psi[,paste0(pats.n,".psi")]] <- "Normal"
psi.sign[,paste0(pats.t.t,".sign")][psi[,paste0(pats.t.t,".psi")] < apply(psi[,paste0(pats.n,".psi")],1,median)] <- "Normal"

## Order
### Normal
psi.order.n <- data.frame(medianPsi=apply(psi[,paste0(pats.n,".psi")],1,median))
psi.order.n$Gene <- unlist(strsplit(rownames(psi),","))[c(T,F)]
psi.order.n <- ddply(psi.order.n,.(Gene), summarise, orderNormal=order(-medianPsi))
rownames(psi.order.n) <- rownames(psi)

### Tumor
x <- by(psi[,paste0(pats.t,".psi")], psi$Gene, function(y){
  apply(y[,colnames(y)!="Gene"],2,function(z){
    v <- rep(F,length(z))
    v[which.max(z)] <- T
    v
    })})
psi.order.t <- do.call("rbind",x)
colnames(psi.order.t) <- paste0(pats.t,".orderTumor")

# Define switches
## Conditions
### transcripts are among the top expressed in the prefered condition
topTumor <- psi.sign=="Tumor" & psi.order.t
topNormal <- psi.sign=="Normal" & psi.order.n$orderNormal==1
### significant deltaPSI
bigChange <- psi.diff.padj[,paste0(pats.t,".psi.diff")] < 0.05
### non-significant change in gene expression
#### CHECK ORDER DOESNT CHANGE
x <- merge(psi.diff.padj,xpr.gene.diff.padj,by="Gene")
noExpressionChange <- x[,paste0(pats.t,".xpr.diff")] > 0.05
### transcripts are expressed
xpr.gene$Gene <- rownames(xpr.gene)
x <- merge(psi.diff.padj,xpr.gene,by="Gene")
expressedTumor <- x[,paste0(pats.t,".xpr")] > 0.1
y <- apply(x[,paste0(pats.n,".xpr")],1,median)
expressedNormal <- cbind(x[,paste0(pats.n,".xpr")], replicate(length(pats.t.t),y)) > 0.1

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
  z$Patient <- gsub(".sign","",rownames(z))
  validCases <- rowSums(is.na(z[,c("Normal","Tumor")])) < 2
  z <- z[validCases,]
  if (nrow(z)){
    z$Gene <- unique(x$Gene)
    z 
  } else
    NA
  })

switches.formatted <- list()
for (g in names(switches)){
  if (is.na(switches[[g]]))
    next
  else{
    switches[[g]]$Gene=g
    switches.formatted[[g]] <- switches[[g]]
  }
}

switches.formatted <- do.call("rbind",switches.formatted)
switches.formatted <- switches.formatted[,c("Gene","Normal","Tumor","Patient")]

write.table(switches.formatted, file=outfile, sep="\t", row.names=F, col.names=F, quote=F)

# stromal.correlation <- read.table(paste0("/projects_rg/TCGA/pipeline/run11/",tumor,"_gene_gsea_full.txt"), check.names=FALSE, header=TRUE)
