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
tpm.nt <- args[1]
tpm.t <- args[2]
tpm.g.nt <- args[3]
tpm.g.t <- args[4]
psi.nt <- args[5]
psi.t <- args[6]
outfile <- args[7]

# Prepare the data
## read psi
psi.nt <- read.table(psi.nt, check.names=FALSE)
psi.t <- read.table(psi.t, check.names=FALSE)
psi <- cbind(psi.nt,psi.t)

rm(psi.nt,psi.t)

## get genes and transcripts
genes <- unlist(strsplit(row.names(psi),","))[c(T,F)]
transcripts <- unlist(strsplit(row.names(psi),","))[c(F,T)]
psi$Gene <- genes

## get sample names
normal <- grep("^.{4}N$",colnames(psi), value=TRUE)
tumor <- grep("^.{4}T$",colnames(psi), value=TRUE)
tumor.paired <- gsub("N$","T",normal)
tumor.unpaired <- setdiff(tumor,tumor.paired)

## read isoform expression
xpr.nt <- read.table(tpm.nt, check.names=FALSE)
xpr.t <- read.table(tpm.t, check.names=FALSE)
xpr <- cbind(xpr.nt,xpr.t)

rm(xpr.nt,xpr.t)

## read gene expression
xpr.gene.nt <- read.table(tpm.g.nt, check.names=FALSE)
xpr.gene.t <- read.table(tpm.g.t, check.names=FALSE)
xpr.gene <- cbind(xpr.gene.nt,xpr.gene.t)

### transform to log
### Some cases with very high variability between normal samples AND orders of
### magnitude above tumor samples where passing this test using raw difference.
### Logarithm differences should alleviate this problem
### Example: ADAMTS8|11095, Normal:uc001qgg.3, Tumor:uc001qgf.2, Sample: A46PT
logxpr.gene <- log2(xpr.gene)

rm(xpr.gene,xpr.gene.nt,xpr.gene.t)

# Calculate switches
## Gene expression
logxpr.gene.diff <- data.frame(matrix(nrow=nrow(logxpr.gene),ncol=length(tumor)))
colnames(logxpr.gene.diff) <- tumor
logxpr.gene.diff[,tumor.paired] <- abs(logxpr.gene[,tumor.paired] - logxpr.gene[,normal])
logxpr.gene.diff[,tumor.unpaired] <- abs(logxpr.gene[,tumor.unpaired] - apply(logxpr.gene[,normal],1,median))

### get p
logxpr.gene.ecdf <- apply(logxpr.gene[,normal],1,getEmpiricalDistribution,0.1)
logxpr.gene.diff.p <- mapply(do.call, logxpr.gene.ecdf, apply(logxpr.gene.diff,1,list))
logxpr.gene.diff.p <- lapply(logxpr.gene.diff.p, function(x) 1 - x)
logxpr.gene.diff.padj <- lapply(logxpr.gene.diff.p,p.adjust)
logxpr.gene.diff.padj <- do.call("rbind",logxpr.gene.diff.padj)

### prepare df
logxpr.gene.diff.padj <- as.data.frame(logxpr.gene.diff.padj)
colnames(logxpr.gene.diff.padj) <- tumor
logxpr.gene.diff.padj$Gene <- rownames(logxpr.gene.diff.padj)

rm(logxpr.gene.diff,logxpr.gene.ecdf,logxpr.gene.diff.p)

## PSI
### calculate differences
psi.diff <- data.frame(matrix(nrow=nrow(psi),ncol=length(tumor)))
colnames(psi.diff) <- tumor
psi.diff[,tumor.paired] <- abs(psi[,tumor.paired] - psi[,normal])
psi.diff[,tumor.unpaired] <- abs(psi[,tumor.unpaired] - apply(psi[,normal],1,median))

### get p
psi.ecdf <- apply(psi[,normal],1,getEmpiricalDistribution,0)
psi.diff.p <- mapply(do.call, psi.ecdf, apply(psi.diff,1,list))
psi.diff.p <- lapply(psi.diff.p, function(x) 1 - x)
psi.diff.padj <- lapply(psi.diff.p,p.adjust)
psi.diff.padj <- do.call("rbind",psi.diff.padj)

### prepare df
psi.diff.padj <- as.data.frame(psi.diff.padj)
colnames(psi.diff.padj) <- tumor
psi.diff.padj$Gene <- genes
psi.diff.padj$Transcript <- transcripts

rm(psi.diff,psi.ecdf,psi.diff.p)

## Sign
psi.sign <- data.frame(matrix(nrow=nrow(psi),ncol=length(tumor)))
colnames(psi.sign) <- tumor
psi.sign[,tumor.paired][psi[,tumor.paired] > psi[,normal]] <- "Tumor"
psi.sign[,tumor.unpaired][psi[,tumor.unpaired] > apply(psi[,normal],1,median)] <- "Tumor"
psi.sign[,tumor.paired][psi[,tumor.paired] < psi[,normal]] <- "Normal"
psi.sign[,tumor.unpaired][psi[,tumor.unpaired] < apply(psi[,normal],1,median)] <- "Normal"

## Order
### Normal
psi.order.n <- data.frame(medianPsi=apply(psi[,normal],1,median))
psi.order.n$Gene <- unlist(strsplit(rownames(psi),","))[c(T,F)]
psi.order.n <- ddply(psi.order.n,.(Gene), summarise, orderNormal=order(-medianPsi))
rownames(psi.order.n) <- rownames(psi)

### Tumor
x <- by(psi[,tumor], psi$Gene, function(y){
  apply(y[,colnames(y)!="Gene"],2,function(z){
    v <- rep(F,length(z))
    v[which.max(z)] <- T
    v
  })
})
psi.order.t <- do.call("rbind",x)
colnames(psi.order.t) <- tumor

# Define switches
## Conditions
### transcripts are among the top expressed in the prefered condition
topTumor <- psi.sign=="Tumor" & psi.order.t
topNormal <- psi.sign=="Normal" & psi.order.n$orderNormal==1
### significant deltaPSI
bigChange <- psi.diff.padj[,tumor] < 0.05
### non-significant change in gene expression
x <- merge(psi.diff.padj,logxpr.gene.diff.padj,by="Gene",suffix=c(".psi.diff",".xpr.diff"))
noExpressionChange <- x[,paste0(tumor,".xpr.diff")] > 0.05
### transcripts are expressed
expressedTumor <- xpr[,tumor] > 0.1
medianN <- apply(xpr[,normal],1,median)
expressedNormal <- cbind(xpr[,normal], replicate(length(tumor.unpaired),medianN)) > 0.1

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
  z$Sample <- rownames(z)
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

switches.df <- list()
for (g in names(switches)){
  if (class(switches[[g]])=="data.frame"){
    switches[[g]]$Gene=g
    switches.df[[g]] <- switches[[g]]
  }
}

switches.df <- do.call("rbind",switches.df)
switches.df <- switches.df[,c("Gene","Normal","Tumor","Sample")]

# we expect the switches to be distributed among paired and
# unpaired samples, those that do not comply with it are
# likely to be caused by a lack of representativeness of the median
# when calculating the deltaPSI
switches.df.byPat <- ddply(switches.df,.(Gene,Normal,Tumor),summarise,
                           Paired=sum(Sample %in% tumor.paired), 
                           Unpaired=sum(Sample %in% tumor.unpaired) )

p <- apply(switches.df.byPat[,c("Paired","Unpaired")],1,
      function(x,p){
        b <- binom.test(x[1],x[1]+x[2],p,"less")
        b$p.value
      },length(tumor.paired)/length(tumor.unpaired))

# remove those genes significantly unbalance towards unpaired patients
# also those that have 0 paired patients, as those are the switches 
# where the method applies
switches.df.byPat.filt <- subset(switches.df.byPat, p >= 0.05 & Paired>0, select=c("Gene","Normal","Tumor"))
switches.df.filt <- merge(switches.df.byPat.filt,switches.df)

switches.df.filt.agg <- ddply(switches.df.filt,.(Gene,Normal,Tumor),summarise,
                              Samples=paste(Sample,collapse=","))

write.table(switches.df.filt, file=outfile, sep="\t", row.names=F, col.names=F, quote=F)

#stromal.correlation <- read.table(paste0("/projects_rg/TCGA/pipeline/run11/",tumor,"_gene_gsea_full.txt"), check.names=FALSE, header=TRUE)