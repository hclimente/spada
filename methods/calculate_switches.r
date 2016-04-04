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

library(plyr)
library(edgeR)

args <- commandArgs(trailingOnly = TRUE)
tpm.nt.file <- args[1]
tpm.t.file <- args[2]
counts.g.nt.file <- args[3]
counts.g.t.file <- args[4]
psi.nt.file <- args[5]
psi.t.file <- args[6]
outfile <- args[7]

# Prepare the data
## read psi
psi.nt <- read.table(psi.nt.file, check.names=FALSE)
psi.t <- read.table(psi.t.file, check.names=FALSE)
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
xpr.nt <- read.table(tpm.nt.file, check.names=FALSE)
xpr.t <- read.table(tpm.t.file, check.names=FALSE)
xpr <- cbind(xpr.nt,xpr.t)

rm(xpr.nt,xpr.t)

## read gene expression
xpr.gene.nt <- read.table(counts.g.nt.file, check.names=FALSE)
xpr.gene.t <- read.table(counts.g.t.file, check.names=FALSE)
xpr.gene <- cbind(xpr.gene.nt,xpr.gene.t)

### tmm normalization
#origin <- factor(substring(colnames(xpr.gene),5,5))
y <- DGEList(counts=xpr.gene)
y <- calcNormFactors(y)
xpr.gene.norm <- cpm(y, normalized.lib.sizes=TRUE)
#xpr.gene.norm <- as.data.frame(y$counts/y$samples$norm.factors)

### transform to log
### Some cases with very high variability between normal samples AND orders of
### magnitude above tumor samples where passing this test using raw difference.
### Logarithm differences should alleviate this problem
### Example: ADAMTS8|11095, Normal:uc001qgg.3, Tumor:uc001qgf.2, Sample: A46PT
logxpr.gene <- log2(xpr.gene.norm)
logxpr.gene[logxpr.gene==-Inf] <- NA

rm(xpr.gene,xpr.gene.nt,xpr.gene.t,xpr.gene.norm)

# Calculate switches
## Gene expression
logxpr.gene.diff <- data.frame(matrix(nrow=nrow(logxpr.gene),ncol=length(tumor)))
colnames(logxpr.gene.diff) <- tumor
logxpr.gene.diff[,tumor.paired] <- abs(logxpr.gene[,tumor.paired] - logxpr.gene[,normal])
logxpr.gene.diff[,tumor.unpaired] <- abs(logxpr.gene[,tumor.unpaired] - apply(logxpr.gene[,normal],1,median))

### get p
# use raw pvalues: we want to minimize false negatives, and multiple
# test correction tends to minimize false positives
logxpr.gene.ecdf <- apply(logxpr.gene[,normal],1,getEmpiricalDistribution,0.1)
logxpr.gene.diff.p <- mapply(do.call, logxpr.gene.ecdf, apply(logxpr.gene.diff,1,list))
logxpr.gene.diff.p <- 1 - logxpr.gene.diff.p

### prepare df
logxpr.gene.diff.p <- as.data.frame(t(logxpr.gene.diff.p))
colnames(logxpr.gene.diff.p) <- tumor
logxpr.gene.diff.p$Gene <- rownames(logxpr.gene.diff.p)

rm(logxpr.gene.diff,logxpr.gene.ecdf)

## PSI
### calculate differences
psi.diff <- data.frame(matrix(nrow=nrow(psi),ncol=length(tumor)))
colnames(psi.diff) <- tumor
psi.diff[,tumor.paired] <- abs(psi[,tumor.paired] - psi[,normal])
psi.diff[,tumor.unpaired] <- abs(psi[,tumor.unpaired] - apply(psi[,normal],1,median))
psi.diff.l <- as.list(as.data.frame(t(psi.diff)))
psi.diff.v <- as.vector(t(psi.diff))

### get p
psi.ecdf <- apply(psi[,normal],1,getEmpiricalDistribution,0)

psi.diff.p <- mapply(function(x,y){x(y)}, psi.ecdf, psi.diff.l)
psi.diff.p <- 1 - psi.diff.p

#### p-values when the difference is 0 are 1, 
#### instead of the value < 1 attributed by ecdf
psi.diff.p[!is.na(psi.diff.v) & psi.diff.v==0] <- 1

### prepare df
psi.diff.p <- as.data.frame(t(psi.diff.p))
colnames(psi.diff.p) <- tumor
psi.diff.p$Gene <- genes
psi.diff.p$Transcript <- transcripts

rm(psi.ecdf,psi.diff.l,psi.diff.v)

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
psi.order.n <- ddply(psi.order.n,.(Gene), summarise, orderNormal=rank(-medianPsi))
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
bigChange <- psi.diff.p[,tumor] < 0.01 & psi.diff > 0.05
### non-significant change in gene expression
x <- merge(psi.diff.p,logxpr.gene.diff.p,by="Gene",suffix=c(".psi.diff",".xpr.diff"),all.x=T)
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
      },length(tumor.paired)/(length(tumor.unpaired)+length(tumor.paired)))

# remove those genes significantly unbalance towards unpaired patients
# also those that have 0 paired patients, as those are the switches 
# where the method applies
switches.df.byPat.filt <- subset(switches.df.byPat, p >= 0.05 & Paired>0, select=c("Gene","Normal","Tumor"))
switches.df.filt <- merge(switches.df.byPat.filt,switches.df)

switches.df.filt.agg <- ddply(switches.df.filt,.(Gene,Normal,Tumor),summarise,
                              Samples=paste(Sample,collapse=","))

write.table(switches.df.filt.agg, file=outfile, sep="\t", row.names=F, col.names=F, quote=F)