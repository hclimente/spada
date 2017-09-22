args <- commandArgs(trailingOnly = TRUE)
tpm.nt <- args[1]
tpm.t <- args[2]
psi.nt <- args[3]
psi.t <- args[4]
outfile <- args[5]

## read isoform expression
xpr.nt <- read.table(tpm.nt, check.names=FALSE)
xpr.t <- read.table(tpm.t, check.names=FALSE)
xpr <- cbind(xpr.nt,xpr.t)

rm(xpr.nt,xpr.t)

## read psi
psi.nt <- read.table(psi.nt, check.names=FALSE)
psi.t <- read.table(psi.t, check.names=FALSE)
psi <- cbind(psi.nt,psi.t)

rm(psi.nt,psi.t)

tumor <- grep("^.{4}T$",colnames(xpr), value=TRUE)
normal <- grep("^.{4}N$",colnames(xpr), value=TRUE)

med.xpr.t <- apply(xpr[,tumor],1,median)
med.xpr.n <- apply(xpr[,normal],1,median)
med.psi.t <- apply(psi[,tumor],1,median)
med.psi.n <- apply(psi[,normal],1,median)

xpr.agg <- data.frame(n.median.xpr=med.xpr.n,t.median.xpr=med.xpr.t,
                      n.median.psi=med.psi.n,t.median.psi=med.psi.t)

xpr.agg$Gene <- unlist(strsplit(row.names(xpr),","))[c(T,F)]
xpr.agg$Transcript <- unlist(strsplit(row.names(xpr),","))[c(F,T)]

write.table(xpr.agg[,c("Gene","Transcript","n.median.xpr","t.median.xpr","n.median.psi","t.median.psi")], 
            file=outfile, sep="\t", row.names=F, col.names=F, quote=F)