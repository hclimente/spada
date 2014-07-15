library(ggplot2)
library(reshape)
library(RColorBrewer)

load("~/SmartAS/Results/TCGA/luad_mE-1.0/RWorkspaces/2_GetCandidates.RData")
setwd("~/SmartAS")
nmb <- as.numeric(inputData["Replicates"])

#Plot simulated ideal cases
for (p in c(0.05, 0.49, 0.95)){
  thing <- data.frame(Normal=rnorm(n=nmb, m=0.3, sd=0.1), Tumor=c(rnorm(n=round(nmb*(1-p)), m=0.3, sd=0.1), rnorm(n=round(nmb*p), m=0.6, sd=0.1)))
  suppressMessages( thing.m <- melt(thing) )

  plt <- ggplot(data=thing.m) + geom_density(alpha=.7, aes(x = value, fill=factor(variable))) + 
         theme_minimal(base_size=20) + ylab("Density") + xlab("PSI") + 
         scale_fill_manual(values = c("firebrick2","steelblue3"), name="Location", breaks=c("Tumor", "Normal"), labels=c("Tumor", "Normal"))
  print(plt)
}

### PLOT SPECIFIC GENE INFORMATION ###
index <- 1
isof <- list()
gene <- as.character(candidateList$Gene[index])
isof[["N"]] <- as.character(candidateList$maxdPSI[index])
isof[["T"]] <- as.character(candidateList$mindPSI[index])

#Plot PSI of the transcripts of the first gene
Isos <- list()
myPalette <- brewer.pal(nrow(interReplicate[["N"]][interReplicate[["N"]]$Gene == gene,]), "Dark2")

for (i in c("N", "T")){
  Isos[[i]] <- interReplicate[[i]][interReplicate[[i]]$Gene == gene, c("Transcript", paste0("PSI_", 1:57) )]
  colnames(Isos[[i]]) <- c("Transcript", 1:57)
  suppressMessages(Isos[[i]] <- melt(Isos[[i]]))
  colnames(Isos[[i]]) <- c("Transcript", "Patient", "PSI")
  
  png(paste0("~/Desktop/", gene, "_", i, ".png"), width=1080, height=1080)
  plt <- ggplot() + geom_bar(data=Isos[[i]], stat="identity", aes(x=Patient, y=PSI, fill=factor(Transcript)), shape=1) +  
                    scale_fill_manual(values = myPalette, name="Transcript") + theme_minimal(base_size=20) +
                    theme(axis.text.x=element_text(angle=90, vjust=0.5, size=15))
  print(plt)
  graphics.off()
}

#Plot PSI of the normal and tumor isoforms
expression_N <- numeric()
expression_N <- c(expression_N, as.numeric(interReplicate[["N"]][interReplicate[["N"]]$Transcript == isof[["N"]], paste0("PSI_", 1:57) ]) )
expression_N <- c(expression_N, as.numeric(interReplicate[["T"]][interReplicate[["T"]]$Transcript == isof[["N"]], paste0("PSI_", 1:57) ]) )

expression_T <- numeric()
expression_T <- c(expression_T, as.numeric(interReplicate[["N"]][interReplicate[["N"]]$Transcript == isof[["T"]], paste0("PSI_", 1:57) ]) )
expression_T <- c(expression_T, as.numeric(interReplicate[["T"]][interReplicate[["T"]]$Transcript == isof[["T"]], paste0("PSI_", 1:57) ]) )

cand <- c(rep("N", nmb), rep("T", nmb) )

kks <- data.frame(Candidate=cand, ExpressionN=expression_N, ExpressionT=expression_T)

png(paste0("~/Dropbox/SmartAS/Presentaciones/", gene, "_PSI.png"), width=1300, height=1080)
ggplot(data=kks, aes(x=ExpressionN, y=ExpressionT, colour=Candidate )) + geom_point(size=8) +
  theme_minimal(base_size=25) + ylab("Tumor PSI") + xlab("Normal PSI") +
  scale_colour_manual(name="Transcript", values=c("N"="steelblue3", "T"="firebrick2"), labels=c("Normal", "Tumor") )
graphics.off()

#Make a plot of the FPR-based method
nPSI <- abs( as.numeric(interReplicate[["N"]][interReplicate[["N"]]$Transcript == isof[["N"]], paste0("PSI_", 1:57) ]) - as.numeric(interReplicate[["T"]][interReplicate[["T"]]$Transcript == isof[["N"]], paste0("PSI_", 1:57) ]))
fok <- data.frame(nPSI=nPSI)
foken <- melt(fok)

diffMatrix <- abs(outer(as.numeric(interReplicate[["N"]][interReplicate[["N"]]$Transcript == isof[["N"]], paste0("PSI_", 1:57) ]),as.numeric(interReplicate[["N"]][interReplicate[["N"]]$Transcript == isof[["N"]], paste0("PSI_", 1:57) ]),"-"))
subtraction <- data.frame(dPSI = diffMatrix[lower.tri(diffMatrix, diag = FALSE)])
fpr5 <- as.numeric( quantile(subtraction, 0.95, na.rm=T) )

signifFok=fok[(fok$nPSI > fpr5) & (fok$tPSI > fpr5),]
signifFoken <- melt(signifFok)

png(paste0("~/Dropbox/SmartAS/Presentaciones/", gene, "_method1.png"), width=1300, height=1080)
plt <- ggplot() + geom_density(data=foken, aes(x=value, fill=variable)) +
                  scale_x_continuous(limits = c(0, 0.4)) + theme_minimal(base_size=25) +
                  scale_y_continuous(limits = c(0,12)) + theme(legend.position="none")
print(plt)
graphics.off()

png(paste0("~/Dropbox/SmartAS/Presentaciones/", gene, "_method2.png"), width=1300, height=1080)
plt <- plt + geom_density(data=subtraction, alpha=.7, fill="darkolivegreen1", colout="", aes(x = dPSI) )
print(plt)
graphics.off()

png(paste0("~/Dropbox/SmartAS/Presentaciones/", gene, "_method3.png"), width=1300, height=1080)
plt1 <- plt + geom_vline(xintercept = fpr5, colour="red", linetype="dashed", size=5) +  guides(fill=FALSE)
print(plt1)
graphics.off()
