#!/soft/R/R-3.0.0/bin/Rscript

load("SmartAS.RData")
setwd(wd)

intraReplicate <- list()
interReplicate <- list()
isoformExpression <- list()

compartment <- inputData[["Compartments"]][1]
reference <- inputData[["Conditions"]][1]
alterated <- inputData[["Conditions"]][2]

printTPMHist <- function(x, xLab, pngName){
  png(paste0(wd,"/Results/DataExploration/", pngName, ".png"), width=960, height=960)
  histogram <- hist(log10(x + 0.0001), 10000)
  histogram$counts <- log10(histogram$counts)
  plot(histogram$mids, histogram$counts, type="h", main=pngName, xlab=xLab, ylab="log10(Frequency)")
  dev.off()
}

printLogFreqHist <- function(x, xLab, pngName){
  png(paste0(wd,"/Results/DataExploration/", pngName,".png"), width=960, height=960)
  histogram <- hist(x, 10000)
  histogram$counts <- log10(histogram$counts)
  plot(histogram$mids, histogram$counts, type="h", main=tag, xlab=xLab, ylab="log10(Frequency)")
  dev.off()
}

plotCorrelations <- function(x, y, lab, pngName){
  xLab=paste0(lab, " Replicate 1")
  yLab=paste0(lab, " Replicate 2")
  png(paste0(wd,"/Results/DataExploration/", pngName, ".png"), width=960, height=960)
  plot(x, y, xlab=xLab, ylab=yLab)
  dev.off()
  cor(x, y, use="complete.obs")
}

simplePlot <- function(x, y, title, xLab, yLab, pngName){
  png(pngName, width=960, height=960)
  plot(x, y, main=title, xlab=xLab, ylab=yLab)
  dev.off()
}

for (replicate in seq(1, inputData[["Replicates"]])){
  for (sample in inputData[["Conditions"]]){

    tag <- paste0(replicate, sample)
    cat("\t* Exploring replicate",replicate,"sample",sample,"\n")
    inputFile=paste0(wd, "/Data/Input/",  paste0(replicate, "_", sample), ".tsv")
    outputFile=paste0(wd, "/Results/", paste0(replicate, "_", sample), ".tsv")
      
    #Read Sailfish table
    isoformExpression[[tag]] <- read.table(inputFile, header=F, sep="\t", stringsAsFactors=F)
    colnames(isoformExpression[[tag]]) <- c("Gene", "Transcript","Genename","TPM")
      
    #Calculate the PSI for each transcript and the total expression of the gene
    vPSI <- as.numeric()
    vtTPM <- as.numeric()
      
    for (thisTranscript in 1:nrow(isoformExpression[[tag]])){
      mask <- isoformExpression[[tag]]$Gene==isoformExpression[[tag]]$Gene[thisTranscript]
      total <- sum(isoformExpression[[tag]]$TPM[mask])
      vtTPM <- c(vtTPM, total)
      if(total!=0){
        vPSI <- c(vPSI, isoformExpression[[tag]]$TPM[thisTranscript]/total)
      } else {
        vPSI <- c(vPSI, NA)
      }
    }
      
    isoformExpression[[tag]]$tTPM <- vtTPM
    isoformExpression[[tag]]$PSI <- vPSI
      
    write.table(isoformExpression[[tag]], file=outputFile, sep="\t", row.names=F)
  }
   
  nTag <- paste0(replicate, "N")
  tTag <- paste0(replicate, "T")
    
  intraReplicate[[replicate]] <- merge(isoformExpression[[nTag]], isoformExpression[[tTag]], by=c("Gene", "Transcript", "Genename"), suffixes=c("_N","_T"), all=T)
  intraReplicate[[replicate]]$deltaPSI <- intraReplicate[[replicate]]$PSI_N - intraReplicate[[replicate]]$PSI_T
  intraReplicate[[replicate]]$la_tTPM <- 0.5 * (log(intraReplicate[[replicate]]$tTPM_N) + log(intraReplicate[[replicate]]$tTPM_T))
    
  #Plot stuff
  printLogFreqHist(intraReplicate[[replicate]]$deltaPSI, "deltaPSI", paste0("deltaPSI_", replicate))
  printTPMHist(intraReplicate[[replicate]]$TPM_N, "log10(TPM_N+0.0001)", paste0("TPM_N_",replicate))
  printLogFreqHist(intraReplicate[[replicate]]$PSI_N, "PSI_N", paste0("PSI_N_",replicate))
  printTPMHist(intraReplicate[[replicate]]$TPM_T, "log10(TPM_T+0.0001)", paste0("TPM_T_",replicate))
  printLogFreqHist(intraReplicate[[replicate]]$PSI_N, "PSI_T", paste0("PSI_T_",replicate))
    
  write.table(intraReplicate[[replicate]], file=paste0(wd,"/Results/IntraReplicate_",replicate,".tsv"), sep="\t", row.names=F)

}

save(isoformExpression, intraReplicate, interReplicate, inputData, wd, file="SmartAS.RData")

for (r1 in set(1,inputData[["Replicates"]])){

  if (!exists("interReplicate[[\"Normal\"]]"))){
    interReplicate[["Normal"]] <- intraReplicate[[r1]]
    interReplicate[["Normal"]] <- subset(interReplicate[["Normal"]], select=-c(paste0("TPM_T_", r1), (paste0"tTPM_T_", r1), paste0("PSI_T_", r1)))
    colnames(interReplicate[["Normal"]]) <- c("Gene", "Transcript", "Genename", paste0("TPM_N_", r1), (paste0"tTPM_N_", r1), paste0("PSI_N_", r1)))
    
  }
  else {
    interReplicate[["Normal"]] <- merge(interReplicate[["Normal"]], intraReplicate[[r2]], by=c("Gene", "Transcript", "Genename"),
                                        suffixes=c("",paste0("_", r2)), all=T)
    interReplicate[["Normal"]] <- subset(interReplicate[["Normal"]], select=-c(paste0("TPM_T_", r1), (paste0"tTPM_T_", r1), paste0("PSI_T_", r1)))
    
  }

  #   thisComp <- paste0(r1, "_", r2)

  #   plotCorrelations(intraReplicate[[r1]]$TPM_N, intraReplicate[[r2]]$TPM_N, "TPM_N", paste0("TPM_N_cor_", thisComp))
  #   plotCorrelations(intraReplicate[[r1]]$TPM_T, intraReplicate[[r2]]$TPM_T, "TPM_T", paste0("TPM_T_cor_", thisComp))
  #   plotCorrelations(intraReplicate[[r1]]$deltaPSI, intraReplicate[[r2]]$deltaPSI, "deltaPSI", paste0("deltaPSI_cor_", thisComp))
  #   plotCorrelations(intraReplicate[[r1]]$PSI_N, intraReplicate[[r2]]$PSI_N, "PSI_N", paste0("PSI_N_cor_", thisComp))
  #   plotCorrelations(intraReplicate[[r1]]$PSI_T, intraReplicate[[r2]]$PSI_T, "PSI_T", paste0("PSI_T_cor_", thisComp))
      
  #   interReplicate[["Normal"]] <- merge(intraReplicate[[r1]], intraReplicate[[r2]], by=c("Gene", "Transcript", "Genename"),
  #                                       suffixes=c(paste0("_", r1),paste0("_", r2)), all=T)
  #   interReplicate[["Normal"]] <- subset(interReplicate[["Normal"]], select=-c(paste0("TPM_T_", r1), (paste0"tTPM_T_", r1), paste0("PSI_T_", r1), 
  #                                                                              paste0("TPM_T_", r2), paste0("tTPM_T_", r2), paste0("PSI_T_", r2)))
  #   interReplicate[["Normal"]]$deltaPSI <- interReplicate[["Normal"]]$PSI_N_1 - interReplicate[["Normal"]]$PSI_N_2
  #   interReplicate[["Normal"]]$la_tTPM <- 0.5 * (log(interReplicate[["Normal"]]$tTPM_N_1) + log(interReplicate[["Normal"]]$tTPM_N_2))
      
  #   interReplicate[["Tumor"]] <- merge(intraReplicate[[r1]], intraReplicate[[r2]], by=c("Gene", "Transcript", "Genename"), 
  #                                      suffixes=c(paste0("_", r1),paste0("_", r2)), all=T)
  #   interReplicate[["Tumor"]] <- subset(interReplicate[["Tumor"]], select=-c(TPM_N_1, tTPM_N_1, PSI_N_1, TPM_N_2, tTPM_N_2, PSI_N_2))
  #   interReplicate[["Tumor"]] <- subset(interReplicate[["Tumor"]], select=-c(paste0("TPM_N_", r1), (paste0"tTPM_N_", r1), paste0("PSI_N_", r1), 
  #                                                                            paste0("TPM_N_", r2), paste0("tTPM_N_", r2), paste0("PSI_N_", r2)))
  #   interReplicate[["Tumor"]]$deltaPSI <- interReplicate[["Tumor"]][,paste0("PSI_T_", r1)] - interReplicate[["Tumor"]][,paste0("PSI_T_", r2)]
  #   interReplicate[["Tumor"]]$la_tTPM <- 0.5 * (log(interReplicate[["Tumor"]][,paste0("tTPM_T_", r2)]) + log(interReplicate[["Tumor"]][,paste0("tTPM_T_", r2)]))
  # }
  simplePlot(intraReplicate[[r1]]$la_tTPM, intraReplicate[[r1]]$deltaPSI, r1, "0.5·(log(sum tTPM_N) + log(sum tTPM_T) )", 
             "deltaPSI", paste0(wd,"/Results/DataExploration/latTPM_PSI_intrarreplicate",r1,".png"))
  simplePlot(interReplicate[["Normal"]]$la_tTPM, interReplicate[["Normal"]]$deltaPSI, paste0(tag, "_N"), "0.5·(log(sum tTPM_1) + log(sum tTPM_2) )", 
            "deltaPSI", paste0(wd,"/Results/DataExploration/latTPM_PSI_interreplicate_N_",tag,".png"))
}

#Estimate the False Positive Rate
FPR <- as.numeric()
psiRange <- as.numeric()
for (psiValue in seq(0, 1, by=0.05)){
  expressionRange <- as.numeric()
  for (expValue in seq(-5, 8, by=0.25)){
    FPR <- c(FPR, sum((abs(interReplicate[["Normal"]]$deltaPSI)>=psiValue & interReplicate[["Normal"]]$la_tTPM>=expValue), na.rm=TRUE)/nrow(interReplicate[["Normal"]]))
    expressionRange <- c(expressionRange,expValue)
  }
  psiRange <- c(psiRange,psiValue)
}

FPRstudy <- matrix(data=FPR, nrow=length(expressionRange), ncol=length(psiRange))
rownames(FPRstudy) <- expressionRange
colnames(FPRstudy) <- psiRange

png(paste0(wd,"/Results/DataExploration/FPR.png"), width=960, height=960)
heatmap(FPRstudy, Rowv=NA, Colv=NA, scale="none", col = heat.colors(256), xlab="|deltaPSI|", ylab="Expression", main="FPR")
dev.off()

save(isoformExpression, intraReplicate, interReplicate, inputData, wd, file="SmartAS.RData")
