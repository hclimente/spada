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

createRow <- function(rawValue){
  newRow <- as.character()
  for (field in strsplit(as.character(rawValue), split="|", fixed=T)){
    newRow <- c(newRow,field)
  }
  if(length(newRow) < 10){
    for (i in seq(length(newRow), 9)){
      newRow <- c(newRow,NA)
    }
  }
  
  return(newRow)
}

for (kmer in inputData[["K-mer"]]){
  for (replicate in inputData[["Replicates"]]){
    tag <- paste0(compartment, replicate, "_", kmer)
    for (sample in inputData[["Conditions"]]){

      thisTag <- paste0(sample, tag)
      cat("\t* Exploring file",thisTag, "\n")
      inputFile=paste0(wd, "/Data/GENCODE/", thisTag, ".filtered.sf")
      outputFile=paste0(wd, "/Results/", thisTag, ".tsv")
      
      #Read Sailfish table
      isoformExpression[[thisTag]] <- read.table(inputFile, header=F, sep="\t", stringsAsFactors=F)
      colnames(isoformExpression[[thisTag]]) <- c("Gene", "Transcript","Genename","TPM")
      
      #Calculate the PSI for each transcript and the total expression of the gene
      vPSI <- as.numeric()
      vtTPM <- as.numeric()
      
      for (thisTranscript in 1:nrow(isoformExpression[[thisTag]])){
        mask <- isoformExpression[[thisTag]]$Gene==isoformExpression[[thisTag]]$Gene[thisTranscript]
        total <- sum(isoformExpression[[thisTag]]$TPM[mask])
        vtTPM <- c(vtTPM, total)
        if(total!=0){
          vPSI <- c(vPSI, isoformExpression[[thisTag]]$TPM[thisTranscript]/total)
        } else {
          vPSI <- c(vPSI, NA)
        }
      }
      
      isoformExpression[[thisTag]]$tTPM <- vtTPM
      isoformExpression[[thisTag]]$PSI <- vPSI
      
      write.table(isoformExpression[[thisTag]], file=outputFile, sep="\t", row.names=F)
    }
   
    refTag <- paste0(reference, tag)
    altTag <- paste0(alterated, tag)
    
    intraReplicate[[tag]] <- merge(isoformExpression[[refTag]], isoformExpression[[altTag]], by=c("Gene", "Transcript", "Genename"), suffixes=c("_ref","_alt"), all=T)
    #intraReplicate[[tag]] <- merge(isoformExpression[[refTag]], isoformExpression[[altTag]], by=c("Gene","Transcript","GENCODE"), suffixes=c("_ref","_alt"), all=T)
    intraReplicate[[tag]]$deltaPSI <- intraReplicate[[tag]]$PSI_ref - intraReplicate[[tag]]$PSI_alt
    intraReplicate[[tag]]$la_tTPM <- 0.5 * (log(intraReplicate[[tag]]$tTPM_ref) + log(intraReplicate[[tag]]$tTPM_alt))
    
    #Plot stuff
    printLogFreqHist(intraReplicate[[tag]]$deltaPSI, "deltaPSI", paste0("deltaPSI_", tag))
    printTPMHist(intraReplicate[[tag]]$TPM_ref, "log10(TPM10+0.0001)", paste0("TPM10_",tag))
    printLogFreqHist(intraReplicate[[tag]]$PSI_ref, "PSI10", paste0("PSI10_",tag))
    printTPMHist(intraReplicate[[tag]]$TPM_alt, "log10(TPM7+0.0001)", paste0("TPM7_",tag))
    printLogFreqHist(intraReplicate[[tag]]$PSI_ref, "PSI7", paste0("PSI7_",tag))
    
    write.table(intraReplicate[[tag]], file=paste0(wd,"/Results/IntraReplicate",tag,".tsv"), sep="\t", row.names=F)

  }

  #Plot more stuff
  tag <- paste0(compartment, kmer)
  tag1 <- paste0(compartment, inputData[["Replicates"]][1], "_", kmer)
  tag2 <- paste0(compartment, inputData[["Replicates"]][2], "_", kmer)
  
  plotCorrelations(intraReplicate[[tag1]]$TPM_ref, intraReplicate[[tag2]]$TPM_ref, "TPM_10", paste0("TPM10_cor_", tag))
  plotCorrelations(intraReplicate[[tag1]]$TPM_alt, intraReplicate[[tag2]]$TPM_alt, "TPM_7", paste0("TPM7_cor_", tag))
  plotCorrelations(intraReplicate[[tag1]]$deltaPSI, intraReplicate[[tag2]]$deltaPSI, "deltaPSI", paste0("deltaPSI_cor_", tag))
  plotCorrelations(intraReplicate[[tag1]]$PSI_ref, intraReplicate[[tag2]]$PSI_ref, "PSI_10", paste0("PSI10_cor_", tag))
  plotCorrelations(intraReplicate[[tag1]]$PSI_alt, intraReplicate[[tag2]]$PSI_alt, "PSI_10", paste0("PSI7_cor_", tag))
  
  interReplicate[["Ref"]] <- merge(intraReplicate[[tag1]], intraReplicate[[tag2]], by=c("Gene", "Transcript", "Genename"), suffixes=c("_1","_2"), all=T)
  #interReplicate[["Ref"]] <- merge(intraReplicate[[tag1]], intraReplicate[[tag2]], by=c("Gene","Transcript","GENCODE"), suffixes=c("_1","_2"), all=T)
  interReplicate[["Ref"]] <- subset(interReplicate[["Ref"]], select=-c(TPM_alt_1, tTPM_alt_1, PSI_alt_1, TPM_alt_2, tTPM_alt_2, PSI_alt_2))
  interReplicate[["Ref"]]$deltaPSI <- interReplicate[["Ref"]]$PSI_ref_1 - interReplicate[["Ref"]]$PSI_ref_2
  interReplicate[["Ref"]]$la_tTPM <- 0.5 * (log(interReplicate[["Ref"]]$tTPM_ref_1) + log(interReplicate[["Ref"]]$tTPM_ref_2))
  
  interReplicate[["Alt"]] <- merge(intraReplicate[[tag1]], intraReplicate[[tag2]], by=c("Gene", "Transcript", "Genename"), suffixes=c("_1","_2"), all=T)
  #interReplicate[["Alt"]] <- merge(intraReplicate[[tag1]], intraReplicate[[tag2]], by=c("Gene","Transcript","GENCODE"), suffixes=c("_1","_2"), all=T)
  interReplicate[["Alt"]] <- subset(interReplicate[["Alt"]], select=-c(TPM_ref_1, tTPM_ref_1, PSI_ref_1, TPM_ref_2, tTPM_ref_2, PSI_ref_2))
  interReplicate[["Alt"]]$deltaPSI <- interReplicate[["Alt"]]$PSI_alt_1 - interReplicate[["Alt"]]$PSI_alt_2
  interReplicate[["Alt"]]$la_tTPM <- 0.5 * (log(interReplicate[["Alt"]]$tTPM_alt_1) + log(interReplicate[["Alt"]]$tTPM_alt_2))
  
  simplePlot(intraReplicate[[tag1]]$la_tTPM, intraReplicate[[tag1]]$deltaPSI, tag1, "0.5·(log(sum tTPM10) + log(sum tTPM7) )", 
             "deltaPSI", paste0(wd,"/Results/DataExploration/latTPM_PSI_intrarreplicate",tag1,".png"))
  simplePlot(intraReplicate[[tag2]]$la_tTPM, intraReplicate[[tag2]]$deltaPSI, tag2, "0.5·(log(sum tTPM10) + log(sum tTPM7) )", 
             "deltaPSI", paste0(wd,"/Results/DataExploration/latTPM_PSI_intrarreplicate",tag2,".png"))
  simplePlot(interReplicate[["Ref"]]$la_tTPM, interReplicate[["Ref"]]$deltaPSI, paste0(tag, "_10"), "0.5·(log(sum tTPM_1) + log(sum tTPM_2) )", 
             "deltaPSI", paste0(wd,"/Results/DataExploration/latTPM_PSI_interreplicate10_",tag,".png"))
  simplePlot(interReplicate[["Alt"]]$la_tTPM, interReplicate[["Alt"]]$deltaPSI, paste0(tag, "_7"), "0.5·(log(sum tTPM_1) + log(sum tTPM_2) )", 
             "deltaPSI", paste0(wd,"/Results/DataExploration/latTPM_PSI_interreplicate7_",tag,".png"))
}

#Estimate the False Positive Rate
FPR <- as.numeric()
psiRange <- as.numeric()
for (psiValue in seq(0, 1, by=0.05)){
  expressionRange <- as.numeric()
  for (expValue in seq(-5, 8, by=0.25)){
    FPR <- c(FPR, sum((abs(interReplicate[["Ref"]]$deltaPSI)>=psiValue & interReplicate[["Ref"]]$la_tTPM>=expValue), na.rm=TRUE)/nrow(interReplicate[["Ref"]]))
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