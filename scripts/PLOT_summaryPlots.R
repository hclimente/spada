library(plyr)
library(ggplot2)
library(reshape2)
library(directlabels)

cancerTypes <- c("brca","coad","hnsc","kich","kirc","kirp","lihc","luad","lusc","prad","thca")
workingDir <- "/genomics/users/hector/TCGA_analysis"
source("~/SmartAS/Pipeline/scripts/smartas_theme.R")
setwd(workingDir)

colorPalette <- c("#CC79A7","#636363","#F0E442","#006D2C","#31A354","#74C476","#FC8D59","#08519C","#3182BD","#D55E00","#5E3C99","#000000","#969696","#EDF8FB","#B3CDE3","#8C96C6","#88419D","#D94701","#2171B5","#FD8D3C","#6BAED6","#FD8D3C","#33A02C","#E31A1C","#FF7F00","#6A3D9A","#B15928","#377EB8","#E41A1C")
names(colorPalette) <- c("brca","coad","hnsc","kich","kirc","kirp","lihc","luad","lusc","prad","thca","coad-hyper","coad-hypo","brca-basal","brca-her2","brca-luminala","brca-luminalb","expression-up","expression-down","psi-up","psi-down","odds","a3","a5","mx","ri","se","normal","tumor")

nPatients <- list()
nPatients[["brca"]] <- 1036
nPatients[["coad"]] <- 262
nPatients[["hnsc"]] <- 422
nPatients[["kich"]] <- 62
nPatients[["kirc"]] <- 505
nPatients[["kirp"]] <- 195
nPatients[["lihc"]] <- 197
nPatients[["luad"]] <- 488
nPatients[["lusc"]] <- 483
nPatients[["prad"]] <- 295
nPatients[["thca"]] <- 497

nPatientsDf <- as.data.frame(do.call("rbind",nPatients))
nPatientsDf$Cancer <- rownames(nPatientsDf)
colnames(nPatientsDf) <- c("TotalPatients","Cancer")

################ CANDIDATES STUDY ################
candidateList <- list()
for (cancer in cancerTypes){
	candidateList[[cancer]] <- read.delim(paste0(cancer,".candidateList.tsv"))
	candidateList[[cancer]] <- candidateList[[cancer]][candidateList[[cancer]]$IsModel==1 & candidateList[[cancer]]$IsRelevant==1 & candidateList[[cancer]]$NotNoise==1,]
	candidateList[[cancer]]$NumPatients <- unlist(lapply(strsplit(as.character(candidateList[[cancer]]$Patients_affected),",",fixed=T),length))
	candidateList[[cancer]]$Percentage <- candidateList[[cancer]]$NumPatients/nPatients[[cancer]]
	candidateList[[cancer]]$Origin <- cancer
	#plotcmd <- "plot(log(candidateList_v3$Patient_percentage),-log(candidateList_v3$Sensitivity))"
	#save.plot(plotcmd, file=paste0(cancer,"_patient-sensitivity.png"), dir=getwd(),w=1000, h=1000, format="png")
	#plotcmd <- "plot(sqrt(candidateList_v3$Patient_percentage),candidateList_v3$Precision)"
	#save.plot(plotcmd, file=paste0(cancer,"_patient-precision.png"), dir=getwd(),w=1000, h=1000, format="png")
}
candidateList2 <- do.call("rbind", candidateList)

# calculate entropy
tempTable <- merge(candidateList2,nPatientsDf,by.x="Origin",by.y="Cancer")
unbalaceMeasure <- ddply(tempTable,.(GeneId,Symbol,Normal_transcript,Tumor_transcript), summarise, unbalanceP=min(apply(cbind(NumPatients,sum(NumPatients)-NumPatients,TotalPatients-NumPatients,sum(TotalPatients)-TotalPatients-(sum(NumPatients)-NumPatients)),1,function(x) { if (x[2]!=0){fish <- fisher.test(matrix(x[1:4],ncol=2),alternative="greater"); fish$p.value} else { 0 } } )))

candidateList2 <- merge(unbalaceMeasure,candidateList2,by=c("GeneId","Symbol","Normal_transcript","Tumor_transcript"))

#most frequent genes, ordered by sum of percentajes in cancer types
for (annotation in c("All","Driver","d1")){
  ngenes <- 50
  
  if (annotation=="All"){
    selectVector <- matrix(T,nrow(candidateList2),1)
  } else {
    selectVector <- candidateList2$Annotation==annotation
  }
  
  accCandidateList <- ddply(candidateList2[selectVector,],.(GeneId,Symbol), summarise, cumsum=sum(Percentage),patients=sum(NumPatients))
  accCandidateList <- accCandidateList[order(-accCandidateList$cumsum),]
  
  candsMatrix <- list()
  for (gene in accCandidateList$GeneId[1:ngenes]){
    symbol = accCandidateList$Symbol[accCandidateList$GeneId==gene]
    candsMatrix[[symbol]] <- list()
    candsMatrix[[symbol]][["symbol"]] <- as.character(symbol)
    for (cancer in cancerTypes){
      if (gene %in% candidateList[[cancer]]$GeneId){
        candsMatrix[[symbol]][[cancer]] <- as.numeric(candidateList[[cancer]]$Percentage[candidateList[[cancer]]$GeneId==gene])
      } else{
        candsMatrix[[symbol]][[cancer]] <- 0
      }
    }
    candsMatrix[[symbol]] <- do.call("cbind", candsMatrix[[symbol]])
  }
  candsMatrix2 <- do.call("rbind", candsMatrix)
  rownames(candsMatrix2) <- candsMatrix2[,c("symbol")]
  
  candsMatrix3 <- as.data.frame(matrix(as.numeric(as.character(candsMatrix2[,!(colnames(candsMatrix2) %in% c("symbol"))])),nrow=ngenes,ncol=11))
  rownames(candsMatrix3) <- rownames(candsMatrix2)
  colnames(candsMatrix3) <- cancerTypes
  candsMatrix3 <- candsMatrix3[order(-rowSums(candsMatrix3)),]
  candsMatrix3 <- candsMatrix3*100/11
  candsMatrix3$Gene <- rownames(candsMatrix3)
  
  candsMatrix3.m <- melt(candsMatrix3)
  colnames(candsMatrix3.m) <- c("Gene","Cancer","patients")
  candsMatrix3.m$Gene <- factor(candsMatrix3.m$Gene,levels=as.factor(rownames(candsMatrix3)))
  
  p <- ggplot(candsMatrix3.m) + geom_bar(stat="identity",aes(x=Gene,y=patients,fill=Cancer))
  p <- p + scale_fill_manual(values=colorPalette) + ylab("Accumulated percentage")
  p <- p + theme(text = element_text(size=20),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour="black"))
  
  png(paste0("freq_",annotation,"_accumulatedFrequency.png"), width=1000, height=800)
  print(p)
  graphics.off()
  
}

#most frequent genes, ordered by sum of patients in cancer types
for (annotation in c("All","Driver","d1")){
  ngenes <- 50
  
  if (annotation=="All"){
    selectVector <- matrix(T,nrow(candidateList2),1)
  } else {
    selectVector <- candidateList2$Annotation==annotation
  }
  
  accCandidateList <- ddply(candidateList2[selectVector,],.(GeneId,Symbol), summarise, cumsum=sum(NumPatients))
  accCandidateList <- accCandidateList[order(-accCandidateList$cumsum),]
  
  candsMatrix <- list()
  for (gene in accCandidateList$GeneId[1:ngenes]){
    symbol = accCandidateList$Symbol[accCandidateList$GeneId==gene]
    candsMatrix[[symbol]] <- list()
    candsMatrix[[symbol]][["symbol"]] <- as.character(symbol)
    for (cancer in cancerTypes){
      if (gene %in% candidateList[[cancer]]$GeneId){
        candsMatrix[[symbol]][[cancer]] <- as.numeric(candidateList[[cancer]]$NumPatients[candidateList[[cancer]]$GeneId==gene])
      } else{
        candsMatrix[[symbol]][[cancer]] <- 0
      }
    }
    candsMatrix[[symbol]] <- do.call("cbind", candsMatrix[[symbol]])
  }
  candsMatrix2 <- do.call("rbind", candsMatrix)
  rownames(candsMatrix2) <- candsMatrix2[,c("symbol")]
  
  candsMatrix3 <- as.data.frame(matrix(as.numeric(as.character(candsMatrix2[,!(colnames(candsMatrix2) %in% c("symbol"))])),nrow=ngenes,ncol=11))
  rownames(candsMatrix3) <- rownames(candsMatrix2)
  colnames(candsMatrix3) <- cancerTypes
  candsMatrix3 <- candsMatrix3[order(-rowSums(candsMatrix3)),]
  candsMatrix3 <- candsMatrix3*100/11
  candsMatrix3$Gene <- rownames(candsMatrix3)
  
  candsMatrix3.m <- melt(candsMatrix3)
  colnames(candsMatrix3.m) <- c("Gene","Cancer","patients")
  candsMatrix3.m$Gene <- factor(candsMatrix3.m$Gene,levels=as.factor(rownames(candsMatrix3)))
  
  p <- ggplot(candsMatrix3.m) + geom_bar(stat="identity",aes(x=Gene,y=patients,fill=Cancer))
  p <- p + scale_fill_manual(values=colorPalette) + ylab("Number of patients")
  p <- p + theme(text = element_text(size=20),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour="black"))
  
  png(paste0("freq_",annotation,"_numberOfPatients.png"), width=1000, height=800)
  print(p)
  graphics.off()
  
}


################ CDS STUDY ################
CDS_study <- read.delim("CDS_study.tsv", header=FALSE)
colnames(CDS_study) <- c("Cancer","Analysis","Both","OnlyN","OnlyT","None","Random_Both","Random_OnlyN","Random_OnlyT","Random_None")

################ CDS CHANGE ################
CDS_change <- read.delim("CDS_change.tsv", header=FALSE)
colnames(CDS_change) <- c("Cancer","Analysis","Random_Change","Random_NoChange","Change","NoChange","p","OR")

################ UTR CHANGE ################
UTR_change <- read.delim("UTR_change.tsv", header=FALSE)
colnames(UTR_change) <- c("Cancer","Analysis","Random_Change","Random_NoChange","Change","NoChange","p","OR")

################ DRIVER AFFECTION ################
Driver_D0 <- read.delim("Driver_D0_enrichment.tsv", header=FALSE)
colnames(Driver_D0) <- c("Cancer","Analysis","Driver_Switch","Driver_NoSwitch","NonDriver_Switch","NonDriver_NoSwitch","p","OR")

############# EXONS #################
exons <- read.delim("exons.tsv",header=FALSE)
colnames(exons) <- c("Cancer","Random","Switch","Origin","Type","Length","CDSLength","CDSRelativeSize","Position","KeepOrf")

# cds relative size
cdsRelativeSize <- list()
cdsPosition <- list()
exonLength <- list()
orfChange <- data.frame(cancer=c(),p=c(),oddsRatio=c())
exonOrigin <- data.frame(cancer=c(),p=c(),oddsRatio=c())
for (knsur in cancerTypes){
  cancer.exons = subset(exons,Cancer==knsur)
  
  # cds relative size
  p <- ggplot() + ggtitle(knsur) + theme_bw() + ylab("")
  p <- p + stat_ecdf(data=subset(cancer.exons,Random=="Random"), aes(x=CDSRelativeSize,colour="green"),show_guide = FALSE)
  p <- p + stat_ecdf(data=subset(cancer.exons,Random=="NonRandom"), aes(x=CDSRelativeSize,colour="red"),show_guide = FALSE)
  
  cdsRelativeSize[[knsur]] <- p
  
  # cds position
  p <- ggplot() + ggtitle(knsur) + theme_bw()
  p <- p + stat_ecdf(data=subset(cancer.exons,Random=="Random"), aes(x=Position,colour="green"),show_guide = FALSE)
  p <- p + stat_ecdf(data=subset(cancer.exons,Random=="NonRandom"), aes(x=Position,colour="red"),show_guide = FALSE)
  
  cdsPosition[[knsur]] <- p
  
  # length
  p <- ggplot() + ggtitle(knsur) + theme_bw()
  p <- p + stat_ecdf(data=subset(cancer.exons,Random=="Random"), aes(x=Length,colour="green"),show_guide = FALSE)
  p <- p + stat_ecdf(data=subset(cancer.exons,Random=="NonRandom"), aes(x=Length,colour="red"),show_guide = FALSE)
  
  exonLength[[knsur]] <- p
  
  # orf
  switches <- table(cancer.exons$KeepOrf[cancer.exons$Random=="NonRandom"])
  randomSwitches <- table(cancer.exons$KeepOrf[cancer.exons$Random=="Random"])
  
  cTable <- rbind(switches,randomSwitches)
  f <- fisher.test(cTable)
  this.Data <- data.frame(cancer=knsur,p=f$p.value,oddsRatio=f$estimate)
  orfChange <- rbind(orfChange,this.Data)
  
  # type
  switches <- table(cancer.exons$Type[cancer.exons$Random=="NonRandom" ])
  randomSwitches <- table(cancer.exons$Type[cancer.exons$Random=="Random" ])
  
  cTable <- rbind(switches,randomSwitches)
  f <- fisher.test(cTable)
  this.Data <- data.frame(cancer=knsur,p=f$p.value,oddsRatio=f$estimate)
  exonOrigin <- rbind(exonOrigin,this.Data)
}

png("exon_cds_relative_size.png", width=1000, height=800)
grid.arrange(cdsRelativeSize[["brca"]],cdsRelativeSize[["coad"]],cdsRelativeSize[["hnsc"]],cdsRelativeSize[["kich"]],cdsRelativeSize[["kirc"]],cdsRelativeSize[["kirp"]],cdsRelativeSize[["lihc"]],cdsRelativeSize[["luad"]],cdsRelativeSize[["lusc"]],cdsRelativeSize[["prad"]],cdsRelativeSize[["thca"]])
graphics.off()

png("exon_cds_position.png", width=1000, height=800)
grid.arrange(cdsPosition[["brca"]],cdsPosition[["coad"]],cdsPosition[["hnsc"]],cdsPosition[["kich"]],cdsPosition[["kirc"]],cdsPosition[["kirp"]],cdsPosition[["lihc"]],cdsPosition[["luad"]],cdsPosition[["lusc"]],cdsPosition[["prad"]],cdsPosition[["thca"]])
graphics.off()

png("exonlength.png", width=1000, height=800)
grid.arrange(exonLength[["brca"]],exonLength[["coad"]],exonLength[["hnsc"]],exonLength[["kich"]],exonLength[["kirc"]],exonLength[["kirp"]],exonLength[["lihc"]],exonLength[["luad"]],exonLength[["lusc"]],exonLength[["prad"]],exonLength[["thca"]])
graphics.off()

exonsNew <- read.delim("exons_new.tsv",header=FALSE)
colnames(exonsNew) <- c("Cancer","Random","nExon","tExon","Position")

nrExonsNew <- exonsNew[exonsNew$Random=="NonRandom" ,]
sum(nrExonsNew$nExon%%3==nrExonsNew$tExon%%3)/nrow(nrExonsNew)
sum(nrExonsNew$nExon[nrExonsNew$Position=="BEGINING"]%%3==nrExonsNew$tExon[nrExonsNew$Position=="BEGINING"]%%3)/nrow(nrExonsNew[nrExonsNew$Position=="BEGINING",])
sum(nrExonsNew$nExon[nrExonsNew$Position=="ENDING"]%%3==nrExonsNew$tExon[nrExonsNew$Position=="ENDING"]%%3)/nrow(nrExonsNew[nrExonsNew$Position=="ENDING",])
sum(nrExonsNew$nExon[nrExonsNew$Position=="MIDDLE"]%%3==nrExonsNew$tExon[nrExonsNew$Position=="MIDDLE"]%%3)/nrow(nrExonsNew[nrExonsNew$Position=="MIDDLE",])
################ exons per switch ################

exonsPerSwitch_rel <- read.delim("exonsPerSwitch_relevant.tsv", header=FALSE)
exonsPerSwitch_rel <- cbind(exonsPerSwitch_rel,"Relevant")
colnames(exonsPerSwitch_rel) <- c("Counts","Cancer","Switch","Relevance")

exonsPerSwitch_nrel <- read.delim("exonsPerSwitch_nonrelevant.tsv", header=FALSE)
exonsPerSwitch_nrel <- cbind(exonsPerSwitch_nrel,"NonRelevant")
colnames(exonsPerSwitch_nrel) <- c("Counts","Cancer","Switch","Relevance")

exonsPerSwitch <- rbind(exonsPerSwitch_rel,exonsPerSwitch_nrel)

png("exons_per_switch.png", width=1000, height=800)
par(mar=c(9,3,5,2))
boxplot(Counts~interaction(Relevance,Cancer),data=exonsPerSwitch,outline=F,col=c("gray40","gray90"),las=2)
graphics.off()

barplot(kk1,names.arg=unique(kk$Cancer),beside=T)
barplot(scale(kk1, center=FALSE, scale=colSums(kk1)),names.arg=unique(kk$Cancer),beside=T)

################ CENTRALITY ANALYSIS ################
#centrality_all <- data.frame()
centrality_rel <- data.frame()
centrality_nrel <- data.frame()

for (cancer in cancerTypes){
	basePath <- ""
	#cancer_all <- read.delim(paste0(basePath,cancer,".tIso_length.tsv"), header=FALSE)
	#tiso_length_all <- rbind(tiso_length_all,cbind(cancer_all,cancer))
	
	cancer_rel <- read.delim(paste0(basePath,cancer,".protein_centrality_relevant.tsv"), header=FALSE)
	cancer_rel <- cbind(cancer_rel,"Relevant",cancer)
	colnames(cancer_rel) <- c("Cancer","Degree","Relevance")
	centrality_rel <- rbind(centrality_rel,cancer_rel)

	cancer_nrel <- read.delim(paste0(basePath,cancer,".protein_centrality_nonrelevant.tsv"), header=FALSE)
	cancer_nrel <- cbind(cancer_nrel,"NonRelevant",cancer)
	colnames(cancer_nrel) <- c("Cancer","Degree","Relevance")
	centrality_nrel <- rbind(centrality_nrel,cancer_nrel)
}

centrality <- rbind(centrality_rel,centrality_nrel)

png("centrality.png", width=1000, height=800)
par(mar=c(9,3,5,2))
boxplot(Degree~interaction(Relevance,Cancer),data=centrality,outline=F,col=c("gray40","gray90"),las=2)
graphics.off()

################ NISO ANALYSIS ################
#niso_length_all <- data.frame()
niso_length_rel <- data.frame()
niso_length_nrel <- data.frame()

for (cancer in cancerTypes){
	basePath <- ""
	#cancer_all <- read.delim(paste0(basePath,cancer,".tIso_length.tsv"), header=FALSE)
	#tiso_length_all <- rbind(tiso_length_all,cbind(cancer_all,cancer))
	
	cancer_rel <- read.delim(paste0(basePath,cancer,".nIso_length_relevant.tsv"), header=FALSE)
	cancer_rel <- cbind(cancer_rel,"Relevant",cancer)
	colnames(cancer_rel) <- c("Length","Relevance","Cancer")
	niso_length_rel <- rbind(niso_length_rel,cancer_rel)

	cancer_nrel <- read.delim(paste0(basePath,cancer,".nIso_length_nonrelevant.tsv"), header=FALSE)
	cancer_nrel <- cbind(cancer_nrel,"NonRelevant",cancer)
	colnames(cancer_nrel) <- c("Length","Relevance","Cancer")
  niso_length_nrel <- rbind(niso_length_nrel,cancer_nrel)
}

niso_length <- rbind(niso_length_rel,niso_length_nrel)

png("niso_length.png", width=1000, height=800)
par(mar=c(9,3,5,2))
boxplot(Length~interaction(Relevance,Cancer),data=niso_length,outline=F,col=c("gray40","gray90"),las=2)
graphics.off()

################ TISO ANALYSIS ################
#tiso_length_all <- data.frame()
tiso_length_rel <- data.frame()
tiso_length_nrel <- data.frame()

for (cancer in cancerTypes){
	basePath <- ""
	#cancer_all <- read.delim(paste0(basePath,cancer,".tIso_length.tsv"), header=FALSE)
	#tiso_length_all <- rbind(tiso_length_all,cbind(cancer_all,cancer))
	
	cancer_rel <- read.delim(paste0(basePath,cancer,".tIso_length_relevant.tsv"), header=FALSE)
	cancer_rel <- cbind(cancer_rel,"Relevant",cancer)
	colnames(cancer_rel) <- c("Length","Relevance","Cancer")
	tiso_length_rel <- rbind(tiso_length_rel,cancer_rel)

	cancer_nrel <- read.delim(paste0(basePath,cancer,".tIso_length_nonrelevant.tsv"), header=FALSE)
	cancer_nrel <- cbind(cancer_nrel,"NonRelevant",cancer)
	colnames(cancer_nrel) <- c("Length","Relevance","Cancer")
  	tiso_length_nrel <- rbind(tiso_length_nrel,cancer_nrel)
}

tiso_length <- rbind(tiso_length_rel,tiso_length_nrel)

png("tiso_length.png", width=1000, height=800)
par(mar=c(9,3,5,2))
boxplot(Length~interaction(Relevance,Cancer),data=tiso_length,outline=F,col=c("gray40","gray90"),las=2)
graphics.off()

################ FEATURES AFFECTED ################
Feats_affected_rel <- read.delim("structural_summary_relevant.tsv", header=FALSE)
Feats_affected_rel <- cbind(Feats_affected_rel,"Relevant")
colnames(Feats_affected_rel) <- c("Cancer","Switch","Pfam","PRINTS","ProSitePatterns","IUPREDLong","IUPREDShort","I3D","Driver","Relevance")

Feats_affected_nrel <- read.delim("structural_summary_nonrelevant.tsv", header=FALSE)
Feats_affected_nrel <- cbind(Feats_affected_nrel,"NonRelevant")
colnames(Feats_affected_nrel) <- c("Cancer","Switch","Pfam","PRINTS","ProSitePatterns","IUPREDLong","IUPREDShort","I3D","Driver","Relevance")

Feats_affected <- rbind(Feats_affected_rel,Feats_affected_nrel)
Feats_affected$Driver <- ifelse(Feats_affected$Driver=="True","Driver","NonDriver")

png("ProSitePatterns.png", width=1000, height=800)
par(mar=c(9,3,5,2))
boxplot(ProSitePatterns~interaction(Driver,Cancer),data=Feats_affected[Feats_affected$Relevance=="Relevant" & Feats_affected$ProSitePatterns!=0,],col=c("gray40","gray90"),las=3,outline=F)
graphics.off()

png("Pfam.png", width=1000, height=800)
par(mar=c(9,3,5,2))
boxplot(Pfam~interaction(Driver,Cancer),data=Feats_affected[Feats_affected$Relevance=="Relevant" & Feats_affected$Pfam!=0,],col=c("gray40","gray90"),las=3,outline=F)
graphics.off()

png("PRINTS.png", width=1000, height=800)
par(mar=c(9,3,5,2))
boxplot(PRINTS~interaction(Driver,Cancer),data=Feats_affected[Feats_affected$Relevance=="Relevant" & Feats_affected$PRINTS!=0,],col=c("gray40","gray90"),las=3,outline=F)
graphics.off()

png("I3D.png", width=1000, height=800)
par(mar=c(9,3,5,2))
boxplot(I3D~interaction(Driver,Cancer),data=Feats_affected[Feats_affected$Relevance=="Relevant" & Feats_affected$I3D!=0,],col=c("gray40","gray90"),las=3,outline=F)
graphics.off()

png("IUPREDShort.png", width=1000, height=800)
par(mar=c(9,3,5,2))
boxplot(IUPREDShort~interaction(Driver,Cancer),data=Feats_affected[Feats_affected$Relevance=="Relevant" & Feats_affected$IUPREDShort!=0,],col=c("gray40","gray90"),las=3,outline=F)
graphics.off()

png("IUPREDLong.png", width=1000, height=800)
par(mar=c(9,3,5,2))
boxplot(IUPREDLong~interaction(Driver,Cancer),data=Feats_affected[Feats_affected$Relevance=="Relevant" & Feats_affected$IUPREDLong!=0,],col=c("gray40","gray90"),las=3,outline=F)
graphics.off()

################ STRUCTURAL FEATURES ################
#exons_all <- read.delim("exons.tsv", header=FALSE)
#exons_all <- cbind(exons_all,"All")
#colnames(exons_all) <- c("Cancer","switch","length","type","keepORF","position","Relevance")

structural_features <- read.delim("structural_features.onlyModels.tsv", header=FALSE)
colnames(structural_features) <- c("Cancer","Gene","Symbol","nIso","tIso","Random","Analysis","WhatsHappening","Feature","Driver","ASDriver","DriverType")
nrStructural_features <- structural_features[structural_features$Random == "NonRandom",]

PfamSomethingHappening = structural_features[structural_features$Analysis=='Pfam' & structural_features$WhatsHappening!="Nothing",]
PfamNothingHappening = structural_features[structural_features$Analysis=='Pfam' & structural_features$WhatsHappening!="Nothing",]
a = sum(PfamSomethingHappening$Random == "NonRandom")
b = sum(PfamSomethingHappening$Random == "Random")
c = sum(PfamNothingHappening$Random == "NonRandom")
d = sum(PfamNothingHappening$Random == "Random")

fisher.test(matrix(c(a,b,c,d),nrow=2,ncol=2))

PfamDriver <- sort(table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$Driver==1 & structural_features$WhatsHappening!="Nothing"]),decreasing=TRUE)
PfamDriver <- PfamDriver[PfamDriver>0]
PfamDriver <- 100*PfamDriver/length(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$Driver==1 & structural_features$WhatsHappening!="Nothing"])

PfamNonDriver <-sort(table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$Driver==0 & structural_features$WhatsHappening!="Nothing"]),decreasing=TRUE)
PfamNonDriver <- PfamNonDriver[PfamNonDriver>0]
PfamNonDriver <- 100*PfamNonDriver/length(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$Driver==0 & structural_features$WhatsHappening!="Nothing"])

PfamDriverGained <- table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$Driver==1 & structural_features$WhatsHappening=="Gained in tumor"])
PfamDriverLost <- table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$Driver==1  & structural_features$WhatsHappening=="Lost in tumor"])
PfamDriverDifferences <- PfamDriverGained - PfamDriverLost 
PfamDriverDifferences <- sort(PfamDriverDifferences,decreasing=TRUE)
PfamDriverDifferences <- PfamDriverDifferences[PfamDriverDifferences!=0]
write.table(PfamDriverDifferences,'PfamDriverDifferences.txt',quote=F,col.names=F)

PfamOncogeneGained <- table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$DriverType=="oncogene" & structural_features$WhatsHappening=="Gained in tumor"])
PfamOncogeneLost <- table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$DriverType=="oncogene"  & structural_features$WhatsHappening=="Lost in tumor"])
PfamOncogeneDifferences <- PfamOncogeneGained - PfamOncogeneLost
PfamOncogeneDifferences <- sort(PfamOncogeneDifferences,decreasing=TRUE)
PfamOncogeneDifferences <- PfamOncogeneDifferences[PfamOncogeneDifferences!=0]
write.table(PfamOncogeneDifferences,'PfamOncogeneDifferences.txt',quote=F,col.names=F)

PfamSuppressorGained <- table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$DriverType=="suppressor" & structural_features$WhatsHappening=="Gained in tumor"])
PfamSuppressorLost <- table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$DriverType=="suppressor"  & structural_features$WhatsHappening=="Lost in tumor"])
PfamSuppressorDifferences <- PfamSuppressorGained - PfamSuppressorLost
PfamSuppressorDifferences <- sort(PfamSuppressorDifferences,decreasing=TRUE)
PfamSuppressorDifferences <- PfamSuppressorDifferences[PfamSuppressorDifferences!=0]
write.table(PfamSuppressorDifferences,'PfamSuppressorDifferences.txt',quote=F,col.names=F)


## ---- to review ---- ##

ProSiteDriverGained <- sort(table(structural_features$Feature[structural_features$Analysis=='ProSitePatterns' & structural_features$Driver==1 & structural_features$WhatsHappening=="Gained in tumor"]),decreasing=TRUE)
ProSiteDriverLost <- sort(table(structural_features$Feature[structural_features$Analysis=='ProSitePatterns' & structural_features$Driver==1  & structural_features$WhatsHappening=="Lost in tumor"]),decreasing=TRUE)
ProSiteDifference <- ProSiteDriverGained - ProSiteDriverLost
ProSiteDifference <- sort(ProSiteDifference,decreasing=TRUE)
ProSiteDifference <- ProSiteDifference[ProSiteDifference!=0]
write.table(ProSiteDifference,'ProSiteDifference.txt',quote=F,col.names=F)


ProSitePatternsDriver <-sort(table(structural_features$Feature[structural_features$Analysis=='ProSitePatterns' & structural_features$Driver=="True"]),decreasing=TRUE)
ProSitePatternsDriver <- ProSitePatternsDriver[ProSitePatternsDriver>0]
ProSitePatternsDriver <- 100*ProSitePatternsDriver/length(structural_features$Feature[structural_features$Analysis=='ProSitePatterns' & structural_features$Driver=="True"])

ProSitePatternsNonDriver <-sort(table(structural_features$Feature[structural_features$Analysis=='ProSitePatterns' & structural_features$Driver=="False"]),decreasing=TRUE)
ProSitePatternsNonDriver <- ProSitePatternsNonDriver[ProSitePatternsNonDriver>0]
ProSitePatternsNonDriver <- 100*ProSitePatternsNonDriver/length(structural_features$Feature[structural_features$Analysis=='ProSitePatterns' & structural_features$Driver=="False"])



ProSitePatternsDriverLost <- ProSitePatternsDriverLost[ProSitePatternsDriverLost>0]
ProSitePatternsDriverLost <- 100*ProSitePatternsDriverLost/length(structural_features$Feature[structural_features$Analysis=='ProSitePatterns' & structural_features$Driver=="True"  & structural_features$WhatsHappening=="Lost in tumor"])

write.table(ProSitePatternsDriver,'ProSitePatternsDriver.txt',quote=F,col.names=F)
write.table(ProSitePatternsNonDriver,'ProSitePatternsNonDriver.txt',quote=F,col.names=F)
write.table(ProSitePatternsDriverLost,'ProSitePatternsDriverLost.txt',quote=F,col.names=F)
write.table(ProSitePatternsDriverGained,'ProSitePatternsDriverGained.txt',quote=F,col.names=F)

featureCounts <- by(structural_features, structural_features$Analysis, table(x[,5))

################ DOMAINS ENRICHMENT ################
for (act in c('Pfam','prosite')){
  for (action in c("Lost in tumor","Gained in tumor")){
    uniqStuff <- unique( nrStructural_features[nrStructural_features$Analysis==act & nrStructural_features$WhatsHappening!="Nothing" & nrStructural_features$WhatsHappening==action,c("Feature")] )
    actionTag <- tolower(unlist(strsplit(action,split = " "))[1])
    
    driverTest <- list()
    oncogeneTest <- list()
    suppressorTest <- list()
    for (feat in uniqStuff){
      
      usedFeature <- nrStructural_features[nrStructural_features$WhatsHappening==action & nrStructural_features$Feature==feat,]
      nonUsedFeature <- nrStructural_features[nrStructural_features$WhatsHappening==action | nrStructural_features$Feature!=feat,]
      
      driverTest[[feat]] <- data.frame(Feature=feat, WhatsHappening=action,
                                    FDriver=sum(usedFeature$Driver==1),
                                    FNonDriver=sum(usedFeature$Driver==0),
                                    NFDriver=sum(nonUsedFeature$Driver==1),
                                    NFNonDriver=sum(nonUsedFeature$Driver==0)
      )
      oncogeneTest[[feat]] <- data.frame(Feature=feat, WhatsHappening=action,
                                      FOncogene=sum(usedFeature$DriverType=='oncogene'),
                                      FNonOncogene=sum(usedFeature$DriverType!='oncogene'),
                                      NFOncogene=sum(nonUsedFeature$DriverType=='oncogene'),
                                      NFNonOncogene=sum(nonUsedFeature$DriverType!='oncogene')  )
      
      suppressorTest[[feat]] <- data.frame(Feature=feat, WhatsHappening=action,
                                        FSuppressor=sum(usedFeature$DriverType=='suppressor'),
                                        FNonSuppressor=sum(usedFeature$DriverType!='suppressor'),
                                        NFSuppressor=sum(nonUsedFeature$DriverType=='suppressor'),
                                        NFNonSuppressor=sum(nonUsedFeature$DriverType!='suppressor')  )
      
    }
    
    driverTestDf <- do.call('rbind',driverTest)
    oncogeneTestDf <- do.call('rbind',oncogeneTest)
    suppressorTestDf <- do.call('rbind',suppressorTest)
    
    driverTestDf$p <- apply(driverTestDf,1, function(x){ 
      k <- fisher.test(x=matrix(as.numeric(x[3:6]),nrow=2,ncol=2),alternative="greater")
      k$p.value } )
    driverTestDf$p.adj <- p.adjust(driverTestDf$p)
    write.table(driverTestDf,file=paste0("driverTest.",act,".",actionTag,".tsv"),sep="\t", row.names=F, col.names=F, quote=F)
    
    oncogeneTestDf$p <- apply(oncogeneTestDf,1, function(x){ 
      k <- fisher.test(x=matrix(as.numeric(x[3:6]),nrow=2,ncol=2),alternative="greater")
      k$p.value } )
    oncogeneTestDf$p.adj <- p.adjust(oncogeneTestDf$p)
    write.table(oncogeneTestDf,file=paste0("oncogeneTest.",act,".",actionTag,".tsv"),sep="\t", row.names=F, col.names=F, quote=F)
    
    suppressorTestDf$p <- apply(suppressorTestDf,1, function(x){ 
      k <- fisher.test(x=matrix(as.numeric(x[3:6]),nrow=2,ncol=2),alternative="greater")
      k$p.value } )
    suppressorTestDf$p.adj <- p.adjust(suppressorTestDf$p)
    write.table(suppressorTestDf,file=paste0("suppressorTest.",act,".",actionTag,".tsv"),sep="\t", row.names=F, col.names=F, quote=F)
  }
}

############ ISOFORM LENGTH ############

isoLengths <- read.delim("isoform_length.tsv", header=T)
p <- c()

for (cancer in cancerTypes){
  this.isoLengths <- isoLengths[isoLengths$Cancer==cancer & isoLengths$Random=="NonRandom",]
  t <- t.test(this.isoLengths$nIsoLength,this.isoLengths$tIsoLength)
  p <- c(p,t$p.value)
}

isoLengthsSummary <- ddply(isoLengths[isoLengths$Random=="NonRandom",],.(Cancer,Random), summarise, nIso=mean(nIsoLength),tIso=mean(tIsoLength))
isoLengthsSummary <- cbind(isoLengthsSummary,p)
isoLengthsSummary$Diff <- isoLengthsSummary$nIso - isoLengthsSummary$tIso
isoLengthsSummary <- isoLengthsSummary[,colnames(isoLengthsSummary)!="Random"]

write.table(isoLengthsSummary,file="isoform_length_p.tsv",sep="\t", row.names=F, col.names=F, quote=F)

############ STUDY FEATURES ############
analysis <- c("interpro","anchor","iupred","prosite")
df.features <- data.frame(cancer=c(),analysis=c(),what=c(), p=c())
for (cancer in cancerTypes){
  for (tag in analysis){
    features <- read.delim(paste0(cancer,".",tag,"_analysis.tsv"),header=T)
    random_data <- read.delim(paste0(cancer,".",tag,"_analysis_random.tsv"),header=T)
    
    features <- features[as.logical(features$Significant),]
    random_data <- random_data[as.logical(random_data$Significant),]
    
    switches <- table(features$What)
    randomSwitches <- table(random_data$What)
    counts.features <- rbind(switches,randomSwitches)
    lostVsNothing <- fisher.test(counts.features[,-1])
    gainVsNothing <- fisher.test(counts.features[,-2])
    gainVsLost <- fisher.test(counts.features[,-3])
    
    this.features <- data.frame(cancer=cancer,analysis=tag,what=c("lostVsNothing","gainVsNothing","gainVsLost"), p=c(lostVsNothing$p.value,gainVsNothing$p.value,gainVsLost$p.value),oddsratio=c(lostVsNothing$estimate,gainVsNothing$estimate,gainVsLost$estimate))
    df.features <- rbind(df.features,this.features)
  }
}

################ MUTATIONS ################
muts <- c("all_mutations","functional_mutations")
swis <- c("all_switches","functional_switches")
analysis <- c("gene","geneset","pannegative")
usedSet <- c("","onlyDrivers")

for (a in analysis){
  for (m in muts){
    for (s in swis){
      for (u in usedSet){
        if (u=="onlyDrivers" && a %in% c("gene","pannegative") ){
          next
        }
        
        if (u=="onlyDrivers"){
          analysisFile <- paste0(paste(a,m,s,u,sep="_"),'.txt')
          outFile <- paste(a,m,s,u,'allCancers.txt',sep="_")
        } else {
          analysisFile <- paste0(paste(a,m,s,sep="_"),'.txt')
          outFile <- paste(a,m,s,'allCancers.txt',sep="_")
        }
        
        comp <- read.delim(analysisFile, header=TRUE)
        
        if (a %in% c("gene","pannegative")){
          comp$Hallmark <- ""
        }
        
        comp <- merge(comp,nPatientsDf)
        
        comp$ms[comp$ms==0] <- 0.000001
        comp$s[comp$s==0] <- 0.000001
        big_comp <- ddply(comp,.(Gene,Symbol,Hallmark), summarise, MS=sum(ms),M=sum(m),S=sum(s),N=sum(n),H=-sum(((ms+s)/TotalPatients)/sum((ms+s)/TotalPatients)*log2(((ms+s)/TotalPatients)/sum((ms+s)/TotalPatients)))/log2(length(ms)))
        
        big_comp$MS <- round(big_comp$MS)
        big_comp$S <- round(big_comp$S)
        #big_comp$H[is.na(big_comp$H)] <- 0
        
        big_comp$p_me <- apply(big_comp,1, function(x){ 
          k <- fisher.test(x=matrix(as.numeric(x[4:7]),nrow=2,ncol=2),alternative="less")
          k$p.value } )
        big_comp$p.adj_me <- p.adjust(big_comp$p_me)
        big_comp$p_o <- apply(big_comp,1, function(x){ 
          k <- fisher.test(x=matrix(as.numeric(x[4:7]),nrow=2,ncol=2),alternative="greater")
          k$p.value } )
        big_comp$p.adj_o <- p.adjust(big_comp$p_o)
        big_comp <- big_comp[order(big_comp$p_me),]
        
        if (a %in% c("gene","pannegative")){
          big_comp <- big_comp[,!colnames(big_comp) %in% c("Hallmark")]
        }
        
        write.table(big_comp,outFile,quote=F,col.names=T,sep="\t",row.names=FALSE)
      }
    }
  }
}

aggAnalysis <- list.files(pattern=".+allCancers.txt", full.names=TRUE)
file.copy(aggAnalysis, "/projects_rg/TCGA/users/hector/mutation_comparison",overwrite=T)

# calculate gene score
gn_all <- read.delim(paste("gene","all_mutations","functional_switches",'allCancers.txt',sep="_"), header=TRUE)
gn_fun <- read.delim(paste("gene","functional_mutations","functional_switches",'allCancers.txt',sep="_"), header=TRUE)

gn_merged <- merge(gn_fun,gn_all,by=c("Gene","Symbol"))
gn_merged$Score <- -log10(gn_merged$p_me.x)+log10(gn_merged$p_o.y)

gn_merged$Sign[gn_merged$Score > 0] <- "Mutual exclusion"
gn_merged$Sign[gn_merged$Score < 0] <- "Coincidence"
gn_merged$Sign[gn_merged$Score == 0] <- "Nothing"

gn_merged$NewScore <- ifelse(gn_merged$Score > 0,sqrt(gn_merged$Score),-sqrt(-gn_merged$Score))
gn_merged$NewScore <- ifelse(gn_merged$Score == 0,0,gn_merged$NewScore)
gn_merged <- gn_merged[order(-gn_merged$NewScore),]
gn_merged$Symbol <- factor(gn_merged$Symbol,levels=as.factor(unique(gn_merged$Symbol)))
gn_merged$Order=1:nrow(gn_merged)

#gn_merged <- gn_merged[abs(gn_merged$Score)>0,]
write.table(gn_merged[,c("Gene","Symbol","Score","Sign" )],'scoredMEGenes.txt',quote=F,col.names=T,sep="\t",row.names=FALSE)
file.copy('scoredMEGenes.txt', "/projects_rg/TCGA/users/hector/mutation_comparison",overwrite=T)

plotThreshold <- 2
gn_merged_plot <- gn_merged[!abs(gn_merged$NewScore)>plotThreshold,]

p <- ggplot(gn_merged_plot, aes(x=Order,y=NewScore,fill = Sign)) + geom_area(alpha=0.75)
p <- p + smartas_theme() + labs(title="Genes ranked by score", x="Genes", y="signed sqrt(Score)")
p <- p + geom_hline(yintercept=0, size=0.4, color="black") + scale_y_continuous(breaks=seq(-plotThreshold,plotThreshold,by=1))
p <- p + scale_fill_manual(values=c("Mutual exclusion"="red","Nothing"="black","Coincidence"="darkblue"))
p <- p + theme(axis.text.x=element_blank())
ggsave("meScores.png",p)
file.copy("meScores.png","/projects_rg/TCGA/users/hector/mutation_comparison",overwrite=T)

# calculate geneset score
gnset_all <- read.delim(paste("geneset","all_mutations","functional_switches",'onlyDrivers','allCancers.txt',sep="_"), header=TRUE)
gnset_fun <- read.delim(paste("geneset","functional_mutations","functional_switches",'onlyDrivers','allCancers.txt',sep="_"), header=TRUE)

gnset_merged <- merge(gnset_fun,gnset_all,by=c("Gene","Symbol","Hallmark"))
gnset_merged$Score <- -log10(gnset_merged$p_me.x)+log10(gnset_merged$p_o.y)

gnset_merged$Sign[gnset_merged$Score > 0] <- "Mutual exclusion"
gnset_merged$Sign[gnset_merged$Score < 0] <- "Coincidence"
gnset_merged$Sign[gnset_merged$Score == 0] <- "Nothing"

gnset_merged$NewScore <- ifelse(gnset_merged$Score > 0,sqrt(gnset_merged$Score),-sqrt(-gnset_merged$Score))
gnset_merged$NewScore <- ifelse(gnset_merged$Score == 0,0,gnset_merged$NewScore)
gnset_merged <- gnset_merged[order(-gnset_merged$NewScore),]
gnset_merged$Symbol <- factor(gnset_merged$Symbol,levels=as.factor(unique(gnset_merged$Symbol)))
gnset_merged$Order=1:nrow(gnset_merged)

write.table(gnset_merged[,c("Gene","Symbol","Hallmark","Score","Sign" )],'scoredMEGenesets.txt',quote=F,col.names=T,sep="\t",row.names=FALSE)
file.copy('scoredMEGenesets.txt', "/projects_rg/TCGA/users/hector/mutation_comparison",overwrite=T)

plotThreshold <- 6
gnset_merged_plot <- gnset_merged[!abs(gnset_merged$NewScore)>plotThreshold,]

p <- ggplot(gnset_merged_plot, aes(x=Order,y=NewScore,fill = Sign)) + geom_area(alpha=0.75)
p <- p + smartas_theme() + labs(title="Genes ranked by score", x="Genes", y="signed sqrt(Score)")
p <- p + geom_hline(yintercept=0, size=0.4, color="black") + scale_y_continuous(breaks=seq(-plotThreshold,plotThreshold,by=1))
p <- p + scale_fill_manual(values=c("Mutual exclusion"="red","Nothing"="black","Coincidence"="darkblue"))
p <- p + theme(axis.text.x=element_blank())
ggsave("meGenesetScores.png",p)
file.copy("meGenesetScores.png","/projects_rg/TCGA/users/hector/mutation_comparison",overwrite=T)

################ MUTATION OVERLAP ################
mutation_feature_overlap <- read.delim("mutation_switch_feature_overlap.txt", header=TRUE)

pvals <- apply(mutation_feature_overlap[,c("MutationsInFeature","TotalMutations","FeatureSize")],1, function(x){ 
  if (x[2]==0){
    p <- 1
  } else {
    p <- binom.test(x[1],x[2],x[3]/100,"greater")
    p <- p$p.value
  }
  p
  } )
adjpvalues <- p.adjust(pvals)

mutation_feature_overlap$p <- pvals
mutation_feature_overlap$adjp <- adjpvalues

mutation_feature_overlap_sum <- ddply(mutation_feature_overlap,.(Gene,Symbol,Normal_transcript,Tumor_transcript,What,FeatureType,Feature,Driver,FeatureSize), summarise, inMut=sum(MutationsInFeature), totalMut=sum(TotalMutations) )
mutation_feature_overlap_sum$Ratio = 100 * mutation_feature_overlap_sum$inMut/mutation_feature_overlap_sum$totalMut

pvals <- apply(mutation_feature_overlap_sum[,c("inMut","totalMut","FeatureSize")],1, function(x){ 
  if (x[2]==0){
    p <- 1
  } else {
    p <- binom.test(x[1],x[2],x[3]/100,"greater")
    p <- p$p.value
  }
  p
} )
adjpvalues <- p.adjust(pvals)

mutation_feature_overlap_sum$p_mutation_overlap <- pvals
mutation_feature_overlap_sum$adjp_mutation_overlap <- adjpvalues

gn_fun <- read.delim(paste("gene","functional_mutations","functional_switches",'allCancers.txt',sep="_"), header=TRUE)
gn_fun <- gn_fun[,c("Gene","Symbol","p_me")]

mutation_feature_overlap_sum <- merge(mutation_feature_overlap_sum,gn_fun,by=c("Gene","Symbol"))

# drivers
p <- ggplot(subset(mutation_feature_overlap_sum, Driver=="True"),aes(-log10(p_me),-log10(p_mutation_overlap))) + geom_point()
p <- p + geom_point(data=subset(mutation_feature_overlap_sum,Driver=="True" & p_mutation_overlap < 0.5 & p_me < 0.5),aes(-log10(p_me),-log10(p_mutation_overlap),color=Symbol))
p <- p + theme_minimal()
p <- direct.label(p)

# all
p <- ggplot(mutation_feature_overlap_sum,aes(-log10(p_me),-log10(p_mutation_overlap))) + geom_point()
p <- p + geom_point(data=subset(mutation_feature_overlap_sum, p_mutation_overlap < 0.2 & p_me < 0.2),aes(-log10(p_me),-log10(p_mutation_overlap),color=Symbol))
p <- p + theme_minimal()
p <- direct.label(p)

# filter by positive MEScore
kk <- mutation_feature_overlap_sum
kk <- kk[kk$Gene %in% geneScores$Gene[geneScores$Score>0],]

p <- ggplot(kk,aes(-log10(p_me),-log10(p_mutation_overlap))) + geom_point()
p <- p + geom_point(data=subset(kk, p_mutation_overlap < 0.2 & p_me < 0.2),aes(-log10(p_me),-log10(p_mutation_overlap),color=Symbol))
p <- p + theme_minimal()
p <- direct.label(p)


p <- p + geom_text(data=subset(mutation_feature_overlap_sum,MutationScore >= 2 & MEScore > 0), aes(MEScore,MutationScore,label=Symbol))

annotate("text", x=c(-0.2,0.9), y=c(33,32), label = c("my label", "label 2"))

mutation_feature_overlap_sum <- mutation_feature_overlap_sum[order(-mutation_feature_overlap_sum$Ratio),]

write.table(mutation_feature_overlap_sum,'mutation_switch_feature_overlap_allCancers.txt',quote=F,col.names=T,sep="\t",row.names=FALSE)
file.copy('mutation_switch_feature_overlap_allCancers.txt',"/projects_rg/TCGA/users/hector/mutation_comparison",overwrite=T)

################ I3D ################
i3d <- read.delim("i3d_broken.tsv", header=TRUE,row.names=NULL)
colnames(i3d) <- c("Cancer","Symbol","Annotation","Partner","PartnerAnnotation","WhatsHappening","InteractionAffection","SequenceCover","PartnerCover","Gene","normalTranscript","tumorTranscript","Uniprot","PartnerUniprot","SequenceInformation","IsoformSpecific","PDB","pymolCommands")

i3d$Structure <- NULL
i3d$Structure[grepl('-MDL-',i3d$PDB)] <- "MODEL"
i3d$Structure[grepl('-EXP-',i3d$PDB)] <- "EXPERIMENTAL"
i3d$Structure[grepl('-MDD-',i3d$PDB)] <- "DOMAIN_MODEL"

i3d$DriverRelationship <- "No driver involved"
i3d$DriverRelationship[grep("Driver",i3d$Annotation)] <- "Driver"
i3d$DriverRelationship[grep("Driver",i3d$PartnerAnnotation)] <- "d1"
i3d$DriverRelationship[Reduce(intersect,  list(grep("Driver",i3d$Annotation),grep("Driver",i3d$PartnerAnnotation)))] <- "Driver-Driver"

p <- ggplot(i3d,aes(InteractionAffection,(SequenceCover+PartnerCover)/2)) + geom_point(aes(shape=DriverRelationship))
p <- p + geom_point(data=subset(i3d,InteractionAffection > 60 & (SequenceCover+PartnerCover)/2 > 60),aes(InteractionAffection,(SequenceCover+PartnerCover)/2,color=Symbol,shape=DriverRelationship ))
p <- p + theme_minimal()
p <- direct.label(p)

# take only interactions where both partners are sufficiently
# described and the interaction is actually affected
i3d_complete <- i3d[i3d$SequenceCover>80 & i3d$InteractionAffection>0 & i3d$PartnerCover>80,!colnames(i3d) %in% c("SequenceInformation","IsoformSpecific","PDB")]
i3d_complete <- i3d_complete[order(-i3d_complete$InteractionAffection),]

# unique interactions
length(unique(i3d$PDB[i3d$SequenceCover>80 & i3d$InteractionAffection>0 & i3d$PartnerCover>80]))
table(unique(i3d[i3d$SequenceCover>80 & i3d$InteractionAffection>0 & i3d$PartnerCover>80,c("PDB","Structure","InteractionAffection")])$Structure)

df <- unique(i3d_complete[,c("Symbol","Partner","normalTranscript","tumorTranscript","InteractionAffection","Annotation","PartnerAnnotation")])
nrow(df)

p <- ggplot(df, aes(InteractionAffection)) + geom_histogram(binwidth=5,fill="#c0392b", alpha=0.75)
p <- p + smartas_theme() + labs(title="Distribution of the percentage of the interaction affected", x="Percentage of the interaction affected", y="Frequency")
p <- p + geom_hline(yintercept=0, size=0.4, color="black") + scale_x_continuous(breaks=seq(0,100, by=5))
ggsave("i3d_interaction_affected.png",p)

freqIntx = as.data.frame(table(df[,c("Annotation","PartnerAnnotation")]))
freqIntx <- freqIntx[order(-freqIntx$Freq),]
write.table(freqIntx,'intxAnnotationFrequency.txt',quote=F,col.names=T,sep="\t",row.names=FALSE)

nrow(i3d_complete)
table(i3d_complete$WhatsHappening)

scoredTable <- merge(i3d,gn_merged[,colnames(gn_merged) %in% c("Gene","Score")],by=c("Gene"))
scoredTable <- rename(scoredTable,c("Score"="GeneScore"))
scoredTable <- merge(scoredTable,gnset_merged[,colnames(gnset_merged) %in% c("Gene","Score")],by=c("Gene"))
scoredTable <- rename(scoredTable,c("Score"="GenesetScore"))

filtScoredTable <- unique(scoredTable[scoredTable$SequenceCover>50 & scoredTable$PartnerCover>50 & scoredTable$InteractionAffection > 10,c("Symbol","Partner","InteractionAffection","Annotation","PartnerAnnotation","SequenceCover","PartnerCover","WhatsHappening","GeneScore","GenesetScore","Structure")])
filtScoredTable <- filtScoredTable[order(-filtScoredTable$GeneScore),]

filtScoredTable <- filtScoredTable[order(-filtScoredTable$GenesetScore),]

################ NEIGHBORHOODS ################
genesetFiles <- c("canonical_pathways","hallmarks","go_biological_process")

for (file in genesetFiles){
  sets_raw <- read.delim(paste0(file,'.txt'),header=TRUE)
  set_counts <- table(sets_raw$GeneSet[sets_raw$qval<0.05])
  set_counts <- set_counts[order(-set_counts)]
  
  write.table(as.data.frame(set_counts),paste0(file,'_counts.txt'),quote=F,col.names=F,sep="\t")
}

kk <- read.delim(paste0("hallmark_info/","HALLMARK_ADIPOGENESIS_onlyDrivers.tsv"))
