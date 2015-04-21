library(plyr)
library(ggplot2)
library(reshape2)

cancerTypes <- c("brca","coad","hnsc","kich","kirc","kirp","lihc","luad","lusc","prad","thca")
workingDir <- "/genomics/users/hector/TCGA_analysis"
setwd(workingDir)

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
accCandidateList <- ddply(candidateList2,.(GeneId,Symbol), summarise, cumsum=sum(Percentage))
accCandidateList <- accCandidateList[order(-accCandidateList$cumsum),]

candsMatrix <- list()
for (gene in accCandidateList$GeneId[1:20]){
  symbol = accCandidateList$Symbol[accCandidateList$GeneId==gene]
  candsMatrix[[symbol]] <- list()
  candsMatrix[[symbol]][["symbol"]] <- as.character(symbol)
  for (cancer in cancerTypes){
    if (gene %in% candidateList[[cancer]]$GeneId){
      candsMatrix[[symbol]][[cancer]] <- candidateList[[cancer]]$Percentage[candidateList[[cancer]]$GeneId==gene]
    } else{
      candsMatrix[[symbol]][[cancer]] <- 0
    }
  }
  candsMatrix[[symbol]] <- do.call("cbind", candsMatrix[[symbol]])
}
candsMatrix2 <- do.call("rbind", candsMatrix)
rownames(candsMatrix2) <- candsMatrix2[,c("symbol")]
candsMatrix2 <- candsMatrix2[,-c("symbol")]
candsMatrix2.m <- melt(candsMatrix2)

colors <- c("red","blue","red","blue","red","blue","red","blue","red","blue","red")

ggplot(candsMatrix2.m,aes(x=Var1,y=value,fill=colors)) + geom_bar(stat="identity")

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
colnames(structural_features) <- c("Cancer","Gene","Symbol","nIso","tIso","Random","Analysis","WhatsHappenning","Feature","Driver","ASDriver","DriverType")

PfamSomethingHappening = structural_features[structural_features$Analysis=='Pfam' & structural_features$WhatsHappenning!="Nothing",]
PfamNothingHappening = structural_features[structural_features$Analysis=='Pfam' & structural_features$WhatsHappenning!="Nothing",]
a = sum(PfamSomethingHappening$Random == "NonRandom")
b = sum(PfamSomethingHappening$Random == "Random")
c = sum(PfamNothingHappening$Random == "NonRandom")
d = sum(PfamNothingHappening$Random == "Random")

fisher.test(matrix(c(a,b,c,d),nrow=2,ncol=2))

PfamDriver <- sort(table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$Driver==1 & structural_features$WhatsHappenning!="Nothing"]),decreasing=TRUE)
PfamDriver <- PfamDriver[PfamDriver>0]
PfamDriver <- 100*PfamDriver/length(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$Driver==1 & structural_features$WhatsHappenning!="Nothing"])

PfamNonDriver <-sort(table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$Driver==0 & structural_features$WhatsHappenning!="Nothing"]),decreasing=TRUE)
PfamNonDriver <- PfamNonDriver[PfamNonDriver>0]
PfamNonDriver <- 100*PfamNonDriver/length(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$Driver==0 & structural_features$WhatsHappenning!="Nothing"])

PfamDriverGained <- table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$Driver==1 & structural_features$WhatsHappenning=="Gained in tumor"])
PfamDriverLost <- table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$Driver==1  & structural_features$WhatsHappenning=="Lost in tumor"])
PfamDriverDifferences <- PfamDriverGained - PfamDriverLost 
PfamDriverDifferences <- sort(PfamDriverDifferences,decreasing=TRUE)
PfamDriverDifferences <- PfamDriverDifferences[PfamDriverDifferences!=0]
write.table(PfamDriverDifferences,'PfamDriverDifferences.txt',quote=F,col.names=F)

PfamOncogeneGained <- table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$DriverType=="oncogene" & structural_features$WhatsHappenning=="Gained in tumor"])
PfamOncogeneLost <- table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$DriverType=="oncogene"  & structural_features$WhatsHappenning=="Lost in tumor"])
PfamOncogeneDifferences <- PfamOncogeneGained - PfamOncogeneLost
PfamOncogeneDifferences <- sort(PfamOncogeneDifferences,decreasing=TRUE)
PfamOncogeneDifferences <- PfamOncogeneDifferences[PfamOncogeneDifferences!=0]
write.table(PfamOncogeneDifferences,'PfamOncogeneDifferences.txt',quote=F,col.names=F)

PfamSuppressorGained <- table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$DriverType=="suppressor" & structural_features$WhatsHappenning=="Gained in tumor"])
PfamSuppressorLost <- table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$DriverType=="suppressor"  & structural_features$WhatsHappenning=="Lost in tumor"])
PfamSuppressorDifferences <- PfamSuppressorGained - PfamSuppressorLost
PfamSuppressorDifferences <- sort(PfamSuppressorDifferences,decreasing=TRUE)
PfamSuppressorDifferences <- PfamSuppressorDifferences[PfamSuppressorDifferences!=0]
write.table(PfamSuppressorDifferences,'PfamSuppressorDifferences.txt',quote=F,col.names=F)


## ---- to review ---- ##

ProSiteDriverGained <- sort(table(structural_features$Feature[structural_features$Analysis=='ProSitePatterns' & structural_features$Driver==1 & structural_features$WhatsHappenning=="Gained in tumor"]),decreasing=TRUE)
ProSiteDriverLost <- sort(table(structural_features$Feature[structural_features$Analysis=='ProSitePatterns' & structural_features$Driver==1  & structural_features$WhatsHappenning=="Lost in tumor"]),decreasing=TRUE)
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
ProSitePatternsDriverLost <- 100*ProSitePatternsDriverLost/length(structural_features$Feature[structural_features$Analysis=='ProSitePatterns' & structural_features$Driver=="True"  & structural_features$Action=="Lost in tumor"])

write.table(ProSitePatternsDriver,'ProSitePatternsDriver.txt',quote=F,col.names=F)
write.table(ProSitePatternsNonDriver,'ProSitePatternsNonDriver.txt',quote=F,col.names=F)
write.table(ProSitePatternsDriverLost,'ProSitePatternsDriverLost.txt',quote=F,col.names=F)
write.table(ProSitePatternsDriverGained,'ProSitePatternsDriverGained.txt',quote=F,col.names=F)

featureCounts <- by(structural_features, structural_features$Analysis, table(x[,5))

################ DOMAINS ENRICHMENT ################
domains <- read.delim("~/Desktop/structural_features.onlyModels.tsv", header=FALSE)
colnames(domains) <- c("Cancer","Gene","Analysis","Symbol","nTx","tTx","Action","Feature","Relevant","Model","Noise","Driver","ASDriver","DriverType")
# once the new feature selection is run, this shouldn't be a problem
domains <- domains[!is.na(domains$Feature),]

uniqStuff <- unique( domains[domains$Analysis=='Pfam',c("Feature","Action")] )

driverTest <- list()
oncogeneTest <- list()
suppressorTest <- list()
for (i in 1:nrow(uniqStuff)){
  feat <- uniqStuff[i,1]
  action <- uniqStuff[i,2]
  
  
  
  usedFeature <- domains[domains$Action==action & domains$Feature==feat,]
  nonUsedFeature <- domains[domains$Action!=action | domains$Feature!=feat,]
  
  driverTest[[i]] <- data.frame(Feature=feat, Action=action,
                                FDriver=sum(usedFeature$Driver==1),
                                FNonDriver=sum(usedFeature$Driver==0),
                                NFDriver=sum(nonUsedFeature$Driver==1),
                                NFNonDriver=sum(nonUsedFeature$Driver==0)
  )
  oncogeneTest[[i]] <- data.frame(Feature=feat, Action=action,
                                  FOncogene=sum(usedFeature$DriverType=='oncogene'),
                                  FNonOncogene=sum(usedFeature$DriverType!='oncogene'),
                                  NFOncogene=sum(nonUsedFeature$DriverType=='oncogene'),
                                  NFNonOncogene=sum(nonUsedFeature$DriverType!='oncogene')  )
  
  suppressorTest[[i]] <- data.frame(Feature=feat, Action=action,
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
write.table(driverTestDf,file="~/Desktop/driverTest.tsv",sep="\t", row.names=F, col.names=F, quote=F)

oncogeneTestDf$p <- apply(oncogeneTestDf,1, function(x){ 
  k <- fisher.test(x=matrix(as.numeric(x[3:6]),nrow=2,ncol=2),alternative="greater")
  k$p.value } )
oncogeneTestDf$p.adj <- p.adjust(oncogeneTestDf$p)
write.table(oncogeneTestDf,file="~/Desktop/oncogeneTest.tsv",sep="\t", row.names=F, col.names=F, quote=F)

suppressorTestDf$p <- apply(suppressorTestDf,1, function(x){ 
  k <- fisher.test(x=matrix(as.numeric(x[3:6]),nrow=2,ncol=2),alternative="greater")
  k$p.value } )
suppressorTestDf$p.adj <- p.adjust(suppressorTestDf$p)
write.table(suppressorTestDf,file="~/Desktop/suppressorTest.tsv",sep="\t", row.names=F, col.names=F, quote=F)

############ STUDY FEATURES ############
analysis <- c("interpro","anchor","iupred","prosite")
df.features <- data.frame(cancer=c(),analysis=c(),what=c(), p=c())
for (cancer in cancerTypes){
  for (tag in analysis){
    features <- read.delim(paste0(cancer,".",tag,"_analysis.tsv"))
    random_data <- read.delim(paste0(cancer,".",tag,"_analysis_random.tsv"))
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