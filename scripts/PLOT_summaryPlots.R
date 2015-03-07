library(plyr)

cancerTypes <- c("brca","coad","hnsc","kich","kirc","kirp","lihc","luad","lusc","prad","thca")

################ CANDIDATES STUDY ################
for (cancer in cancerTypes){
	candidateList <- read.delim(paste0("~/SmartAS/testResults/TCGA/",cancer,"/candidateList_v3.tsv"))
	plotcmd <- "plot(log(candidateList_v3$Patient_percentage),-log(candidateList_v3$Sensitivity))"
	save.plot(plotcmd, file=paste0(cancer,"_patient-sensitivity.png"), dir=getwd(),w=1000, h=1000, format="png")
	plotcmd <- "plot(sqrt(candidateList_v3$Patient_percentage),candidateList_v3$Precision)"
	save.plot(plotcmd, file=paste0(cancer,"_patient-precision.png"), dir=getwd(),w=1000, h=1000, format="png")
}

################ CDS STUDY ################
#CDS_study_all <- read.delim("/genomics/users/hector/TCGA_analysis/CDS_study.tsv", header=FALSE)
#CDS_study_all <- cbind(CDS_study_all,"All")
#colnames(CDS_study_all) <- c("Cancer","Analysis","Both","Only_nIso","Only_tIso","None","Total","Relevance")

CDS_study_rel <- read.delim("/genomics/users/hector/TCGA_analysis/CDS_study_relevant.tsv", header=FALSE)
CDS_study_rel <- cbind(CDS_study_rel,"Relevant")
colnames(CDS_study_rel) <- c("Cancer","Analysis","Both","Only_nIso","Only_tIso","None","Total","Relevance")

CDS_study_nrel <- read.delim("/genomics/users/hector/TCGA_analysis/CDS_study_nonrelevant.tsv", header=FALSE)
CDS_study_nrel <- cbind(CDS_study_nrel,"NonRelevant")
colnames(CDS_study_nrel) <- c("Cancer","Analysis","Both","Only_nIso","Only_tIso","None","Total","Relevance")

CDS_study <- rbind(CDS_study_rel,CDS_study_nrel)

################ CDS CHANGE ################
#CDS_change_all <- read.delim("/genomics/users/hector/TCGA_analysis/CDS_change.tsv", header=FALSE)
#CDS_change_all <- cbind(CDS_change_all,"All")
#colnames(CDS_change_all) <- c("Cancer","Analysis","Yes","No","Total","Relevance")

CDS_change_rel <- read.delim("/genomics/users/hector/TCGA_analysis/CDS_change_relevant.tsv", header=FALSE)
CDS_change_rel <- cbind(CDS_change_rel,"Relevant")
colnames(CDS_change_rel) <- c("Cancer","Analysis","Yes","No","Total","Relevance")

CDS_change_nrel <- read.delim("/genomics/users/hector/TCGA_analysis/CDS_change_nonrelevant.tsv", header=FALSE)
CDS_change_nrel <- cbind(CDS_change_nrel,"NonRelevant")
colnames(CDS_change_nrel) <- c("Cancer","Analysis","Yes","No","Total","Relevance")

CDS_change <- rbind(CDS_change_rel,CDS_change_nrel)
CDS_change$Ratio <- CDS_change$Yes/CDS_change$Total

plotTable_CDS_change <- data.frame()

for (cancer in unique(CDS_change$Cancer)){
  plotTable_CDS_change <- rbind(plotTable_CDS_change, t(CDS_change$Ratio[CDS_change$Cancer==cancer]))
}
colnames(plotTable_CDS_change) <- c("Rel","NonRel")

png("CDS_change.png", width=1000, height=800)
barplot(t(plotTable_CDS_change),beside=T,names.arg=unique(CDS_change$Cancer))
graphics.off()

################ UTR CHANGE ################
#UTR_change_all <- read.delim("/genomics/users/hector/TCGA_analysis/UTR_change.tsv", header=FALSE)
#UTR_change_all <- cbind(UTR_change_all,"All")
#colnames(UTR_change_all) <- c("Cancer","Analysis","Yes","No","Total","Relevance")

UTR_change_rel <- read.delim("/genomics/users/hector/TCGA_analysis/UTR_change_relevant.tsv", header=FALSE)
UTR_change_rel <- cbind(UTR_change_rel,"Relevant")
colnames(UTR_change_rel) <- c("Cancer","Analysis","Yes","No","Total","Relevance")

UTR_change_nrel <- read.delim("/genomics/users/hector/TCGA_analysis/UTR_change_nonrelevant.tsv", header=FALSE)
UTR_change_nrel <- cbind(UTR_change_nrel,"NonRelevant")
colnames(UTR_change_nrel) <- c("Cancer","Analysis","Yes","No","Total","Relevance")

UTR_change <- rbind(UTR_change_rel,UTR_change_nrel)
UTR_change$Ratio <- UTR_change$Yes/UTR_change$Total

plotTable_UTR_change <- data.frame()

for (cancer in unique(UTR_change$Cancer)){
  plotTable_UTR_change <- rbind(plotTable_UTR_change, t(UTR_change$Ratio[UTR_change$Cancer==cancer]))
}
colnames(plotTable_UTR_change) <- c("Rel","NonRel")

png("UTR_change.png", width=1000, height=800)
barplot(t(plotTable_UTR_change),beside=T,names.arg=unique(UTR_change$Cancer))
graphics.off()

################ RELEVANCE ANALYSIS ################
Relevant <- read.delim("/genomics/users/hector/TCGA_analysis/Relevant.tsv", header=FALSE)
Relevant <- cbind(Relevant,"All")
colnames(Relevant) <- c("Cancer","Analysis","Yes","No","Total","Relevance")

Relevant$Ratio <- Relevant$Yes/Relevant$Total

png("Relevant.png", width=1000, height=800)
barplot(Relevant$Ratio,names.arg=unique(Relevant$Cancer))
graphics.off()

################ LOOPS ################
#structural_loops_all <- read.delim("/genomics/users/hector/TCGA_analysis/structural_loops.tsv", header=FALSE)
#structural_loops_all <- cbind(structural_loops_all,"All")
#colnames(structural_loops_all) <- c("Cancer","Different","Same","Only_nIso","Only_tIso","None","Total","Relevance")

structural_loops_rel <- read.delim("/genomics/users/hector/TCGA_analysis/structural_loops_relevant.tsv", header=FALSE)
structural_loops_rel <- cbind(structural_loops_rel,"Relevant")
colnames(structural_loops_rel) <- c("Cancer","Different","Same","Only_nIso","Only_tIso","None","Total","Relevance")

structural_loops_nrel <- read.delim("/genomics/users/hector/TCGA_analysis/structural_loops_nonrelevant.tsv", header=FALSE)
structural_loops_nrel <- cbind(structural_loops_nrel,"NonRelevant")
colnames(structural_loops_nrel) <- c("Cancer","Different","Same","Only_nIso","Only_tIso","None","Total","Relevance")

structural_loops <- rbind(structural_loops_rel,structural_loops_nrel)

################ DRIVER AFFECTION ################
#Driver_affection_all <- read.delim("/genomics/users/hector/TCGA_analysis/Driver_affection.tsv", header=FALSE)
#Driver_affection_all <- cbind(Driver_affection_all,"All")
#colnames(Driver_affection_all) <- c("Cancer","Analysis","Yes","No","Total","Relevance")

Driver_affection_rel <- read.delim("/genomics/users/hector/TCGA_analysis/Driver_affection_relevant.tsv", header=FALSE)
Driver_affection_rel <- cbind(Driver_affection_rel,"Relevant")
colnames(Driver_affection_rel) <- c("Cancer","Analysis","Yes","No","Total","Relevance")

Driver_affection_nrel <- read.delim("/genomics/users/hector/TCGA_analysis/Driver_affection_nonrelevant.tsv", header=FALSE)
Driver_affection_nrel <- cbind(Driver_affection_nrel,"NonRelevant")
colnames(Driver_affection_nrel) <- c("Cancer","Analysis","Yes","No","Total","Relevance")

Driver_affection <- rbind(Driver_affection_rel,Driver_affection_nrel)
Driver_affection$Ratio <- Driver_affection$Yes/Driver_affection$Total

plotTable_Driver_affection <- data.frame()

for (cancer in unique(Driver_affection$Cancer)){
  plotTable_Driver_affection <- rbind(plotTable_Driver_affection, t(Driver_affection$Ratio[Driver_affection$Cancer==cancer]))
}
colnames(plotTable_Driver_affection) <- c("Rel","NonRel")

png("Driver_affection.png", width=1000, height=800)
barplot(t(plotTable_Driver_affection),beside=T,names.arg=unique(Driver_affection$Cancer))
graphics.off()

################ EXON ANALYSIS ################
#exons_all <- read.delim("/genomics/users/hector/TCGA_analysis/exons.tsv", header=FALSE)
#exons_all <- cbind(exons_all,"All")
#colnames(exons_all) <- c("Cancer","switch","length","type","keepORF","position","Relevance")

exons_rel <- read.delim("/genomics/users/hector/TCGA_analysis/exons_relevant.tsv", header=FALSE)
exons_rel <- cbind(exons_rel,"Relevant")
colnames(exons_rel) <- c("Cancer","Switch","length","type","keepORF","position","Relevance")

exons_nrel <- read.delim("/genomics/users/hector/TCGA_analysis/exons_nonrelevant.tsv", header=FALSE)
exons_nrel <- cbind(exons_nrel,"NonRelevant")
colnames(exons_nrel) <- c("Cancer","Switch","length","type","keepORF","position","Relevance")

exons <- rbind(exons_rel,exons_nrel)

png("exon_length.png", width=1000, height=800)
par(mar=c(9,3,5,2))
boxplot(length~interaction(Relevance,Cancer),data=exons,outline=F,col=c("gray40","gray90"),las=2)
graphics.off()

################ exons per switch ################

exonsPerSwitch_rel <- read.delim("/genomics/users/hector/TCGA_analysis/exonsPerSwitch_relevant.tsv", header=FALSE)
exonsPerSwitch_rel <- cbind(exonsPerSwitch_rel,"Relevant")
colnames(exonsPerSwitch_rel) <- c("Counts","Cancer","Switch","Relevance")

exonsPerSwitch_nrel <- read.delim("/genomics/users/hector/TCGA_analysis/exonsPerSwitch_nonrelevant.tsv", header=FALSE)
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
	basePath <- "/genomics/users/hector/TCGA_analysis/"
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
	basePath <- "/genomics/users/hector/TCGA_analysis/"
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
	basePath <- "/genomics/users/hector/TCGA_analysis/"
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
Feats_affected_rel <- read.delim("/genomics/users/hector/TCGA_analysis/structural_summary_relevant.tsv", header=FALSE)
Feats_affected_rel <- cbind(Feats_affected_rel,"Relevant")
colnames(Feats_affected_rel) <- c("Cancer","Switch","Pfam","PRINTS","ProSitePatterns","IUPREDLong","IUPREDShort","I3D","Driver","Relevance")

Feats_affected_nrel <- read.delim("/genomics/users/hector/TCGA_analysis/structural_summary_nonrelevant.tsv", header=FALSE)
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
#exons_all <- read.delim("/genomics/users/hector/TCGA_analysis/exons.tsv", header=FALSE)
#exons_all <- cbind(exons_all,"All")
#colnames(exons_all) <- c("Cancer","switch","length","type","keepORF","position","Relevance")

structural_features_rel <- read.delim("/genomics/users/hector/TCGA_analysis/structural_features_relevant.tsv", header=FALSE)
structural_features_rel <- cbind(structural_features_rel,"Relevant")
colnames(structural_features_rel) <- c("Cancer","Switch","Analysis","Action","Feature","Driver","Relevance")

structural_features_nrel <- read.delim("/genomics/users/hector/TCGA_analysis/structural_features_nonrelevant.tsv", header=FALSE)
structural_features_nrel <- cbind(structural_features_nrel,"NonRelevant")
colnames(structural_features_nrel) <- c("Cancer","Switch","Analysis","Action","Feature","Driver","Relevance")

structural_features <- rbind(structural_features_rel,structural_features_nrel)

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