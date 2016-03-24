# ################ CENTRALITY ANALYSIS ################
# #centrality_all <- data.frame()
# centrality_rel <- data.frame()
# centrality_nrel <- data.frame()
# 
# for (cancer in cancerTypes){
#   basePath <- ""
#   #cancer_all <- read.delim(paste0(basePath,cancer,".tIso_length.tsv"), header=FALSE)
#   #tiso_length_all <- rbind(tiso_length_all,cbind(cancer_all,cancer))
#   
#   cancer_rel <- read.delim(paste0(basePath,cancer,".protein_centrality_functional.tsv"), header=FALSE)
#   cancer_rel <- cbind(cancer_rel,"Functional",cancer)
#   colnames(cancer_rel) <- c("Cancer","Degree","Relevance")
#   centrality_rel <- rbind(centrality_rel,cancer_rel)
#   
#   cancer_nrel <- read.delim(paste0(basePath,cancer,".protein_centrality_nonfunctional.tsv"), header=FALSE)
#   cancer_nrel <- cbind(cancer_nrel,"NonFunctional",cancer)
#   colnames(cancer_nrel) <- c("Cancer","Degree","Relevance")
#   centrality_nrel <- rbind(centrality_nrel,cancer_nrel)
# }
# 
# centrality <- rbind(centrality_rel,centrality_nrel)
# 
# png("centrality.png", width=1000, height=800)
# par(mar=c(9,3,5,2))
# boxplot(Degree~interaction(Relevance,Cancer),data=centrality,outline=F,col=c("gray40","gray90"),las=2)
# graphics.off()
# 
# ################ FEATURES AFFECTED ################
# Feats_affected_rel <- read.delim("structural_summary_functional.tsv", header=FALSE)
# Feats_affected_rel <- cbind(Feats_affected_rel,"Functional")
# colnames(Feats_affected_rel) <- c("Cancer","Switch","Pfam","PRINTS","ProSitePatterns","IUPREDLong","IUPREDShort","I3D","Driver","Relevance")
# 
# Feats_affected_nrel <- read.delim("structural_summary_nonfunctional.tsv", header=FALSE)
# Feats_affected_nrel <- cbind(Feats_affected_nrel,"NonFunctional")
# colnames(Feats_affected_nrel) <- c("Cancer","Switch","Pfam","PRINTS","ProSitePatterns","IUPREDLong","IUPREDShort","I3D","Driver","Relevance")
# 
# Feats_affected <- rbind(Feats_affected_rel,Feats_affected_nrel)
# Feats_affected$Driver <- ifelse(Feats_affected$Driver=="True","Driver","NonDriver")
# 
# png("ProSitePatterns.png", width=1000, height=800)
# par(mar=c(9,3,5,2))
# boxplot(ProSitePatterns~interaction(Driver,Cancer),data=Feats_affected[Feats_affected$Relevance=="Functional" & Feats_affected$ProSitePatterns!=0,],col=c("gray40","gray90"),las=3,outline=F)
# graphics.off()
# 
# png("Pfam.png", width=1000, height=800)
# par(mar=c(9,3,5,2))
# boxplot(Pfam~interaction(Driver,Cancer),data=Feats_affected[Feats_affected$Relevance=="Functional" & Feats_affected$Pfam!=0,],col=c("gray40","gray90"),las=3,outline=F)
# graphics.off()
# 
# png("PRINTS.png", width=1000, height=800)
# par(mar=c(9,3,5,2))
# boxplot(PRINTS~interaction(Driver,Cancer),data=Feats_affected[Feats_affected$Relevance=="Functional" & Feats_affected$PRINTS!=0,],col=c("gray40","gray90"),las=3,outline=F)
# graphics.off()
# 
# png("I3D.png", width=1000, height=800)
# par(mar=c(9,3,5,2))
# boxplot(I3D~interaction(Driver,Cancer),data=Feats_affected[Feats_affected$Relevance=="Functional" & Feats_affected$I3D!=0,],col=c("gray40","gray90"),las=3,outline=F)
# graphics.off()
# 
# png("IUPREDShort.png", width=1000, height=800)
# par(mar=c(9,3,5,2))
# boxplot(IUPREDShort~interaction(Driver,Cancer),data=Feats_affected[Feats_affected$Relevance=="Functional" & Feats_affected$IUPREDShort!=0,],col=c("gray40","gray90"),las=3,outline=F)
# graphics.off()
# 
# png("IUPREDLong.png", width=1000, height=800)
# par(mar=c(9,3,5,2))
# boxplot(IUPREDLong~interaction(Driver,Cancer),data=Feats_affected[Feats_affected$Relevance=="Functional" & Feats_affected$IUPREDLong!=0,],col=c("gray40","gray90"),las=3,outline=F)
# graphics.off()
# 
# ################ STRUCTURAL FEATURES ################
# #exons_all <- read.delim("exons.tsv", header=FALSE)
# #exons_all <- cbind(exons_all,"All")
# #colnames(exons_all) <- c("Cancer","switch","length","type","keepORF","position","Relevance")
# 
# structural_features <- read.delim("structural_features.onlyModels.tsv", header=FALSE)
# colnames(structural_features) <- c("Cancer","Gene","Symbol","nIso","tIso","Random","Analysis","WhatsHappening","Feature","Driver","ASDriver","DriverType")
# nrStructural_features <- structural_features[structural_features$Random == "NonRandom",]
# 
# PfamSomethingHappening = structural_features[structural_features$Analysis=='Pfam' & structural_features$WhatsHappening!="Nothing",]
# PfamNothingHappening = structural_features[structural_features$Analysis=='Pfam' & structural_features$WhatsHappening!="Nothing",]
# a = sum(PfamSomethingHappening$Random == "NonRandom")
# b = sum(PfamSomethingHappening$Random == "Random")
# c = sum(PfamNothingHappening$Random == "NonRandom")
# d = sum(PfamNothingHappening$Random == "Random")
# 
# fisher.test(matrix(c(a,b,c,d),nrow=2,ncol=2))
# 
# PfamDriver <- sort(table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$Driver==1 & structural_features$WhatsHappening!="Nothing"]),decreasing=TRUE)
# PfamDriver <- PfamDriver[PfamDriver>0]
# PfamDriver <- 100*PfamDriver/length(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$Driver==1 & structural_features$WhatsHappening!="Nothing"])
# 
# PfamNonDriver <-sort(table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$Driver==0 & structural_features$WhatsHappening!="Nothing"]),decreasing=TRUE)
# PfamNonDriver <- PfamNonDriver[PfamNonDriver>0]
# PfamNonDriver <- 100*PfamNonDriver/length(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$Driver==0 & structural_features$WhatsHappening!="Nothing"])
# 
# PfamDriverGained <- table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$Driver==1 & structural_features$WhatsHappening=="Gained_in_tumor"])
# PfamDriverLost <- table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$Driver==1  & structural_features$WhatsHappening=="Lost_in_tumor"])
# PfamDriverDifferences <- PfamDriverGained - PfamDriverLost 
# PfamDriverDifferences <- sort(PfamDriverDifferences,decreasing=TRUE)
# PfamDriverDifferences <- PfamDriverDifferences[PfamDriverDifferences!=0]
# write.table(PfamDriverDifferences,'PfamDriverDifferences.txt',quote=F,col.names=F)
# 
# PfamOncogeneGained <- table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$DriverType=="oncogene" & structural_features$WhatsHappening=="Gained_in_tumor"])
# PfamOncogeneLost <- table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$DriverType=="oncogene"  & structural_features$WhatsHappening=="Lost_in_tumor"])
# PfamOncogeneDifferences <- PfamOncogeneGained - PfamOncogeneLost
# PfamOncogeneDifferences <- sort(PfamOncogeneDifferences,decreasing=TRUE)
# PfamOncogeneDifferences <- PfamOncogeneDifferences[PfamOncogeneDifferences!=0]
# write.table(PfamOncogeneDifferences,'PfamOncogeneDifferences.txt',quote=F,col.names=F)
# 
# PfamSuppressorGained <- table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$DriverType=="suppressor" & structural_features$WhatsHappening=="Gained_in_tumor"])
# PfamSuppressorLost <- table(structural_features$Feature[structural_features$Analysis=='Pfam' & structural_features$DriverType=="suppressor"  & structural_features$WhatsHappening=="Lost_in_tumor"])
# PfamSuppressorDifferences <- PfamSuppressorGained - PfamSuppressorLost
# PfamSuppressorDifferences <- sort(PfamSuppressorDifferences,decreasing=TRUE)
# PfamSuppressorDifferences <- PfamSuppressorDifferences[PfamSuppressorDifferences!=0]
# write.table(PfamSuppressorDifferences,'PfamSuppressorDifferences.txt',quote=F,col.names=F)
# 
# 
# ## ---- to review ---- ##
# 
# ProSiteDriverGained <- sort(table(structural_features$Feature[structural_features$Analysis=='ProSitePatterns' & structural_features$Driver==1 & structural_features$WhatsHappening=="Gained_in_tumor"]),decreasing=TRUE)
# ProSiteDriverLost <- sort(table(structural_features$Feature[structural_features$Analysis=='ProSitePatterns' & structural_features$Driver==1  & structural_features$WhatsHappening=="Lost_in_tumor"]),decreasing=TRUE)
# ProSiteDifference <- ProSiteDriverGained - ProSiteDriverLost
# ProSiteDifference <- sort(ProSiteDifference,decreasing=TRUE)
# ProSiteDifference <- ProSiteDifference[ProSiteDifference!=0]
# write.table(ProSiteDifference,'ProSiteDifference.txt',quote=F,col.names=F)
# 
# 
# ProSitePatternsDriver <-sort(table(structural_features$Feature[structural_features$Analysis=='ProSitePatterns' & structural_features$Driver=="True"]),decreasing=TRUE)
# ProSitePatternsDriver <- ProSitePatternsDriver[ProSitePatternsDriver>0]
# ProSitePatternsDriver <- 100*ProSitePatternsDriver/length(structural_features$Feature[structural_features$Analysis=='ProSitePatterns' & structural_features$Driver=="True"])
# 
# ProSitePatternsNonDriver <-sort(table(structural_features$Feature[structural_features$Analysis=='ProSitePatterns' & structural_features$Driver=="False"]),decreasing=TRUE)
# ProSitePatternsNonDriver <- ProSitePatternsNonDriver[ProSitePatternsNonDriver>0]
# ProSitePatternsNonDriver <- 100*ProSitePatternsNonDriver/length(structural_features$Feature[structural_features$Analysis=='ProSitePatterns' & structural_features$Driver=="False"])
# 
# 
# 
# ProSitePatternsDriverLost <- ProSitePatternsDriverLost[ProSitePatternsDriverLost>0]
# ProSitePatternsDriverLost <- 100*ProSitePatternsDriverLost/length(structural_features$Feature[structural_features$Analysis=='ProSitePatterns' & structural_features$Driver=="True"  & structural_features$WhatsHappening=="Lost_in_tumor"])
# 
# write.table(ProSitePatternsDriver,'ProSitePatternsDriver.txt',quote=F,col.names=F)
# write.table(ProSitePatternsNonDriver,'ProSitePatternsNonDriver.txt',quote=F,col.names=F)
# write.table(ProSitePatternsDriverLost,'ProSitePatternsDriverLost.txt',quote=F,col.names=F)
# write.table(ProSitePatternsDriverGained,'ProSitePatternsDriverGained.txt',quote=F,col.names=F)
# 
# #featureCounts <- by(structural_features, structural_features$Analysis, table(x[,5))