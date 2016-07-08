#!/soft/R/R-3.0.0/bin/Rscript

### Switch number ###

switchNumber <- read.delim("~/SmartAS/testResults/TCGA/analysis/switchNum.txt", header=FALSE,sep=" ")
colnames(switchNumber) <- c("Cancer","Count")
patientNumber <- read.delim("~/SmartAS/testResults/TCGA/analysis/numPatients.txt", header=FALSE,sep=" ")
colnames(patientNumber) <- c("Cancer","Count")

kk <- merge(switchNumber,patientNumber,by="Cancer",suffixes = c(".switches",".patients"))
dat <- t(matrix(c(kk$Count.switches,kk$Count.patients),ncol=2))

png("~/SmartAS/testResults/TCGA/analysis/switchNumber.png", width=1500, height=1000)
barplot(dat, beside=TRUE, names.arg=switchNumber$Cancer,legend=c("Num switches","Num patients"), cex.names=2)
graphics.off()

### UTR_change ###
UTR_change <- read.delim("~/SmartAS/testResults/TCGA/analysis/UTR_change.txt", header=FALSE)
colnames(UTR_change) <- c("Cancer","Yes","No","Total")
UTR_change$Yes = as.numeric(UTR_change$Yes)/as.numeric(UTR_change$Total)
UTR_change$No = as.numeric(UTR_change$No)/as.numeric(UTR_change$Total)

dat <- t(matrix(c(UTR_change$Yes,UTR_change$No),ncol=2))

png("~/SmartAS/testResults/TCGA/analysis/UTR_change.png", width=1500, height=1000)
barplot(dat,xlab="UTR change", beside=TRUE, names.arg=UTR_change$Cancer,legend=c("Yes","No"), cex.names=2)
graphics.off()

### CDS_study ###
CDS_study <- read.delim("~/SmartAS/testResults/TCGA/analysis/CDS_study.txt", header=FALSE)
colnames(CDS_study) <- c("Cancer","Both","Only_nIso","Only_tIso","None","Total")
CDS_study$Both = as.numeric(CDS_study$Both)/as.numeric(CDS_study$Total)
CDS_study$Only_nIso = as.numeric(CDS_study$Only_nIso)/as.numeric(CDS_study$Total)
CDS_study$Only_tIso = as.numeric(CDS_study$Only_tIso)/as.numeric(CDS_study$Total)
CDS_study$None = as.numeric(CDS_study$None)/as.numeric(CDS_study$Total)

dat <- t(matrix(c(CDS_study$Both,CDS_study$Only_nIso,CDS_study$Only_tIso,CDS_study$None),ncol=4))

png("~/SmartAS/testResults/TCGA/analysis/CDS_study.png", width=1500, height=1000)
barplot(dat,xlab="CDS", beside=TRUE, names.arg=CDS_study$Cancer,legend=c("Both","Only_nIso","Only_tIso","None"), cex.names=2)
graphics.off()

### CDS_change ###
CDS_change <- read.delim("~/SmartAS/testResults/TCGA/analysis/CDS_change.txt", header=FALSE)
colnames(CDS_change) <- c("Cancer","Yes","No","Total")

CDS_change$Yes = as.numeric(CDS_change$Yes)/as.numeric(CDS_change$Total)
CDS_change$No = as.numeric(CDS_change$No)/as.numeric(CDS_change$Total)

dat <- t(matrix(c(CDS_change$Yes,CDS_change$No),ncol=2))

png("~/SmartAS/testResults/TCGA/analysis/CDS_change.png", width=1500, height=1000)
barplot(dat,xlab="CDS change", beside=TRUE, names.arg=CDS_change$Cancer,legend=c("Yes","No"), cex.names=2)
graphics.off()

### Driver_affection ###
Driver_affection <- read.delim("~/SmartAS/testResults/TCGA/analysis/Driver_affection.txt", header=FALSE)
colnames(Driver_affection) <- c("Cancer","Yes","No","Total")
Driver_affection$Yes = as.numeric(Driver_affection$Yes)/as.numeric(Driver_affection$Total)
Driver_affection$No = as.numeric(Driver_affection$No)/as.numeric(Driver_affection$Total)

dat <- t(matrix(c(Driver_affection$Yes,Driver_affection$No),ncol=2))

png("~/SmartAS/testResults/TCGA/analysis/Driver_affection.png", width=1500, height=1000)
barplot(dat,xlab="Driver", beside=TRUE, names.arg=Driver_affection$Cancer,legend=c("Yes","No"), cex.names=2)
graphics.off()

## Structural ##
### Mapped loops ###
loops <- read.delim("~/SmartAS/testResults/TCGA/analysis/loops.txt", header=FALSE)
colnames(loops) <- c("Cancer","Different","Same","Only_nIso","Only_tIso","None","Total")
loops$Different = as.numeric(loops$Different)/as.numeric(loops$Total)
loops$Same = as.numeric(loops$Same)/as.numeric(loops$Total)
loops$Only_nIso = as.numeric(loops$Only_nIso)/as.numeric(loops$Total)
loops$Only_tIso = as.numeric(loops$Only_tIso)/as.numeric(loops$Total)
loops$None = as.numeric(loops$None)/as.numeric(loops$Total)

dat <- t(matrix(c(loops$Mapped_loops,loops$Different,loops$Same,loops$Only_nIso,loops$Only_tIso,loops$None),ncol=5))

png("~/SmartAS/testResults/TCGA/analysis/loops.png", width=1500, height=1000)
barplot(dat,xlab="Mapped loops", beside=TRUE, names.arg=loops$Cancer,legend=c("Different","Same","Only_nIso","Only_tIso","None"), cex.names=2)
graphics.off()

### Functional ###
functional_change  <- read.delim("~/SmartAS/testResults/TCGA/analysis/functional_change.txt", header=FALSE)
colnames(functional_change) <- c("Cancer","Yes","No","Total")
functional_change$Yes = as.numeric(functional_change$Yes)/as.numeric(functional_change$Total)
functional_change$No = as.numeric(functional_change$No)/as.numeric(functional_change$Total)

dat <- t(matrix(c(functional_change$Yes,functional_change$No),ncol=2))

png("~/SmartAS/testResults/TCGA/analysis/functional_change.png", width=1500, height=1000)
barplot(dat,xlab="Feature changes", beside=TRUE, names.arg=functional_change$Cancer,legend=c("Yes","No"), cex.names=2)
graphics.off()

change_type  <- read.delim("~/SmartAS/testResults/TCGA/analysis/change_type.txt", header=FALSE)
colnames(change_type) <- c("Cancer","PRINTS","ProSiteProfiles","Pfam","ProSitePatterns","Total")
change_type$PRINTS = as.numeric(change_type$PRINTS)/as.numeric(change_type$Total)
change_type$ProSiteProfiles = as.numeric(change_type$ProSiteProfiles)/as.numeric(change_type$Total)
change_type$Pfam = as.numeric(change_type$Pfam)/as.numeric(change_type$Total)
change_type$ProSitePatterns = as.numeric(change_type$ProSitePatterns)/as.numeric(change_type$Total)

dat <- t(matrix(c(change_type$PRINTS,change_type$ProSiteProfiles,change_type$ProSitePatterns,change_type$Pfam),ncol=4))

png("~/SmartAS/testResults/TCGA/analysis/functional_change_type.png", width=1500, height=1000)
barplot(dat,xlab="Feature change DB", beside=TRUE, names.arg=change_type$Cancer,legend=c("PRINTS","ProSiteProfiles","ProSitePatterns","Pfam"), cex.names=2)
graphics.off()

specificFeatures <- read.delim("~/SmartAS/testResults/TCGA/analysis/specificFeatures.txt", header=FALSE)
colnames(specificFeatures) <- c("Cancer","Feature","Count")
specificFeatures$Count = as.numeric(specificFeatures$Count)

specificFeatures <- aggregate(Count~Feature,specificFeatures,FUN=sum)
specificFeatures <- specificFeatures[with(specificFeatures, order(-Count)), ]
specificFeatures <- specificFeatures[1:10,]

png("~/SmartAS/testResults/TCGA/analysis/motifs.png", width=1500, height=1000)
par(las=2)
barplot(specificFeatures$Count,xlab="Feature change DB",horiz=TRUE,names.arg=specificFeatures$Feature, cex.names=2)
graphics.off()

## Disorderd
disordered_change <- read.delim("~/SmartAS/testResults/TCGA/analysis/disordered_change.txt", header=FALSE)
colnames(disordered_change) <- c("Cancer","Yes","No","Total")
disordered_change$Yes = as.numeric(disordered_change$Yes)/as.numeric(disordered_change$Total)
disordered_change$No = as.numeric(disordered_change$No)/as.numeric(disordered_change$Total)

dat <- t(matrix(c(disordered_change$Yes,disordered_change$No),ncol=2))

png("~/SmartAS/testResults/TCGA/analysis/disordered_change.png", width=1500, height=1000)
barplot(dat,xlab="Change in disordered region",beside=TRUE,names.arg=disordered_change$Cancer, cex.names=2)
graphics.off()

specificFeatures <- read.delim("~/SmartAS/testResults/TCGA/analysis/disordered_meanLength.txt", header=FALSE)
colnames(specificFeatures) <- c("Cancer","meanLength")
meanMeanLength <- mean(as.numeric(specificFeatures$meanLength))
cat("mean length",meanMeanLength,"\n")

### Accuracy ###
accKmeans <- read.delim("~/SmartAS/testResults/TCGA/analysis/kmeans.txt", header=FALSE)
colnames(accKmeans) <- c("cancer","precision","sensitivity","precision_pval","sensitivity_pval")
meanKmeansPrecision <- mean(accKmeans$precision)
meanKmeansSensitiv <- mean(accKmeans$sensitivity)
meanKmeansPP <- mean(accKmeans$precision_pval)
meanKmeansPS <- mean(accKmeans$sensitivity_pval)

accHclust <- read.delim("~/SmartAS/testResults/TCGA/analysis/hclust.txt", header=FALSE)
colnames(accHclust) <- c("cancer","precision","sensitivity","precision_pval","sensitivity_pval")
meanHclustPrecision <- mean(accHclust$precision)
meanHclustSensitiv <- mean(accHclust$sensitivity)
accHclustPP <- mean(accHclust$precision_pval)
accHclustPS <- mean(accHclust$sensitivity_pval)

cat("Precision: ",meanKmeansPrecision,meanHclustPrecision,"\n")
cat("Sensitivity: ",meanKmeansSensitiv,meanHclustSensitiv,"\n")
cat("Precision-pvalue: ",meanKmeansPP,accHclustPP,"\n")
cat("Sensitivity-pvalue: ",meanKmeansPS,accHclustPS,"\n")

