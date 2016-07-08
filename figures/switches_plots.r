#!/usr/bin/env Rscript

source("~/smartas/pipeline/scripts/variablesAndFunctions.r")

setwd(workingDir)
setwd("switches")

# 1 - CANDIDATES STUDY ---------------------
# 1.1 - read all candidate lists ====
candidateList <- list()
for (cancer in cancerTypes){
  candidateList[[cancer]] <- read.delim(paste0(cancer,".candidateList.tsv"))
  candidateList[[cancer]] <- candidateList[[cancer]][candidateList[[cancer]]$IsModel==1 & candidateList[[cancer]]$NotNoise==1,]
  candidateList[[cancer]]$NumPatients <- unlist(lapply(strsplit(as.character(candidateList[[cancer]]$Patients_affected),",",fixed=T),length))
  candidateList[[cancer]]$Percentage <- candidateList[[cancer]]$NumPatients/nPatients[[cancer]]
  candidateList[[cancer]]$Tumor <- cancer
}
candidatesDf <- do.call("rbind", candidateList)
candidatesDf <- merge(candidatesDf,nPatientsDf,by.x="Tumor",by.y="Cancer")

# join driver type
driverTypesFile <- paste0(workingDir,"/data/Databases/cancer_networks_SuppTables_v7_S7.csv")
driverTypes <- read.delim(driverTypesFile,header=F)
colnames(driverTypes) <- c("Symbol","DriverType")

candidatesDf <- merge(candidatesDf,driverTypes,all.x=TRUE)
candidatesDf$DriverType <- as.character(candidatesDf$DriverType)
candidatesDf$DriverType[is.na(candidatesDf$DriverType)] <- "No"

x <- ddply(candidatesDf,.(Tumor),summarise,
           median=median(NumPatients),
           mad=mad(NumPatients,na.rm=T))

y <- merge(x,candidatesDf)
y$NumPatients.ZScore <- (y$NumPatients - y$median)/(1.486*y$mad)

candidatesDf <- merge(candidatesDf,y[,c("Tumor","GeneId","Normal_transcript","Tumor_transcript","NumPatients.ZScore")])

write.table(candidatesDf,'tables/candidateList_splitByTumor_models_notNoise.txt',quote=F,row.names=F, sep="\t")
rm(driverTypes,driverTypesFile)

# 1.2 - join all candidate lists ====

# calculate aggregated table (sum patients, etc.) and calculate unbalance as the minimun p of enrichment in a fisher test
candidatesDf_agg <- ddply(candidatesDf,.(GeneId,Symbol,Normal_transcript,Tumor_transcript,
                                         Normal_protein,Tumor_protein,Annotation,
                                         DriverAnnotation,IsFunctional,Driver,DriverType,
                                         Druggable,CDS_Normal,CDS_Tumor,CDS_change,UTR_change),
                          summarise, CancerAffected=paste(Tumor,collapse = ","),
                          Patients=paste(Patients_affected,collapse = ","), 
                          PatientNumber=sum(NumPatients), 
                          Percentage=sum(NumPatients)/nPatients[["total"]],
                          p.unbalance=min(apply(
                            cbind(NumPatients,sum(NumPatients)-NumPatients,
                                  TotalPatients-NumPatients,
                                  sum(TotalPatients)-TotalPatients-(sum(NumPatients)-NumPatients)),
                            1,function(x) { 
                              if (x[2]!=0){
                                f <- fisher.test(matrix(x[1:4],ncol=2),alternative="greater"); 
                                f$p.value
                              } else { 0 } } )),
                          Entropy=-sum((NumPatients/TotalPatients)/
                                         sum(NumPatients/TotalPatients)*
                                         log2(
                                           (NumPatients/TotalPatients)/
                                             sum(NumPatients/TotalPatients))
                          )/log2(length(NumPatients))
)

candidatesDf_agg <- candidatesDf_agg[order(-candidatesDf_agg$Percentage),]
candidatesDf_agg <- candidatesDf_agg[,c("GeneId","Symbol","Normal_transcript",
                                        "Tumor_transcript","Normal_protein",
                                        "Tumor_protein","Annotation","DriverAnnotation",
                                        "IsFunctional","Driver","DriverType","Druggable",
                                        "CDS_Normal","CDS_Tumor","CDS_change","UTR_change",
                                        "CancerAffected","PatientNumber","Percentage",
                                        "p.unbalance","Patients","Entropy")]
write.table(candidatesDf_agg,'tables/candidateList_allCancers_models_notNoise.txt',quote=F,row.names=F, sep="\t")

# 1.3 - calculate most frequent genes altered ====

df <- candidatesDf
df$Symbol <- as.character(df$Symbol)
df$Symbol[ df$DriverAnnotation=="Driver" ] <- paste0("*",df$Symbol[ df$DriverAnnotation=="Driver" ],"*")

#most frequent genes
for (thing in c("# patients","% patients")){
  for (effect in c("allGenes","functionalGenes")){
    for (annotation in c("allGenes","driver","d1")){
      ngenes <- 50
      
      if (annotation=="allGenes"){
        genesSelection <- matrix(T,nrow(df),1)
      } else {
        genesSelection <- df$DriverAnnotation==annotation
      }
      
      if (effect=="functionalGenes") { 
        functionalSelection <- as.logical(df$IsFunctional)
      } else {
        functionalSelection <- matrix(T,nrow(df),1)
      }
      
      if (thing=="# patients"){
        accCandidateList <- ddply(df[genesSelection & functionalSelection,],.(GeneId,Symbol), summarise, cumsum=sum(NumPatients))
      } else if (thing=="% patients"){
        accCandidateList <- ddply(df[genesSelection & functionalSelection,],.(GeneId,Symbol), summarise, cumsum=sum(Percentage))
      }
      
      accCandidateList <- accCandidateList[order(-accCandidateList$cumsum),]
      
      candsMatrix <- list()
      for (gene in accCandidateList$GeneId[1:ngenes]){
        symbol = accCandidateList$Symbol[accCandidateList$GeneId==gene]
        candsMatrix[[symbol]] <- list()
        candsMatrix[[symbol]][["symbol"]] <- as.character(symbol)
        for (cancer in cancerTypes){
          if (gene %in% candidateList[[cancer]]$GeneId){
            if (thing=="# patients"){
              candsMatrix[[symbol]][[cancer]] <- as.numeric(candidateList[[cancer]]$NumPatients[candidateList[[cancer]]$GeneId==gene])
            } else if (thing=="% patients"){
              candsMatrix[[symbol]][[cancer]] <- as.numeric(candidateList[[cancer]]$Percentage[candidateList[[cancer]]$GeneId==gene])
            }
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
      if (thing=="% patients"){
        candsMatrix3 <- candsMatrix3*100/11
      }
      candsMatrix3$Gene <- rownames(candsMatrix3)
      
      candsMatrix3.m <- melt(candsMatrix3)
      colnames(candsMatrix3.m) <- c("Gene","Cancer","patients")
      
      candsMatrix3.m$Gene <- factor(candsMatrix3.m$Gene,levels=as.factor(unique(candsMatrix3.m$Gene)))
      
      p <- ggplot(candsMatrix3.m) + 
        geom_bar(stat="identity",aes(x=Gene,y=patients,fill=Cancer)) + 
        scale_fill_manual(values=colorPalette) + 
        labs(x="Genes",y=thing) + 
        theme(text = element_text(size=20),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour="black"))
      
      if (thing=="# patients"){
        ggsave(paste0("figures/",paste("genes","numPatients",annotation,effect,paste0("top",ngenes,"genes"),"allCancers.png",sep="_")),p,width = 12)
      } else if (thing=="% patients"){
        ggsave(paste0("figures/",paste("genes","percPatients",annotation,effect,paste0("top",ngenes,"genes"),"allCancers.png",sep="_")),p,width = 12)
      }
    }
  }
}

#most frequent switches
for (thing in c("# patients","% patients")){
  for (effect in c("allGenes","functionalGenes")){
    for (annotation in c("allGenes","driver","d1")){
      ngenes <- 50
      
      if (annotation=="allGenes"){
        genesSelection <- matrix(T,nrow(df),1)
      } else {
        genesSelection <- df$DriverAnnotation==annotation
      }
      
      if (effect=="functionalGenes") { 
        functionalSelection <- as.logical(df$IsFunctional)
      } else {
        functionalSelection <- matrix(T,nrow(df),1)
      }
      
      if (thing=="# patients"){
        accCandidateList <- ddply(df[genesSelection & functionalSelection,],.(GeneId,Symbol,Normal_transcript,Tumor_transcript), summarise, cumsum=sum(NumPatients))
      } else if (thing=="% patients"){
        accCandidateList <- ddply(df[genesSelection & functionalSelection,],.(GeneId,Symbol,Normal_transcript,Tumor_transcript), summarise, cumsum=sum(Percentage))
      }
      
      accCandidateList <- accCandidateList[order(-accCandidateList$cumsum),]
      accCandidateList$Switch <- apply(subset(accCandidateList,select = c(Symbol,Normal_transcript,Tumor_transcript)),1, function(x){paste(x[1],x[2],x[3],sep=" ")})
      
      candsMatrix <- list()
      for (s in accCandidateList$Switch[1:ngenes]){
        candsMatrix[[s]] <- list()
        candsMatrix[[s]][["switch"]] <- s
        sp <- unlist(strsplit(s," "))
        nTx <- sp[2]
        tTx <- sp[3]
        for (cancer in cancerTypes){
          if (nTx %in% candidateList[[cancer]]$Normal_transcript & tTx %in% candidateList[[cancer]]$Tumor_transcript){
            if (thing=="# patients"){
              candsMatrix[[s]][[cancer]] <- as.numeric(candidateList[[cancer]]$NumPatients[candidateList[[cancer]]$Normal_transcript==nTx & candidateList[[cancer]]$Tumor_transcript==tTx])
            } else if (thing=="% patients"){
              candsMatrix[[s]][[cancer]] <- as.numeric(candidateList[[cancer]]$Percentage[candidateList[[cancer]]$Normal_transcript==nTx & candidateList[[cancer]]$Tumor_transcript==tTx])
            }
          } else{
            candsMatrix[[s]][[cancer]] <- 0
          }
        }
        candsMatrix[[s]] <- do.call("cbind", candsMatrix[[s]])
      }
      candsMatrix2 <- do.call("rbind", candsMatrix)
      rownames(candsMatrix2) <- candsMatrix2[,c("switch")]
      
      candsMatrix3 <- as.data.frame(matrix(as.numeric(as.character(candsMatrix2[,!(colnames(candsMatrix2) %in% c("switch"))])),nrow=ngenes,ncol=11))
      rownames(candsMatrix3) <- rownames(candsMatrix2)
      colnames(candsMatrix3) <- cancerTypes
      candsMatrix3 <- candsMatrix3[order(-rowSums(candsMatrix3)),]
      if (thing=="% patients"){
        candsMatrix3 <- candsMatrix3*100/11
      }
      candsMatrix3$Gene <- rownames(candsMatrix3)
      
      candsMatrix3.m <- melt(candsMatrix3)
      colnames(candsMatrix3.m) <- c("Gene","Cancer","patients")
      candsMatrix3.m$Gene <- factor(candsMatrix3.m$Gene,levels=as.factor(rownames(candsMatrix3)))
      
      p <- ggplot(candsMatrix3.m) + 
        geom_bar(stat="identity",aes(x=Gene,y=patients,fill=Cancer)) + 
        scale_fill_manual(values=colorPalette) + 
        labs(x="Switches",y=thing) + 
        theme(text = element_text(size=20),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour="black"))
      if (thing=="# patients"){
        ggsave(paste0("figures/",paste("switches","numPatients",annotation,effect,paste0("top",ngenes,"switches"),"allCancers.png",sep="_")),p,width = 12)
      } else if (thing=="% patients"){
        ggsave(paste0("figures/",paste("switches","percPatients",annotation,effect,paste0("top",ngenes,"switches"),"allCancers.png",sep="_")),p,width = 12)
      }
      
    }
  }
}

rm(effect,annotation,cancer,functionalSelection,gene,
   genesSelection,ngenes,symbol,p,candsMatrix,candsMatrix2,
   candsMatrix3,candsMatrix3.m,accCandidateList)

# 1.4 - Study exon length ====
# exons.tsv study
exons <- read.delim("exons.tsv")

# cds relative size
cdsRelativeSize <- list()
cdsPosition <- list()
exonLength <- list()
orfChange <- data.frame(cancer=c(),p=c(),oddsRatio=c())
exonOrigin <- data.frame(cancer=c(),p=c(),oddsRatio=c())
for (knsur in cancerTypes){
  cancer.exons = subset(exons,Cancer==knsur & (Origin=="CDS" | Origin=="CDS-UTR") )
  
  # cds relative size
  cdsRelativeSize[[knsur]] <- ggplot() + 
    stat_ecdf(data=subset(cancer.exons,Random=="Random"), aes(x=CDSRelativeSize,colour="green"),show_guide = FALSE) + 
    stat_ecdf(data=subset(cancer.exons,Random=="NonRandom"), aes(x=CDSRelativeSize,colour="red"),show_guide = FALSE) +
    ggtitle(knsur) + 
    theme_bw() + 
    ylab("")
  
  # cds position
  cdsPosition[[knsur]] <- ggplot() + ggtitle(knsur) + theme_bw() + 
    stat_ecdf(data=subset(cancer.exons,Random=="NonRandom"), aes(x=Position,colour="red"),show_guide = FALSE) + 
    stat_ecdf(data=subset(cancer.exons,Random=="Random"), aes(x=Position,colour="green"),show_guide = FALSE)
  
  # length
  exonLength[[knsur]] <- ggplot() + ggtitle(knsur) + theme_bw() + 
    stat_ecdf(data=subset(cancer.exons,Random=="NonRandom"), aes(x=Length,colour="red"),show_guide = FALSE) + 
    stat_ecdf(data=subset(cancer.exons,Random=="Random"), aes(x=Length,colour="green"),show_guide = FALSE)
  
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

# to put a legend common to all plots
# https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
png("figures/exon_cds_relative_size.png", width=800, height=1000)
grid.arrange(cdsRelativeSize[["brca"]],cdsRelativeSize[["coad"]],cdsRelativeSize[["hnsc"]],cdsRelativeSize[["kich"]],cdsRelativeSize[["kirc"]],cdsRelativeSize[["kirp"]],cdsRelativeSize[["lihc"]],cdsRelativeSize[["luad"]],cdsRelativeSize[["lusc"]],cdsRelativeSize[["prad"]],cdsRelativeSize[["thca"]])
graphics.off()

png("figures/exon_cds_position.png", width=800, height=1000)
grid.arrange(cdsPosition[["brca"]],cdsPosition[["coad"]],cdsPosition[["hnsc"]],cdsPosition[["kich"]],cdsPosition[["kirc"]],cdsPosition[["kirp"]],cdsPosition[["lihc"]],cdsPosition[["luad"]],cdsPosition[["lusc"]],cdsPosition[["prad"]],cdsPosition[["thca"]])
graphics.off()

png("figures/exonlength.png", width=800, height=1000)
grid.arrange(exonLength[["brca"]],exonLength[["coad"]],exonLength[["hnsc"]],exonLength[["kich"]],exonLength[["kirc"]],exonLength[["kirp"]],exonLength[["lihc"]],exonLength[["luad"]],exonLength[["lusc"]],exonLength[["prad"]],exonLength[["thca"]])
graphics.off()

write.table(orfChange,'tables/exonOrfChange_testVsRandom.txt',quote=F,row.names=F, sep="\t")
write.table(exonOrigin,'tables/exonOrigin_NorT_testVsRandom.txt',quote=F,row.names=F, sep="\t")

rm(cdsRelativeSize,cdsPosition,exonLength,knsur,cTable,exons,
   this.Data,exonOrigin,f,switches,randomSwitches,orfChange,cancer.exons)

# exons_new.tsv study
exonsNew <- read.delim("exons_new.tsv")
nrExonsNew <- exonsNew[exonsNew$Random=="NonRandom" ,]

sum(nrExonsNew$nExon%%3==nrExonsNew$tExon%%3)/nrow(nrExonsNew)
sum(nrExonsNew$nExon[nrExonsNew$Tag=="BEGINNING"]%%3==nrExonsNew$tExon[nrExonsNew$Tag=="BEGINNING"]%%3)/nrow(nrExonsNew[nrExonsNew$Tag=="BEGINNING",])
sum(nrExonsNew$nExon[nrExonsNew$Tag=="ENDING"]%%3==nrExonsNew$tExon[nrExonsNew$Tag=="ENDING"]%%3)/nrow(nrExonsNew[nrExonsNew$Tag=="ENDING",])
sum(nrExonsNew$nExon[nrExonsNew$Tag=="MIDDLE"]%%3==nrExonsNew$tExon[nrExonsNew$Tag=="MIDDLE"]%%3)/nrow(nrExonsNew[nrExonsNew$Tag=="MIDDLE",])
sum(nrExonsNew$nExon[nrExonsNew$Tag=="COMMON"]%%3==nrExonsNew$tExon[nrExonsNew$Tag=="COMMON"]%%3)/nrow(nrExonsNew[nrExonsNew$Tag=="COMMON",])

rm(exonsNew,nrExonsNew)

# 1.5 - Study exons per switch ====
# 
# 
# exonsPerSwitch_rel <- read.delim("exonsPerSwitch.tsv", header=FALSE)
# exonsPerSwitch_rel <- cbind(exonsPerSwitch_rel,"Functional")
# colnames(exonsPerSwitch_rel) <- c("Counts","Cancer","Switch","Relevance")
# 
# exonsPerSwitch_nrel <- read.delim("exonsPerSwitch_nonfunctional.tsv", header=FALSE)
# exonsPerSwitch_nrel <- cbind(exonsPerSwitch_nrel,"NonFunctional")
# colnames(exonsPerSwitch_nrel) <- c("Counts","Cancer","Switch","Relevance")
# 
# exonsPerSwitch <- rbind(exonsPerSwitch_rel,exonsPerSwitch_nrel)
# 
# png("figures/exons_per_switch.png", width=1000, height=800)
# par(mar=c(9,3,5,2))
# boxplot(Counts~interaction(Relevance,Cancer),data=exonsPerSwitch,outline=F,col=c("gray40","gray90"),las=2)
# graphics.off()
# 
# barplot(kk1,names.arg=unique(kk$Cancer),beside=T)
# barplot(scale(kk1, center=FALSE, scale=colSums(kk1)),names.arg=unique(kk$Cancer),beside=T)

# 1.7 - Generate plot about relevance of the switches ====

candidatesDf_copy <- candidatesDf
drivers <- candidatesDf_copy$DriverAnnotation=="Driver"
candidatesDf_copy$Driver <- as.character(candidatesDf_copy$Driver)
candidatesDf_copy$Driver[drivers] <- "Unlabeled driver"
candidatesDf_copy$Driver[!drivers] <- "NonDriver"
candidatesDf_copy$Driver[candidatesDf_copy$DriverType == "oncogene"] <- "Oncogene"
candidatesDf_copy$Driver[candidatesDf_copy$DriverType == "suppressor"] <- "Suppressor"
candidatesDf_copy$Driver <- factor(candidatesDf_copy$Driver,c("NonDriver","Unlabeled driver","Oncogene","Suppressor"))

candsStats <- ddply(candidatesDf_copy,.(Tumor,IsFunctional,Driver),summarise,Count=length(GeneId))
colnames(candsStats) <- c("Cancer","Functional","Driver","Count")
func <- as.logical(candsStats$Functional)
candsStats$Functional[func] <- "Functional"
candsStats$Functional[!func] <- "NonFunctional"
candsStats$Functional <- factor(candsStats$Functional,c("Functional","NonFunctional"))

p <- ggplot(candsStats) + 
  geom_bar(stat="identity",aes(x=Cancer,y=Count,fill=Functional,alpha=Driver)) + 
  theme_minimal() + ylab("Number of switches") + scale_alpha_manual(values=c(1,0.2,0.5,0.8)) + 
  theme(text = element_text(size=20),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour="black")) + 
  scale_fill_manual(values=c("#d95f02","#7570b3"))

ggsave(paste0("figures/switchNumber_funcional_drivers.png"),p, width = 8, height = 7)

rm(p,candsStats,func,candidatesDf_copy,drivers)

# 1.8 - Plot number of transcripts vs number of patients affected by a switch ====

txsAndGenes <- read.delim(paste0(workingDir,"/data/ucsc/genesAndTranscripts_splitByGene.txt"))

txCount <- ddply(txsAndGenes,.(GeneId),summarize,nTxs=length(tx))
txCount_withInfo <- merge(candidatesDf_agg,txCount,all.x=T)
txCount_withInfo <- ddply(txCount_withInfo,.(GeneId,nTxs),summarize,PatientNumber=sum(PatientNumber))

txCount_withInfo$nTxs_factor <- cut(txCount_withInfo$nTxs, c(1:10,11))
txCount_withInfo$nTxs_factor <- factor(txCount_withInfo$nTxs_factor)

p <- ggplot(txCount_withInfo,aes(nTxs_factor,PatientNumber)) + 
  geom_boxplot() + 
  smartas_theme() + 
  xlab("Number of transcripts") + 
  ylab("Patients with switch")
ggsave(paste0("figures/gene_numberOfPatientsWithSwitch_vs_numberOfTranscripts.png"),p)

rm(p,txsAndGenes,txCount,txCount_withInfo)

# 1.9 - Check correlation with immune infiltration ====

nonTumorCells <- list()

for (tumor in cancerTypes){
  stromal.correlation <- read.table(paste0("/projects_rg/TCGA/pipeline/run11/",tumor,"_gene_gsea_full.txt"), check.names=FALSE, header=TRUE)
  stromal.correlation.t <- stromal.correlation[grep("T$",stromal.correlation$sample),]
  
  switches.tumor <- candidatesDf[candidatesDf$Tumor==tumor,]
  
  p <- apply(switches.tumor,1,function(z){
    y <- z[19]
    x <- unlist(strsplit(y,","))
    
    p <- wilcox.test(stromal.correlation.t$immune_set[stromal.correlation.t$sample %in% x],
                     stromal.correlation.t$immune_set[! stromal.correlation.t$sample %in% x])
    q <- wilcox.test(stromal.correlation.t$stromal_set[stromal.correlation.t$sample %in% x],
                     stromal.correlation.t$stromal_set[! stromal.correlation.t$sample %in% x])
           
    c(p$p.value,q$p.value)
    
  })
  
  df <- as.data.frame(t(p))
  colnames(df) <- c("immune_set","stromal_set")
  df <- cbind(switches.tumor[,c("Tumor","GeneId","Normal_transcript","Tumor_transcript")],df)
  
  nonTumorCells[[tumor]] <- df
}

nonTumorCells.all <- do.call("rbind",nonTumorCells)

write.table(nonTumorCells.all,'tables/immune_strome_contamination.txt',quote=F,row.names=F, sep="\t")
rm(nonTumorCells,p)