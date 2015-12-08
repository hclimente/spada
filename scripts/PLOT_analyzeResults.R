#!/soft/R/R-3.2.1/bin/Rscript

tryCatch(source("~/SmartAS/Pipeline/scripts/variablesAndFunctions.r"),error=function(e){})
tryCatch(source("~/Pipeline/scripts/variablesAndFunctions.r"),error=function(e){})

setwd(workingDir)

# 1 - CANDIDATES STUDY ---------------------
# 1.1 - read all candidate lists ====
setwd("switches")

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
driverTypesFile <- paste0(workingDir,"/Data/Databases/cancer_networks_SuppTables_v7_S7.csv")
driverTypes <- read.delim(driverTypesFile,header=F)
colnames(driverTypes) <- c("Symbol","DriverType")

candidatesDf <- merge(candidatesDf,driverTypes,all.x=TRUE)
candidatesDf$DriverType <- as.character(candidatesDf$DriverType)
candidatesDf$DriverType[is.na(candidatesDf$DriverType)] <- "No"

write.table(candidatesDf,'tables/candidateList_splitByTumor_models_notNoise.txt',quote=F,row.names=F, sep="\t")
rm(driverTypes,driverTypesFile)

setwd(workingDir)

# 1.2 - join all candidate lists ====
setwd("switches")

# calculate aggregated table (sum patients, etc.) and calculate unbalance as the minimun p of enrichment in a fisher test
candidatesDf_agg <- ddply(candidatesDf,.(GeneId,Symbol,Normal_transcript,Tumor_transcript,
                                         Normal_protein,Tumor_protein,Annotation,DriverAnnotation,IsRelevant,
                                         Driver,DriverType,Druggable,CDS,CDS_change,UTR_change),
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
candidatesDf_agg <- candidatesDf_agg[,c("GeneId","Symbol","Normal_transcript","Tumor_transcript",
                                        "Normal_protein","Tumor_protein","Annotation",
                                        "DriverAnnotation","IsRelevant","Driver","DriverType",
                                        "Druggable","CDS","CDS_change","UTR_change",
                                        "CancerAffected","PatientNumber","Percentage",
                                        "p.unbalance","Patients","Entropy")]
write.table(candidatesDf_agg,'tables/candidateList_allCancers_models_notNoise.txt',quote=F,row.names=F, sep="\t")

setwd(workingDir)

# 1.3 - calculate most frequent genes altered ====
setwd("switches")

df <- candidatesDf
df$Symbol <- as.character(df$Symbol)
df$Symbol[ df$DriverAnnotation=="Driver" ] <- paste0("*",df$Symbol[ df$DriverAnnotation=="Driver" ],"*")

#most frequent genes
for (thing in c("# patients","% patients")){
  for (effect in c("allGenes","functionalGenes")){
    for (annotation in c("allGenes","Driver","d1")){
      ngenes <- 50
      
      if (annotation=="allGenes"){
        genesSelection <- matrix(T,nrow(df),1)
      } else {
        genesSelection <- df$DriverAnnotation==annotation
      }
      
      if (effect=="functionalGenes") { 
        functionalSelection <- as.logical(df$IsRelevant)
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
    for (annotation in c("allGenes","Driver","d1")){
      ngenes <- 50
      
      if (annotation=="allGenes"){
        genesSelection <- matrix(T,nrow(df),1)
      } else {
        genesSelection <- df$DriverAnnotation==annotation
      }
      
      if (effect=="functionalGenes") { 
        functionalSelection <- as.logical(df$IsRelevant)
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
setwd(workingDir)

# 1.4 - Study exon length ====
setwd("switches")

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
    stat_ecdf(data=subset(cancer.exons,Random=="FullRandom"), aes(x=CDSRelativeSize,colour="green"),show_guide = FALSE) + 
    stat_ecdf(data=subset(cancer.exons,Random=="NonRandom"), aes(x=CDSRelativeSize,colour="red"),show_guide = FALSE) +
    ggtitle(knsur) + 
    theme_bw() + 
    ylab("")
  
  # cds position
  cdsPosition[[knsur]] <- ggplot() + ggtitle(knsur) + theme_bw() + 
    stat_ecdf(data=subset(cancer.exons,Random=="NonRandom"), aes(x=Position,colour="red"),show_guide = FALSE) + 
    stat_ecdf(data=subset(cancer.exons,Random=="FullRandom"), aes(x=Position,colour="green"),show_guide = FALSE)
  
  # length
  exonLength[[knsur]] <- ggplot() + ggtitle(knsur) + theme_bw() + 
    stat_ecdf(data=subset(cancer.exons,Random=="NonRandom"), aes(x=Length,colour="red"),show_guide = FALSE) + 
    stat_ecdf(data=subset(cancer.exons,Random=="FullRandom"), aes(x=Length,colour="green"),show_guide = FALSE)
  
  # orf
  switches <- table(cancer.exons$KeepOrf[cancer.exons$Random=="NonRandom"])
  randomSwitches <- table(cancer.exons$KeepOrf[cancer.exons$Random=="FullRandom"])
  
  cTable <- rbind(switches,randomSwitches)
  f <- fisher.test(cTable)
  this.Data <- data.frame(cancer=knsur,p=f$p.value,oddsRatio=f$estimate)
  orfChange <- rbind(orfChange,this.Data)
  
  # type
  switches <- table(cancer.exons$Type[cancer.exons$Random=="NonRandom" ])
  randomSwitches <- table(cancer.exons$Type[cancer.exons$Random=="FullRandom" ])
  
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
setwd(workingDir)

# 1.5 - Study exons per switch ====
# setwd("switches")
# 
# exonsPerSwitch_rel <- read.delim("exonsPerSwitch.tsv", header=FALSE)
# exonsPerSwitch_rel <- cbind(exonsPerSwitch_rel,"Relevant")
# colnames(exonsPerSwitch_rel) <- c("Counts","Cancer","Switch","Relevance")
# 
# exonsPerSwitch_nrel <- read.delim("exonsPerSwitch_nonrelevant.tsv", header=FALSE)
# exonsPerSwitch_nrel <- cbind(exonsPerSwitch_nrel,"NonRelevant")
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
# 
# setwd(workingDir)

# 1.6 - Study isoform length ====
setwd("switches")

# do not use 0 in analyses
isoLengths <- read.delim("isoform_length.tsv", header=T)
isoLengths$nIsoLength[isoLengths$nIsoLength==0] <- NA
isoLengths$tIsoLength[isoLengths$tIsoLength==0] <- NA

for (r in c("NonRandom","FullRandom")){
  
  this.isoLengths <- isoLengths[isoLengths$Random==r,]
  
  p <- c()
  
  for (cancer in cancerTypes){
    cancer.isoLengths <- this.isoLengths[this.isoLengths$Cancer==cancer,]
    t <- t.test(cancer.isoLengths$nIsoLength,cancer.isoLengths$tIsoLength,paired=TRUE,alternative="greater")
    p <- c(p,t$p.value)
  }
  
  isoLengthsSummary <- ddply(subset(this.isoLengths,!is.na(nIsoLength) & !is.na(tIsoLength)),
                             .(Cancer), summarise, 
                             mean_nIso=mean(nIsoLength),
                             mean_tIso=mean(tIsoLength),
                             mean_diff=mean(nIsoLength-tIsoLength))
  isoLengthsSummary <- cbind(isoLengthsSummary,p)
  
  write.table(isoLengthsSummary,file=paste0("tables/isoform_length_",r,"_p.tsv"),sep="\t", row.names=F, col.names=F, quote=F)
  
  isoLengths.Melt <- data.frame(Cancer=c(as.character(this.isoLengths$Cancer),as.character(this.isoLengths$Cancer)),
                                Length=c(this.isoLengths$nIsoLength,this.isoLengths$tIsoLength),
                                IsoformOrigin=c(rep("Normal",nrow(this.isoLengths)),rep("Tumor",nrow(this.isoLengths))))
  
  # mean difference
  mean(this.isoLengths$nIsoLength-this.isoLengths$tIsoLength,na.rm=T)
  
  ranges <- isoLengths.Melt
  colnames(ranges) <- c("Cancer","y","Categories")
  ranges <- ranges[!is.na(ranges$y),]
  
  sigElements <- getBoxplotAsterisks(isoLengthsSummary,ranges)
  
  # plot isoform length
  p <- ggplot(isoLengths.Melt) + 
    geom_boxplot(aes(factor(Cancer),Length,fill=IsoformOrigin),outlier.colour = NA) + 
    smartas_theme() + 
    scale_y_continuous(limits=c(0,1500)) + 
    xlab("Cancer") + 
    scale_fill_manual(values=c("#d95f02","#7570b3")) +
    theme(legend.position="bottom")
  
  if (! empty(sigElements$asterisks) || ! empty(sigElements$arcs) ){
    p <- p + geom_text(data=sigElements$asterisks ,aes(x=x_ast,y=y_ast,label='*'),size=15) + 
      geom_line(data=sigElements$arcs, aes(x_arc,y_arc,group=Cancer))
  }
    
  ggsave(paste0("figures/isoformLength_",r,".png"),p, width = 6.5, height = 5)
  
}

rm(p,t,sigElements,ranges,isoLengths,isoLengths.Melt,r,
   isoLengthsSummary,cancer.isoLengths,this.isoLengths)
setwd(workingDir)

# # 1.6.1 - Study normal isoform length ####
setwd("switches")

# #niso_length_all <- data.frame()
# niso_length_rel <- data.frame()
# niso_length_nrel <- data.frame()
# 
# for (cancer in cancerTypes){
#   basePath <- ""
#   #cancer_all <- read.delim(paste0(basePath,cancer,".tIso_length.tsv"), header=FALSE)
#   #tiso_length_all <- rbind(tiso_length_all,cbind(cancer_all,cancer))
#   
#   cancer_rel <- read.delim(paste0(basePath,cancer,".nIso_length_relevant.tsv"), header=FALSE)
#   cancer_rel <- cbind(cancer_rel,"Relevant",cancer)
#   colnames(cancer_rel) <- c("Length","Relevance","Cancer")
#   niso_length_rel <- rbind(niso_length_rel,cancer_rel)
#   
#   cancer_nrel <- read.delim(paste0(basePath,cancer,".nIso_length_nonrelevant.tsv"), header=FALSE)
#   cancer_nrel <- cbind(cancer_nrel,"NonRelevant",cancer)
#   colnames(cancer_nrel) <- c("Length","Relevance","Cancer")
#   niso_length_nrel <- rbind(niso_length_nrel,cancer_nrel)
# }
# 
# niso_length <- rbind(niso_length_rel,niso_length_nrel)
# 
# png("niso_length.png", width=1000, height=800)
# par(mar=c(9,3,5,2))
# boxplot(Length~interaction(Relevance,Cancer),data=niso_length,outline=F,col=c("gray40","gray90"),las=2)
# graphics.off()
# 
setwd(workingDir)

# # 1.6.2 - Study tumor isoform length ####
setwd("switches")
# #tiso_length_all <- data.frame()
# tiso_length_rel <- data.frame()
# tiso_length_nrel <- data.frame()
# 
# for (cancer in cancerTypes){
#   basePath <- ""
#   #cancer_all <- read.delim(paste0(basePath,cancer,".tIso_length.tsv"), header=FALSE)
#   #tiso_length_all <- rbind(tiso_length_all,cbind(cancer_all,cancer))
#   
#   cancer_rel <- read.delim(paste0(basePath,cancer,".tIso_length_relevant.tsv"), header=FALSE)
#   cancer_rel <- cbind(cancer_rel,"Relevant",cancer)
#   colnames(cancer_rel) <- c("Length","Relevance","Cancer")
#   tiso_length_rel <- rbind(tiso_length_rel,cancer_rel)
#   
#   cancer_nrel <- read.delim(paste0(basePath,cancer,".tIso_length_nonrelevant.tsv"), header=FALSE)
#   cancer_nrel <- cbind(cancer_nrel,"NonRelevant",cancer)
#   colnames(cancer_nrel) <- c("Length","Relevance","Cancer")
#   tiso_length_nrel <- rbind(tiso_length_nrel,cancer_nrel)
# }
# 
# tiso_length <- rbind(tiso_length_rel,tiso_length_nrel)
# 
# png("tiso_length.png", width=1000, height=800)
# par(mar=c(9,3,5,2))
# boxplot(Length~interaction(Relevance,Cancer),data=tiso_length,outline=F,col=c("gray40","gray90"),las=2)
# graphics.off()

setwd(workingDir)

# 1.7 - Generate plot about relevance of the switches ====
setwd("switches")

candidatesDf_copy <- candidatesDf
drivers <- candidatesDf_copy$DriverAnnotation=="Driver"
candidatesDf_copy$Driver <- as.character(candidatesDf_copy$Driver)
candidatesDf_copy$Driver[drivers] <- "Unlabeled driver"
candidatesDf_copy$Driver[!drivers] <- "NonDriver"
candidatesDf_copy$Driver[candidatesDf_copy$DriverType == "oncogene"] <- "Oncogene"
candidatesDf_copy$Driver[candidatesDf_copy$DriverType == "suppressor"] <- "Suppressor"
candidatesDf_copy$Driver <- factor(candidatesDf_copy$Driver,c("NonDriver","Unlabeled driver","Oncogene","Suppressor"))

candsStats <- ddply(candidatesDf_copy,.(Tumor,IsRelevant,Driver),summarise,Count=length(GeneId))
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

rm(p,candsStats,func,candidatesDf_copy,drivers,cancer)
setwd(workingDir)

# 1.8 - Plot number of transcripts vs number of patients affected by a switch ====
setwd("switches")

txsAndGenes <- read.delim(paste0(workingDir,"/Data/TCGA/genesAndTranscripts_splitByGene.txt"))

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
setwd(workingDir)

# 1.9 - Print proportion of switches in drivers and non-drivers ====
setwd("switches")

for (a in c("d0_enrichment","d0_relevant_enrichment","d1_enrichment","d1_relevant_enrichment")){

  candsStats <- read.delim(paste0(a,".tsv"))
  colnames(candsStats) <- c("Cancer","Analysis","SwitchedDriver","NonSwitchedDriver","SwitchedNonDrivers","NonSwitchedNonDrivers","p","OddsRatio")
  candsStats_driver <- melt(candsStats[,c("Cancer","SwitchedDriver","NonSwitchedDriver")],value.name="Counts")
  candsStats_driver$variable <- as.character(candsStats_driver$variable)
  candsStats_driver$variable[grep("NonSwitched",candsStats_driver$variable)] <- "No switch"
  candsStats_driver$variable[grep("Switched",candsStats_driver$variable)] <- "Switch"
  candsStats_driver$Driver <- "Driver"
  
  candsStats_nonDriver <- melt(candsStats[,c("Cancer","SwitchedNonDrivers","NonSwitchedNonDrivers")],value.name="Counts")
  candsStats_nonDriver$variable <- as.character(candsStats_nonDriver$variable)
  candsStats_nonDriver$variable[grep("NonSwitched",candsStats_nonDriver$variable)] <- "No switch"
  candsStats_nonDriver$variable[grep("Switched",candsStats_nonDriver$variable)] <- "Switch"
  candsStats_nonDriver$Driver <- "NonDriver"
  
  geneInfo <- rbind(candsStats_driver,candsStats_nonDriver)
  
  totalGenes <- ddply(geneInfo,.(Cancer,Driver),summarise,Total=sum(Counts))
  
  geneInfo <- merge(geneInfo,totalGenes)
  geneInfo$NormCounts <- geneInfo$Counts/geneInfo$Total
  geneInfo$Switched <- factor(geneInfo$variable,c("Switch","No switch"))
  geneInfo$Cancer2 <- paste(geneInfo$Cancer,geneInfo$Driver,sep="_")
  geneInfo$Driver <- factor(geneInfo$Driver,c("Driver","NonDriver"))
  
  geneInfo <- geneInfo[order(geneInfo$Switched),]
  geneInfo <- geneInfo[order(geneInfo$Driver),]
  geneInfo$y <- geneInfo$NormCounts
  geneInfo$variable <- geneInfo$Driver
  
  sigElements <- getBarplotAsterisks(candsStats,subset(geneInfo,Switched=="Switch"))
  
  p <- ggplot(subset(geneInfo,Switched=="Switch")) + 
    geom_bar(stat="identity",aes(x=Cancer2,y=NormCounts,fill=Cancer,alpha=Driver)) + 
    smartas_theme() + 
    ylab("Proportion of total with a switch") + 
    xlab("Tumor type") + 
    theme(text = element_text(size=20),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour="black")) + 
    scale_fill_manual(values=colorPalette) + 
    scale_alpha_manual(values=c(1,0.7)) + 
    scale_x_discrete(labels=c("Driver","NonDriver")) + 
    theme(axis.text.x = element_blank()) + 
    geom_text(data=sigElements$asterisks ,aes(x=x_ast,y=y_ast,label='*'),size=15) + 
    geom_line(data=sigElements$arcs, aes(x_arc,y_arc,group=Cancer))

  ggsave(paste0("figures/",a,".png"),p, width=8, height=7)
  
}

candidatesDf_copy <- candidatesDf
relevantSwitches <- as.logical(candidatesDf_copy$IsRelevant)
candidatesDf_copy$IsRelevant[relevantSwitches] <- "Functional"
candidatesDf_copy$IsRelevant[!relevantSwitches] <- "Non-functional"
candidatesDf_copy$IsRelevant <- factor(candidatesDf_copy$IsRelevant,c("Functional","Non-functional"))

p <- numeric()
for (cancer in cancerTypes){
  thisDf <- candidatesDf_copy[candidatesDf_copy$Tumor==cancer,]
  selectRelevant <- thisDf$IsRelevant == "Functional"
  k <- ks.test(thisDf$NumPatients[selectRelevant],thisDf$NumPatients[!selectRelevant],alternative="less")
  p <- c(p,k$p.value)
}

signif <- data.frame(Cancer=cancerTypes,p=p)

plotInfo <- candidatesDf_copy[,c("Tumor","NumPatients","IsRelevant")]
colnames(plotInfo) <- c("Cancer","y","Categories")

sigElements <- getBoxplotAsterisks(signif,plotInfo)

# plot isoform length
p <- ggplot(candidatesDf_copy) + 
  geom_boxplot(aes(factor(Tumor),NumPatients,fill=IsRelevant),outlier.colour = NA) + 
  scale_fill_manual(values=c("#d95f02","#7570b3")) + 
  xlab("Cancer") + ylab("Patient number") + scale_y_continuous(limits=c(0,175)) +
  smartas_theme() + theme(legend.position="bottom") +
  geom_text(data=sigElements$asterisks ,aes(x=x_ast,y=y_ast,label='*'),size=15) +
  geom_line(data=sigElements$arcs, aes(x_arc,y_arc,group=Cancer))

ggsave("figures/length_comparison_functionalVsNonFunctional.png",p, width=8, height=7)

rm(a,p,k,relevantSwitches,selectRelevant,cancer,sigElements,
   totalGenes,thisDf,signif,plotInfo,geneInfo,candsStats,
   candsStats_driver,candsStats_nonDriver,candidatesDf_copy)
setwd(workingDir)

# 2 - NEIGHBORHOODS ---------------------

# 2.1 - Count number of times a geneset is significantly altered in different cancer types ====
setwd('neighborhood_analysis')
genesetType <- c("canonical_pathways","hallmarks","go_biological_process","oncogenic_signatures")

for (gnset in c("relevant","all")){
  for (type in genesetType){
    file <- paste(type,gnset,sep='_')
    sets_raw <- read.delim(paste0(file,'.txt'),header=TRUE)
    set_counts <- table(sets_raw$GeneSet[sets_raw$qval<0.05])
    set_counts <- set_counts[order(-set_counts)]
    set_counts_df <- as.data.frame(set_counts)
    set_counts_df$Geneset <- row.names(set_counts_df)
    set_counts_df$Counts <- set_counts_df$set_counts
    set_counts_df <- set_counts_df[,c("Geneset","Counts")]
    
    write.table(set_counts_df,paste0('tables/',file,'_counts.txt'),row.names=F, quote=F,sep="\t")
    
    # plot affected in more than 8
    set_counts_df$Geneset <- splitTextInLines(set_counts_df$Geneset)
    set_counts_df$Geneset <- factor(set_counts_df$Geneset, levels=set_counts_df$Geneset)
    p <- ggplot(subset(set_counts_df,Counts>8)) + 
      geom_bar(stat="identity",aes(x=Geneset,y=Counts)) + 
      smartas_theme() + 
      theme(axis.text.x=element_text(size=5,angle=90,hjust=1,vjust=0.5,colour="black"))
    ggsave(paste0("figures/",type,"_",gnset,"_moreThan8TumorTypes.png"),p, width = 6.5, height = 5)
    
    set_counts_df$simpleGeneset <- set_counts_df$Geneset
    set_counts_df$simpleGeneset <- gsub("\n"," ",set_counts_df$simpleGeneset)
    set_counts_df$simpleGeneset <- gsub("REACTOME ","",set_counts_df$simpleGeneset)
    set_counts_df$simpleGeneset <- gsub("BIOCARTA ","",set_counts_df$simpleGeneset)
    set_counts_df$simpleGeneset <- gsub("KEGG ","",set_counts_df$simpleGeneset)
    set_counts_df$simpleGeneset <- gsub("PID ","",set_counts_df$simpleGeneset)
    set_counts_df$simpleGeneset <- gsub("HALLMARK ","",set_counts_df$simpleGeneset)
    
    pal <- brewer.pal(9, "BuGn")
    png(paste0("figures/wordcloud_",type,"_",gnset,".png"), width=3000,height=2500)
    wordcloud(set_counts_df$simpleGeneset,set_counts_df$Counts,colors=pal,scale=c(4,.2))
    dev.off()
    
  }
}

rm(p,genesetType,gnset,type,file,sets_raw,set_counts,set_counts_df)
setwd(workingDir)

# kk <- read.delim(paste0("hallmark_info/","HALLMARK_ADIPOGENESIS_onlyDrivers.tsv"))

# 3 - MUTATION OVERLAP ---------------------

# 3.1 - calculate allCancers p.values for each combination of mutation and switch types ====
setwd('mutations')

muts <- c("all_mutations","functional_mutations")
swis <- c("all_switches","functional_switches")
analyses <- c("gene","geneset","pannegative")
usedSet <- c("","onlyDrivers")

for (a in analyses){
  for (m in muts){
    for (s in swis){
      for (u in usedSet){
        if (u=="onlyDrivers" && a %in% c("gene","pannegative") ){
          next
        }
        
        if (u=="onlyDrivers"){
          analysisFile <- paste0(paste(a,m,s,u,sep="_"),'.txt')
          outFile <- paste("tables",paste(a,m,s,u,'allCancers.txt',sep="_"),sep="/")
        } else {
          analysisFile <- paste0(paste(a,m,s,sep="_"),'.txt')
          outFile <- paste("tables",paste(a,m,s,'allCancers.txt',sep="_"),sep="/")
        }
        
        analysis <- read.delim(analysisFile, header=TRUE)
        
        if (a %in% c("gene","pannegative")){
          analysis$Hallmark <- ""
        }
        
        analysis <- merge(analysis,nPatientsDf)
        
        analysis$ms[analysis$ms==0] <- 0.000001
        analysis$s[analysis$s==0] <- 0.000001
        analysis_agg <- ddply(analysis,.(Gene,Symbol,nTx,tTx,Hallmark), summarise, MS=sum(ms),M=sum(m),S=sum(s),N=sum(n),H=-sum(((ms+s)/TotalPatients)/sum((ms+s)/TotalPatients)*log2(((ms+s)/TotalPatients)/sum((ms+s)/TotalPatients)))/log2(length(ms)))
        
        analysis_agg$MS <- round(analysis_agg$MS)
        analysis_agg$S <- round(analysis_agg$S)
        
        me <- apply(analysis_agg,1, function(x){ 
          f <- fisher.test(x=matrix(as.numeric(x[6:9]),nrow=2,ncol=2),alternative="less")
          data.frame(p=f$p.value,OR=f$estimate) } )
        
        me <- do.call("rbind",me)
        analysis_agg$p_me <- me$p
        analysis_agg$p.adj_me <- p.adjust(analysis_agg$p_me)

        o <- apply(analysis_agg,1, function(x){ 
          f <- fisher.test(x=matrix(as.numeric(x[6:9]),nrow=2,ncol=2),alternative="greater")
          data.frame(p=f$p.value,OR=f$estimate) } )
        
        o <- do.call("rbind",o)
        analysis_agg$p_o <- o$p
        analysis_agg$p.adj_o <- p.adjust(analysis_agg$p_o)
        
        analysis_agg$OR <- o$OR
        
        analysis_agg <- analysis_agg[order(analysis_agg$p_me),]
        
        if (a %in% c("gene","pannegative")){
          analysis_agg <- analysis_agg[,!colnames(analysis_agg) %in% c("Hallmark")]
        } else if (a == "geneset") {
          geneset_score <- ddply(analysis_agg,.(Hallmark), summarise, Score=max(p_me),AdjustedScore=max(p.adj_me))
        }
        
        write.table(analysis_agg,outFile,quote=F,col.names=T,sep="\t",row.names=FALSE)
      }
    }
  }
}

setwd(workingDir)

# 3.2 - calculate gene meScore ====
setwd('mutations')

gn_all <- read.delim(paste("tables/gene","all_mutations","functional_switches",'allCancers.txt',sep="_"), header=TRUE)
gn_fun <- read.delim(paste("tables/gene","functional_mutations","functional_switches",'allCancers.txt',sep="_"), header=TRUE)

gn_merged <- merge(gn_fun,gn_all,by=c("Gene","Symbol","nTx","tTx"))
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
write.table(gn_merged[,c("Gene","Symbol","nTx","tTx","Score","Sign" )],'tables/scoredMEGenes.txt',quote=F,col.names=T,sep="\t",row.names=FALSE)

p <- ggplot(gn_merged, aes(x=Order,y=NewScore,fill = Sign)) + 
  geom_area(alpha=0.75) + smartas_theme() + 
  labs(title="Genes ranked by score", x="Genes", y="signed sqrt(Score)") + 
  geom_hline(yintercept=0, size=0.4, color="black") + 
  scale_fill_manual(values=c("Mutual exclusion"="red","Nothing"="black","Coincidence"="darkblue")) + 
  theme(axis.text.x=element_blank())
ggsave("figures/meScores.png",p)

setwd(workingDir)

# 3.3 - calculate geneset meScore ====
setwd('mutations')

gnset_all <- read.delim(paste("tables/geneset","all_mutations","functional_switches",'onlyDrivers','allCancers.txt',sep="_"), header=TRUE)
gnset_fun <- read.delim(paste("tables/geneset","functional_mutations","functional_switches",'onlyDrivers','allCancers.txt',sep="_"), header=TRUE)

gnset_merged <- merge(gnset_fun,gnset_all,by=c("Gene","Symbol","nTx","tTx","Hallmark"))
gnset_merged$Score <- -log10(gnset_merged$p_me.x)+log10(gnset_merged$p_o.y)

gnset_merged$Sign[gnset_merged$Score > 0] <- "Mutual exclusion"
gnset_merged$Sign[gnset_merged$Score < 0] <- "Coincidence"
gnset_merged$Sign[gnset_merged$Score == 0] <- "Nothing"

gnset_merged$NewScore <- ifelse(gnset_merged$Score > 0,sqrt(gnset_merged$Score),-sqrt(-gnset_merged$Score))
gnset_merged$NewScore <- ifelse(gnset_merged$Score == 0,0,gnset_merged$NewScore)
gnset_merged <- gnset_merged[order(-gnset_merged$NewScore),]
gnset_merged$Symbol <- factor(gnset_merged$Symbol,levels=as.factor(unique(gnset_merged$Symbol)))
gnset_merged$Order=1:nrow(gnset_merged)

write.table(gnset_merged[,c("Gene","Symbol","Hallmark","Score","Sign" )],'tables/scoredMEGenesets.txt',quote=F,col.names=T,sep="\t",row.names=FALSE)

p <- ggplot(gnset_merged, aes(x=Order,y=NewScore,fill = Sign)) + 
  geom_area(alpha=0.75) + smartas_theme() + 
  labs(title="Genes ranked by score", x="Genes", y="signed sqrt(Score)") + 
  geom_hline(yintercept=0, size=0.4, color="black") + 
  scale_fill_manual(values=c("Mutual exclusion"="red","Nothing"="black","Coincidence"="darkblue")) + 
  theme(axis.text.x=element_blank())
ggsave("figures/meGenesetScores.png",p)

setwd(workingDir)

# 3.4 - calculate p.value of overlap between mutation in a feature and switches affecting it ====
setwd('mutations')

mut_feat_overlap <- read.delim("mutation_switch_feature_overlap.txt", header=TRUE)
pvals <- apply(mut_feat_overlap[,c("MutationsInFeature","TotalMutations","FeatureSize")],1, function(x){ 
  if (x[2]==0){
    p <- 1
  } else {
    p <- binom.test(x[1],x[2],x[3],"greater")
    p <- p$p.value
  }
  p
} )
adjpvalues <- p.adjust(pvals)

mut_feat_overlap$p <- pvals
mut_feat_overlap$adjp <- adjpvalues
mut_feat_overlap <- mut_feat_overlap[order(mut_feat_overlap$p),]

write.table(mut_feat_overlap,'tables/mutation_switch_feature_overlap_withPVals.txt',quote=F,col.names=T,sep="\t",row.names=FALSE)

# create the allCancers table
mut_feat_overlap_agg <- ddply(mut_feat_overlap,
                              .(Gene,Symbol,Normal_transcript,Tumor_transcript,What,
                                FeatureType,Feature,Driver,FeatureSize), 
                              summarise, inMut=sum(MutationsInFeature), 
                              totalMut=sum(TotalMutations) )
mut_feat_overlap_agg$Ratio = 100 * mut_feat_overlap_agg$inMut/mut_feat_overlap_agg$totalMut

pvals <- apply(mut_feat_overlap_agg[,c("inMut","totalMut","FeatureSize")],1, function(x){ 
  if (x[2]==0){
    p <- 1
  } else {
    p <- binom.test(x[1],x[2],x[3],"greater")
    p <- p$p.value
  }
  p
} )
adjpvalues <- p.adjust(pvals)

mut_feat_overlap_agg$p_mutation_feature_overlap <- pvals
mut_feat_overlap_agg$adjp_mutation_feature_overlap <- adjpvalues
mut_feat_overlap_agg <- mut_feat_overlap_agg[order(mut_feat_overlap_agg$p_mutation_feature_overlap),]

write.table(mut_feat_overlap_agg,'tables/mutation_switch_feature_overlap_allCancers_withPVals.txt',quote=F,col.names=T,sep="\t",row.names=FALSE)

p <- ggplot() + 
  geom_point(data=mut_feat_overlap_agg,
             aes(FeatureSize*100,Ratio)) + 
  geom_point(data=subset(mut_feat_overlap_agg,p_mutation_feature_overlap<0.05),
             aes(FeatureSize*100,Ratio),size=3,color="#7570b3") + 
  geom_point(data=subset(mut_feat_overlap_agg,adjp_mutation_feature_overlap<0.05),
             aes(FeatureSize*100,Ratio),size=4.5,color="#d95f02") + 
  smartas_theme() + 
  xlab("Feature size (as % of the total length)") + 
  ylab("#mutations (as % of the total in the gene)")
ggsave("figures/affectedFeatures_mutations_vs_switches.png",p)

setwd(workingDir)

# 3.5 - Plot p_me in a gene against the p_overlap of mutations and features ====
setwd('mutations')

gn_fun <- read.delim(paste("tables/gene","functional_mutations","functional_switches",'allCancers.txt',sep="_"), header=TRUE)
gn_fun <- gn_fun[,c("Gene","Symbol","p_me")]

mut_feat_overlap_agg <- merge(mut_feat_overlap_agg,gn_fun,by=c("Gene","Symbol"))
mut_feat_overlap_agg$PlotId <- paste(mut_feat_overlap_agg$Symbol,mut_feat_overlap_agg$Feature,sep="-")

# drivers
p <- ggplot(subset(mut_feat_overlap_agg, Driver=="True"),aes(-log10(p_me),-log10(p_mutation_feature_overlap))) + 
  geom_point() + 
  geom_point(data=subset(mut_feat_overlap_agg,Driver=="True" & p_mutation_feature_overlap < 0.5 & p_me < 0.5),aes(-log10(p_me),-log10(p_mutation_feature_overlap),color=PlotId)) + 
  smartas_theme()
p <- direct.label(p)
ggsave("figures/ME_vs_Overlap_onlyDrivers.png",p)

# all
p <- ggplot(mut_feat_overlap_agg,aes(-log10(p_me),-log10(p_mutation_feature_overlap))) + 
  geom_point() + 
  geom_point(data=subset(mut_feat_overlap_agg, p_mutation_feature_overlap < 0.2 & p_me < 0.2),aes(-log10(p_me),-log10(p_mutation_feature_overlap),color=PlotId)) + 
  smartas_theme()
p <- direct.label(p)
ggsave("figures/ME_vs_Overlap.png",p)

# filter by positive MEScore
mut_feat_overlap_agg_filt <- mut_feat_overlap_agg[mut_feat_overlap_agg$Gene %in% gn_merged$Gene[gn_merged$Sign == "Mutual exclusion"],]

p <- ggplot(mut_feat_overlap_agg_filt,aes(-log10(p_me),-log10(p_mutation_feature_overlap))) + 
  geom_point() + 
  geom_point(data=subset(mut_feat_overlap_agg_filt, p_mutation_feature_overlap < 0.2 & p_me < 0.2),aes(-log10(p_me),-log10(p_mutation_feature_overlap),color=PlotId)) + 
  smartas_theme()
p <- direct.label(p)
ggsave("figures/ME_vs_Overlap_onlyPositiveMescore.png",p)

# p <- p + geom_text(data=subset(mut_feat_overlap_agg,MutationScore >= 2 & MEScore > 0), aes(MEScore,MutationScore,label=Symbol))
# 
# annotate("text", x=c(-0.2,0.9), y=c(33,32), label = c("my label", "label 2"))

mut_feat_overlap_agg <- mut_feat_overlap_agg[order(-mut_feat_overlap_agg$Ratio),]

write.table(mut_feat_overlap_agg,'tables/mutation_switch_feature_overlap_allCancers2.txt',quote=F,col.names=T,sep="\t",row.names=FALSE)

setwd(workingDir)

# 3.6 - feature overlap ====
setwd('mutations')

feat_enrich <- read.delim('feature_enrichment.txt', header=TRUE)
feat_enrich <- feat_enrich[feat_enrich$DomainFrequency != 0,]
feat_enrich$MutTotal <- feat_enrich$MutIn + feat_enrich$MutOut
feat_enrich$SwitchesTotal <- feat_enrich$SwitchesIn + feat_enrich$SwitchesOut
feat_enrich$MutFreq <- feat_enrich$MutIn/feat_enrich$MutTotal
feat_enrich$SwitchFreq <- feat_enrich$SwitchesIn/feat_enrich$SwitchesTotal

for (cancer in cancerTypes){
  
  cancer.feat_enrich <- feat_enrich[feat_enrich$Cancer==cancer,]
  
  # plot normalized mutation frequency and switched frequency
  minMut <- quantile(cancer.feat_enrich$MutFreq/cancer.feat_enrich$DomainFrequency,0.95)
  minSwitch <- quantile(cancer.feat_enrich$SwitchFreq/cancer.feat_enrich$DomainFrequency,0.95)
  

  # calculate regression coefficient
  r.df <- cancer.feat_enrich[cancer.feat_enrich$MutFreq > 0 & cancer.feat_enrich$SwitchFreq > 0,]
  r <- cor(log2(r.df$MutFreq/r.df$DomainFrequency),log2(r.df$SwitchFreq/r.df$DomainFrequency))
  
  p <- ggplot(cancer.feat_enrich,aes(log2(MutFreq/DomainFrequency),log2(SwitchFreq/DomainFrequency))) + 
    geom_point() + 
    geom_point(data=subset(cancer.feat_enrich, MutFreq/DomainFrequency > minMut & SwitchFreq/DomainFrequency > minSwitch),aes(log2(MutFreq/DomainFrequency),log2(SwitchFreq/DomainFrequency),color=Domain)) + 
    smartas_theme() + 
    geom_smooth(data=subset(cancer.feat_enrich,MutFreq>0 & SwitchFreq>0),method=lm) + 
    geom_text(x=4,y=5,label=paste0("R = ",round(r,2)))
  p <- direct.label(p)
  
  ggsave(paste0("figures/",cancer,".domain_enrichment.png"),p)
  
}

# aggregate cancers
feat_enrich_agg <- ddply(feat_enrich,.(Domain), summarise, 
                        MutIn=sum(MutIn), MutTotal=sum(MutTotal),
                        SwitchesIn=sum(SwitchesIn), SwitchesTotal=sum(SwitchesTotal),
                        MeanDomainFrequency=mean(DomainFrequency))
feat_enrich_agg$MutFreq <- feat_enrich_agg$MutIn/feat_enrich_agg$MutTotal
feat_enrich_agg$SwitchFreq <- feat_enrich_agg$SwitchesIn/feat_enrich_agg$SwitchesTotal

minMut <- quantile(feat_enrich_agg$MutFreq/feat_enrich_agg$MeanDomainFrequency,0.95)
minSwitch <- quantile(feat_enrich_agg$SwitchFreq/feat_enrich_agg$MeanDomainFrequency,0.95)

# calculate regression coefficient
r.df <- feat_enrich_agg[feat_enrich_agg$MutFreq > 0 & feat_enrich_agg$SwitchFreq > 0,]
r <- cor(log2(r.df$MutFreq/r.df$MeanDomainFrequency),log2(r.df$SwitchFreq/r.df$MeanDomainFrequency))

p <- ggplot(data=feat_enrich_agg,
            aes(log2(MutFreq/MeanDomainFrequency),log2(SwitchFreq/MeanDomainFrequency))) + 
  geom_point() + 
  geom_point(data=subset(feat_enrich_agg,
                         MutFreq/MeanDomainFrequency > minMut & 
                         SwitchFreq/MeanDomainFrequency > minSwitch),
             aes(log2(MutFreq/MeanDomainFrequency),
                 log2(SwitchFreq/MeanDomainFrequency),
                 color=Domain)) + 
  smartas_theme() + 
  geom_smooth(data=subset(feat_enrich_agg,MutFreq>0 & SwitchFreq>0),method=lm) + 
  geom_text(x=4,y=5,label=paste0("R = ",round(r,2)))
p <- direct.label(p)

ggsave("figures/domain_enrichment_allCancers.png",p)

setwd(workingDir)

# 3.7 - Search for something that makes the meScore+ special :)  ====
setwd('mutations')

x <- gn_merged[gn_merged$NewScore > 1,c("Gene","nTx","tTx")]
colnames(x) <- c("GeneId","Normal_transcript","Tumor_transcript")
y <- gn_merged[gn_merged$NewScore < -1,c("Gene","nTx","tTx")]
colnames(y) <- c("GeneId","Normal_transcript","Tumor_transcript")

studyGroups(x,y,candidatesDf_agg)

feat_enrich <- read.delim('feature_enrichment.txt', header=TRUE)
feat_enrich <- feat_enrich[feat_enrich$DomainFrequency != 0,]
feat_enrich$MutTotal <- feat_enrich$MutIn + feat_enrich$MutOut
feat_enrich$SwitchesTotal <- feat_enrich$SwitchesIn + feat_enrich$SwitchesOut
feat_enrich$MutFreq <- feat_enrich$MutIn/feat_enrich$MutTotal
feat_enrich$SwitchFreq <- feat_enrich$SwitchesIn/feat_enrich$SwitchesTotal

structural_features <- read.delim(paste0(workingDir,"/structural_analysis/structural_features.onlyModels.tsv"))
structural_features <- structural_features[structural_features$Random=="NonRandom",]
structural_features$Feature2 <- gsub(" ","_",structural_features$Feature)

kk <- merge(structural_features,feat_enrich,by.x=c("Cancer","Feature2"),by.y=c("Cancer","Domain"))
zz <- merge(kk,gn_merged)

setwd(workingDir)

# 4 - STRUCTURAL ANALYSIS ---------------------
# 4.1 - i3d plots ====
setwd("structural_analysis")

i3d <- read.delim("i3d_broken.tsv", header=TRUE,row.names=NULL)

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
ggsave("figures/i3d_InteractionVsAffection.png",p)

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
p <- p + smartas_theme() + labs(x="Percentage of the interface lost", y="Frequency")
p <- p + geom_hline(yintercept=0, size=0.4, color="black") + scale_x_continuous(breaks=seq(0,100, by=5))
ggsave("figures/i3d_interactionAffected_completeStructures.png",p)

freqIntx = as.data.frame(table(df[,c("Annotation","PartnerAnnotation")]))
freqIntx <- freqIntx[order(-freqIntx$Freq),]
write.table(freqIntx,'tables/i3d_AnnotationFrequency.txt',quote=F,col.names=T,sep="\t",row.names=FALSE)

nrow(i3d_complete)
table(i3d_complete$WhatsHappening)

scoredTable <- merge(i3d,gn_merged[,colnames(gn_merged) %in% c("Gene","Score")],by=c("Gene"))
scoredTable <- rename(scoredTable,c("Score"="GeneScore"))
scoredTable <- merge(scoredTable,gnset_merged[,colnames(gnset_merged) %in% c("Gene","Score")],by=c("Gene"))
scoredTable <- rename(scoredTable,c("Score"="GenesetScore"))

filtScoredTable <- unique(scoredTable[scoredTable$SequenceCover>50 & scoredTable$PartnerCover>50 & scoredTable$InteractionAffection > 10,c("Symbol","Partner","InteractionAffection","Annotation","PartnerAnnotation","SequenceCover","PartnerCover","WhatsHappening","GeneScore","GenesetScore","Structure")])
filtScoredTable <- filtScoredTable[order(-filtScoredTable$GeneScore),]

filtScoredTable <- filtScoredTable[order(-filtScoredTable$GenesetScore),]

setwd(workingDir)

# 4.2 - domain enrichment in drivers,oncogenes,suppressors ====
setwd("structural_analysis")

structural_features <- read.delim("structural_features.onlyModels.tsv")
nrStructural_features <- structural_features[structural_features$Random == "NonRandom",]

for (featureType in c('Pfam','prosite')){
  
  annotation <- read.delim(paste0(workingDir,"/Data/Databases/",featureType,"2go.clean.txt"),row.names = NULL,header=F)
  colnames(annotation) <- c("id","GO")
  
  z <- gsub(";GO:[[:digit:]]+", "", annotation$GO)
  z <- gsub("^GO:", "", z)
  z <- paste0(toupper(substring(z, 1,1)),substring(z, 2))
  
  annotation$GO <- z
    
  for (action in c("Lost_in_tumor","Gained_in_tumor")){
    uniqStuff <- unique( nrStructural_features$Feature[nrStructural_features$Analysis==featureType & nrStructural_features$WhatsHappenning!="Nothing" & nrStructural_features$WhatsHappenning==action] )
    actionTag <- tolower(unlist(strsplit(action,split = " "))[1])
    
    for ( geneset in c("oncogene",'suppressor',"driver")){
      affectedFeats <- list()
      
      for (feat in uniqStuff){
        
        usedFeature <- nrStructural_features[nrStructural_features$WhatsHappenning==action & nrStructural_features$Feature==feat,]
        nonUsedFeature <- nrStructural_features[nrStructural_features$WhatsHappenning==action | nrStructural_features$Feature!=feat,]
        
        if (geneset=="driver"){
          FG=sum(usedFeature$Driver==1)
          FNG=sum(usedFeature$Driver==0)
          NFG=sum(nonUsedFeature$Driver==1)
          NFNG=sum(nonUsedFeature$Driver==0)
        } else {
          FG=sum(usedFeature$DriverType==geneset)
          FNG=sum(usedFeature$DriverType!=geneset)
          NFG=sum(nonUsedFeature$DriverType==geneset)
          NFNG=sum(nonUsedFeature$DriverType!=geneset)
        }
        
        affectedFeats[[feat]] <- data.frame(Feature=feat, WhatsHappenning=action,
                                            FG=FG,FNG=FNG,NFG=NFG,NFNG=NFNG)
      }
      
      affectedFeatsDf <- do.call('rbind',affectedFeats)
      
      affectedFeatsDf$p <- apply(affectedFeatsDf,1, function(x){ 
        k <- fisher.test(x=matrix(as.numeric(x[3:6]),nrow=2,ncol=2),alternative="greater")
        k$p.value } )
      affectedFeatsDf$p.adj <- p.adjust(affectedFeatsDf$p)
      affectedFeatsDf <- affectedFeatsDf[order(affectedFeatsDf$p),]
      kk <- merge(affectedFeatsDf,annotation,all.x=T)
      
      write.table(affectedFeatsDf,file=paste("tables",paste(featureType,actionTag,paste0(geneset,".tsv"),sep="_"),sep="/"),sep="\t", row.names=F, quote=F)
      
      # plot the significant groups
      affectedFeatsDf_s <- affectedFeatsDf[affectedFeatsDf$p.adj < 0.05,]
      
      if (nrow(affectedFeatsDf_s)==0){
        next
      }
      
      affectedFeatsDf_s$id <- unlist(strsplit(as.character(affectedFeatsDf_s$Feature),"|",fixed=T))[c(T,F)]
      affectedFeatsDf_s <- merge(affectedFeatsDf_s,annotation,all.x=T)
      
      annotationCount <- sort(table(as.character(affectedFeatsDf_s$GO)),decreasing=T)
      topGroups <- names(head(annotationCount[annotationCount > 1],4))
      
      plotAnnotatedFeats <- ddply(affectedFeatsDf_s,.(id,Feature,p.adj),function(y){ 
        x <- intersect(y$GO,topGroups)
        if (length(x) > 0){names(sort(-annotationCount[x])[1])} else {"Other"}})
      plotAnnotatedFeats$Annotation <- plotAnnotatedFeats$V1
      plotAnnotatedFeats$log_p <- -log10(plotAnnotatedFeats$p.adj)
      
      plotAnnotatedFeats <- plotAnnotatedFeats[order(-plotAnnotatedFeats$log_p,plotAnnotatedFeats$Feature),]
      plotAnnotatedFeats$Feature <- gsub("_"," ",unlist(strsplit(as.character(plotAnnotatedFeats$Feature),"|",fixed=T))[c(F,T)])
      
      plotAnnotatedFeats$Feature <- splitTextInLines(plotAnnotatedFeats$Feature)
      plotAnnotatedFeats$Feature <- factor(plotAnnotatedFeats$Feature,levels=as.factor(plotAnnotatedFeats$Feature))
      plotAnnotatedFeats$Annotation <- factor(plotAnnotatedFeats$Annotation,levels=as.factor(c(topGroups,"Other")))
      
      annotationPalette <- c("#d7191c","#fdae61","#ffffbf","#abdda4","#2b83ba")
      names(annotationPalette) <- c(topGroups,"Other")
      
      p <- ggplot(plotAnnotatedFeats,aes(x=Feature,y=log_p,fill=Annotation)) + 
        geom_bar(stat="identity",position="dodge") + 
        smartas_theme() + 
        scale_fill_manual(values=annotationPalette) + 
        xlab(paste0(featureType," features")) + ylab("-log10(adjusted p)") + 
        facet_wrap(~Annotation,drop=TRUE,scale="free_x") + 
        theme(axis.text.x=element_text(size=7,angle=90,hjust=1,vjust=0.5,colour="black"))
      ggsave(paste0("figures/",paste(featureType,actionTag,paste0(geneset,".png"),sep="_")),p)
    }
  }
}
  
setwd(workingDir)

# 4.3 - Study the enrichment in a particular set of features being gained or lost ====
setwd("structural_analysis")

analysis <- c("interpro","anchor","iupred","prosite")
analysisEnrichment <- data.frame(Cancer=c(),Analysis=c(),What=c(), p=c())
for (cancer in cancerTypes){
  for (tag in analysis){
    features <- read.delim(paste0(cancer,".",tag,"_analysis.tsv"),header=T)
    random_data <- read.delim(paste0(cancer,".",tag,"_analysis_random.tsv"),header=T)
    
    if ( tag %in% c("anchor","iupred") ){
      features <- features[as.logical(features$Significant),]
      random_data <- random_data[as.logical(random_data$Significant),]
    }    
    
    switches <- table(features$What)
    randomSwitches <- table(random_data$What)
    counts.features <- rbind(switches,randomSwitches)
    lostVsNothing <- fisher.test(counts.features[,-1])
    gainVsNothing <- fisher.test(counts.features[,-2])
    gainVsLost <- fisher.test(counts.features[,-3])
    
    this.features <- data.frame(Cancer=cancer,Analysis=tag,What=c("lostVsNothing","gainVsNothing","gainVsLost"), p=c(lostVsNothing$p.value,gainVsNothing$p.value,gainVsLost$p.value),oddsratio=c(lostVsNothing$estimate,gainVsNothing$estimate,gainVsLost$estimate))
    analysisEnrichment <- rbind(analysisEnrichment,this.features)
  }
}

write.table(analysisEnrichment,file="tables/analysis_enrichment.tsv",sep="\t", row.names=F, quote=F)

p <- ggplot(analysisEnrichment,aes(oddsratio,-log10(p),shape=What)) + 
  geom_point(color="gray",size=4) + 
  geom_point(data=subset(analysisEnrichment,p<0.05),aes(oddsratio,-log10(p),shape=What,color=Analysis),size=4) + 
  smartas_theme() + xlab("Odds ratio") + 
  scale_x_continuous(limits=c(0,2)) + 
  theme(legend.position="bottom")
ggsave("figures/featuresGainedAndLostVsRandom.png",p)

setwd(workingDir)

# 4.4 - Analyze Eduard Porta's interaction results ====
setwd("structural_analysis")

for (cancer in cancerTypes){

  intxFile <- paste0("EduardPorta-Interactions/report/",cancer,"_switches.txt")
  no_col <- max(count.fields(intxFile,sep = "\t"))
  interactions <- read.table(intxFile,header=F,fill=T,col.names=1:no_col)
  intxCols <- paste0("Interaction",1:(no_col-6))
  colnames(interactions) <- c("GeneId","Symbol","Normal_transcript","Tumor_transcript","nIso","tIso",intxCols)
  
  thisCancerCands <- candidatesDf[candidatesDf$Tumor==cancer,]
  thisI3d <- i3d[i3d$Cancer==cancer,]
  
  interactions <- merge(interactions,thisCancerCands)
  
  interactions$InteractionsAltered <- colSums(apply(interactions[,intxCols],1,function(x) {grepl("Lost.+",x) | grepl("Switched.+",x)} ))
  interactions$Driver[interactions$Driver==1] <- "Driver"
  interactions$Driver[interactions$Driver==0] <- "NotDriver"
  
  # plot num of patients
  p <- ggplot(interactions,aes(InteractionsAltered,NumPatients,color=Driver)) + geom_point()
  p <- p + theme_minimal()
  ggsave(paste0("figures/",cancer,".EP_InteractionsLost_vs_NumPatients.png"),p)
  
  # plot num of patients
  thisI3d.numIntx <- ddply(subset(thisI3d,InteractionAffection>10),.(Gene,Symbol,normalTranscript,tumorTranscript),summarise,nI3dBrokenIntxs=length(PartnerUniprot))
  interactions.thisI3d.numIntx <- merge(interactions,thisI3d.numIntx,by.x=c("GeneId","Symbol","Normal_transcript","Tumor_transcript"),by.y=c("Gene","Symbol","normalTranscript","tumorTranscript"))  
  
  p <- ggplot(interactions.thisI3d.numIntx,aes(InteractionsAltered,nI3dBrokenIntxs,color=Driver)) + geom_point()
  p <- p + theme_minimal()
  ggsave(paste0("figures/",cancer,".EP_InteractionsLost_vs_I3DInteractionsLost.png"),p)
}

setwd(workingDir)

# 4.5 -  ====
setwd("structural_analysis")



setwd(workingDir)

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
#   cancer_rel <- read.delim(paste0(basePath,cancer,".protein_centrality_relevant.tsv"), header=FALSE)
#   cancer_rel <- cbind(cancer_rel,"Relevant",cancer)
#   colnames(cancer_rel) <- c("Cancer","Degree","Relevance")
#   centrality_rel <- rbind(centrality_rel,cancer_rel)
#   
#   cancer_nrel <- read.delim(paste0(basePath,cancer,".protein_centrality_nonrelevant.tsv"), header=FALSE)
#   cancer_nrel <- cbind(cancer_nrel,"NonRelevant",cancer)
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
# Feats_affected_rel <- read.delim("structural_summary_relevant.tsv", header=FALSE)
# Feats_affected_rel <- cbind(Feats_affected_rel,"Relevant")
# colnames(Feats_affected_rel) <- c("Cancer","Switch","Pfam","PRINTS","ProSitePatterns","IUPREDLong","IUPREDShort","I3D","Driver","Relevance")
# 
# Feats_affected_nrel <- read.delim("structural_summary_nonrelevant.tsv", header=FALSE)
# Feats_affected_nrel <- cbind(Feats_affected_nrel,"NonRelevant")
# colnames(Feats_affected_nrel) <- c("Cancer","Switch","Pfam","PRINTS","ProSitePatterns","IUPREDLong","IUPREDShort","I3D","Driver","Relevance")
# 
# Feats_affected <- rbind(Feats_affected_rel,Feats_affected_nrel)
# Feats_affected$Driver <- ifelse(Feats_affected$Driver=="True","Driver","NonDriver")
# 
# png("ProSitePatterns.png", width=1000, height=800)
# par(mar=c(9,3,5,2))
# boxplot(ProSitePatterns~interaction(Driver,Cancer),data=Feats_affected[Feats_affected$Relevance=="Relevant" & Feats_affected$ProSitePatterns!=0,],col=c("gray40","gray90"),las=3,outline=F)
# graphics.off()
# 
# png("Pfam.png", width=1000, height=800)
# par(mar=c(9,3,5,2))
# boxplot(Pfam~interaction(Driver,Cancer),data=Feats_affected[Feats_affected$Relevance=="Relevant" & Feats_affected$Pfam!=0,],col=c("gray40","gray90"),las=3,outline=F)
# graphics.off()
# 
# png("PRINTS.png", width=1000, height=800)
# par(mar=c(9,3,5,2))
# boxplot(PRINTS~interaction(Driver,Cancer),data=Feats_affected[Feats_affected$Relevance=="Relevant" & Feats_affected$PRINTS!=0,],col=c("gray40","gray90"),las=3,outline=F)
# graphics.off()
# 
# png("I3D.png", width=1000, height=800)
# par(mar=c(9,3,5,2))
# boxplot(I3D~interaction(Driver,Cancer),data=Feats_affected[Feats_affected$Relevance=="Relevant" & Feats_affected$I3D!=0,],col=c("gray40","gray90"),las=3,outline=F)
# graphics.off()
# 
# png("IUPREDShort.png", width=1000, height=800)
# par(mar=c(9,3,5,2))
# boxplot(IUPREDShort~interaction(Driver,Cancer),data=Feats_affected[Feats_affected$Relevance=="Relevant" & Feats_affected$IUPREDShort!=0,],col=c("gray40","gray90"),las=3,outline=F)
# graphics.off()
# 
# png("IUPREDLong.png", width=1000, height=800)
# par(mar=c(9,3,5,2))
# boxplot(IUPREDLong~interaction(Driver,Cancer),data=Feats_affected[Feats_affected$Relevance=="Relevant" & Feats_affected$IUPREDLong!=0,],col=c("gray40","gray90"),las=3,outline=F)
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