#!/soft/R/R-3.2.1/bin/Rscript

tryCatch(source("~/smartas/pipeline/scripts/variablesAndFunctions.r"),error=function(e){})

setwd(workingDir)
setwd("structural_analysis")

# 4.1 - i3d plots ====
# i3d <- read.delim("i3d_broken.tsv", header=TRUE,row.names=NULL)
# 
# i3d$Structure <- NULL
# i3d$Structure[grepl('-MDL-',i3d$PDB)] <- "MODEL"
# i3d$Structure[grepl('-EXP-',i3d$PDB)] <- "EXPERIMENTAL"
# i3d$Structure[grepl('-MDD-',i3d$PDB)] <- "DOMAIN_MODEL"
# 
# i3d$DriverRelationship <- "No driver involved"
# i3d$DriverRelationship[grep("Driver",i3d$Annotation)] <- "Driver"
# i3d$DriverRelationship[grep("Driver",i3d$PartnerAnnotation)] <- "d1"
# i3d$DriverRelationship[Reduce(intersect,  list(grep("Driver",i3d$Annotation),grep("Driver",i3d$PartnerAnnotation)))] <- "Driver-Driver"
# 
# p <- ggplot(i3d,aes(InteractionAffection,(SequenceCover+PartnerCover)/2)) + geom_point(aes(shape=DriverRelationship))
# p <- p + geom_point(data=subset(i3d,InteractionAffection > 60 & (SequenceCover+PartnerCover)/2 > 60),aes(InteractionAffection,(SequenceCover+PartnerCover)/2,color=Symbol,shape=DriverRelationship ))
# p <- p + theme_minimal()
# p <- direct.label(p)
# ggsave("figures/i3d_InteractionVsAffection.png",p)
# 
# # take only interactions where both partners are sufficiently
# # described and the interaction is actually affected
# i3d_complete <- i3d[i3d$SequenceCover>80 & i3d$InteractionAffection>0 & i3d$PartnerCover>80,!colnames(i3d) %in% c("SequenceInformation","IsoformSpecific","PDB")]
# i3d_complete <- i3d_complete[order(-i3d_complete$InteractionAffection),]
# 
# # unique interactions
# length(unique(i3d$PDB[i3d$SequenceCover>80 & i3d$InteractionAffection>0 & i3d$PartnerCover>80]))
# table(unique(i3d[i3d$SequenceCover>80 & i3d$InteractionAffection>0 & i3d$PartnerCover>80,c("PDB","Structure","InteractionAffection")])$Structure)
# 
# df <- unique(i3d_complete[,c("Symbol","Partner","normalTranscript","tumorTranscript","InteractionAffection","Annotation","PartnerAnnotation")])
# nrow(df)
# 
# p <- ggplot(df, aes(InteractionAffection)) + geom_histogram(binwidth=5,fill="#c0392b", alpha=0.75)
# p <- p + smartas_theme() + labs(x="Percentage of the interface lost", y="Frequency")
# p <- p + geom_hline(yintercept=0, size=0.4, color="black") + scale_x_continuous(breaks=seq(0,100, by=5))
# ggsave("figures/i3d_interactionAffected_completeStructures.png",p)
# 
# freqIntx = as.data.frame(table(df[,c("Annotation","PartnerAnnotation")]))
# freqIntx <- freqIntx[order(-freqIntx$Freq),]
# write.table(freqIntx,'tables/i3d_AnnotationFrequency.txt',quote=F,col.names=T,sep="\t",row.names=FALSE)
# 
# nrow(i3d_complete)
# table(i3d_complete$WhatsHappening)
# 
# scoredTable <- merge(i3d,gn_merged[,colnames(gn_merged) %in% c("Gene","Score")],by=c("Gene"))
# scoredTable <- rename(scoredTable,c("Score"="GeneScore"))
# scoredTable <- merge(scoredTable,gnset_merged[,colnames(gnset_merged) %in% c("Gene","Score")],by=c("Gene"))
# scoredTable <- rename(scoredTable,c("Score"="GenesetScore"))
# 
# filtScoredTable <- unique(scoredTable[scoredTable$SequenceCover>50 & scoredTable$PartnerCover>50 & scoredTable$InteractionAffection > 10,c("Symbol","Partner","InteractionAffection","Annotation","PartnerAnnotation","SequenceCover","PartnerCover","WhatsHappening","GeneScore","GenesetScore","Structure")])
# filtScoredTable <- filtScoredTable[order(-filtScoredTable$GeneScore),]
# 
# filtScoredTable <- filtScoredTable[order(-filtScoredTable$GenesetScore),]

# 4.2 - domain enrichment in drivers,oncogenes,suppressors ====
structural_features <- read.delim("structural_features.onlyModels.tsv")
nrStructural_features <- structural_features[structural_features$Random == "NonRandom",]

for (featureType in c('Pfam','prosite')){
  
  annotation <- read.delim(paste0(workingDir,"/data/Databases/",featureType,"2go.clean.txt"),row.names = NULL,header=F)
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

# 4.3 - Study the enrichment in a particular set of features being gained or lost ====
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

# 4.4 - Analyze Eduard Porta's interaction results ====
# for (cancer in cancerTypes){
#   
#   intxFile <- paste0("EduardPorta-Interactions/report/",cancer,"_switches.txt")
#   no_col <- max(count.fields(intxFile,sep = "\t"))
#   interactions <- read.table(intxFile,header=F,fill=T,col.names=1:no_col)
#   intxCols <- paste0("Interaction",1:(no_col-6))
#   colnames(interactions) <- c("GeneId","Symbol","Normal_transcript","Tumor_transcript","nIso","tIso",intxCols)
#   
#   thisCancerCands <- candidatesDf[candidatesDf$Tumor==cancer,]
#   thisI3d <- i3d[i3d$Cancer==cancer,]
#   
#   interactions <- merge(interactions,thisCancerCands)
#   
#   interactions$InteractionsAltered <- colSums(apply(interactions[,intxCols],1,function(x) {grepl("Lost.+",x) | grepl("Switched.+",x)} ))
#   interactions$Driver[interactions$Driver==1] <- "Driver"
#   interactions$Driver[interactions$Driver==0] <- "NotDriver"
#   
#   # plot num of patients
#   p <- ggplot(interactions,aes(InteractionsAltered,NumPatients,color=Driver)) + geom_point()
#   p <- p + theme_minimal()
#   ggsave(paste0("figures/",cancer,".EP_InteractionsLost_vs_NumPatients.png"),p)
#   
#   # plot num of patients
#   thisI3d.numIntx <- ddply(subset(thisI3d,InteractionAffection>10),.(Gene,Symbol,normalTranscript,tumorTranscript),summarise,nI3dBrokenIntxs=length(PartnerUniprot))
#   interactions.thisI3d.numIntx <- merge(interactions,thisI3d.numIntx,by.x=c("GeneId","Symbol","Normal_transcript","Tumor_transcript"),by.y=c("Gene","Symbol","normalTranscript","tumorTranscript"))  
#   
#   p <- ggplot(interactions.thisI3d.numIntx,aes(InteractionsAltered,nI3dBrokenIntxs,color=Driver)) + geom_point()
#   p <- p + theme_minimal()
#   ggsave(paste0("figures/",cancer,".EP_InteractionsLost_vs_I3DInteractionsLost.png"),p)
# }