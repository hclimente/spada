#!/usr/bin/env Rscript

tryCatch(source("~/smartas/pipeline/scripts/variablesAndFunctions.r"),error=function(e){})

setwd(workingDir)
setwd('mutations')

# 3 - MUTATION OVERLAP ---------------------

# 3.1 - calculate allCancers p.values for each combination of mutation and switch types ====
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

# 3.2 - calculate gene meScore ====
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

# 3.3 - calculate geneset meScore ====
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

# 3.4 - calculate p.value of overlap between mutation in a feature and switches affecting it ====
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

# 3.5 - Plot p_me in a gene against the p_overlap of mutations and features ====
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

# 3.6 - feature overlap ====
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

# 3.7 - Search for something that makes the meScore+ special :)  ====
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