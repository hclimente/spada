#!/usr/bin/env Rscript

tryCatch(source("~/smartas/pipeline/scripts/variablesAndFunctions.r"),error=function(e){})

setwd(workingDir)
setwd('neighborhood_analysis')

# 2 - NEIGHBORHOODS ---------------------
# 2.1 - Count number of times a geneset is significantly altered in different cancer types ====
genesetType <- c("canonical_pathways","hallmarks","go_biological_process","oncogenic_signatures")

for (gnset in c("functional","all")){
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
    
    if (sum(set_counts_df$Counts>8)){
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
}

rm(p,genesetType,gnset,type,file,sets_raw,set_counts,set_counts_df)