library(gplots)

kk <- NULL
patientOrigin = data.frame(Patient=character(),Cancer=character())
for (folder in list.files(path="~/SmartAS/Results/TCGA")){
  if (folder=="analysis"){next}
  
  relevantSwitches <- read.delim(paste0("~/SmartAS/Results/TCGA/",folder,"/result_summary/relevantSwitches.tsv"))
  if (is.null(kk)){
    kk <- relevantSwitches
  } else {
    kk <- merge(kk,relevantSwitches,all=TRUE) 
  }
  pp <- cbind(colnames(relevantSwitches)[2:ncol(relevantSwitches)],unlist(strsplit(folder,"_"))[1])
  patientOrigin <- rbind(patientOrigin,pp)
  
}

kk[is.na(kk)] <- 0
rownames(kk) <- kk$Genes
kk <- kk[,-1]

cancerColor<-read.delim("known_cancer_colour")  #I am reading here the seperate file that I sent you for cancer colours.
cancerColor$color<-as.character(cancerColor$color)

cancer = unlist(lapply(strsplit(colnames(kk),"\\."),"[[",1))

color <- character()
for (can in cancer){
  color <- c(color,switch(can, "brca"={cancerColor[1,2]},"coad"={cancerColor[2,2]},
                  "hnsc"={cancerColor[3,2]},"kich"={cancerColor[4,2]},
                  "kirc"={cancerColor[5,2]},"kirp"={cancerColor[6,2]},
                  "lihc"={cancerColor[7,2]},"luad"={cancerColor[8,2]},            
                  "lusc"={cancerColor[9,2]},"prad"={cancerColor[10,2]},
                  "thca"={cancerColor[11,2]}))
  
}

kk <- as.matrix(kk)

heatmap.2(as.matrix(kk),col=colorpanel(2,"white", "blue"), trace="none", Rowv=FALSE,Colv=FALSE,na.rm=TRUE,
          labRow = rownames(kk), labCol = colnames(kk),
          xlab = NULL, ylab = NULL,dendrogram='none',
          ColSideColors = color,cexRow=0.6, 
          key = FALSE,
          main= "Relevant genes affection")

legend(x="bottomleft",legend=cancerColor$cancer,fill=cancerColor$color, cex=1.3, border=FALSE)