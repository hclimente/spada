library(gplots)

all_drivers <- read.delim("~/Downloads/all_drivers.tsv")
dgidb_export_all_drivers_bygene_results <- read.delim("~/Downloads/dgidb_export_all_drivers_bygene_results.tsv", dec=",")

kk <- NULL
patientOrigin = data.frame(Patient=character(),Cancer=character())
for (folder in list.files(path="~/SmartAS/testResults/TCGA")){
  if (folder=="analysis"){next}
  
  relevantSwitches <- read.delim(paste0("~/SmartAS/testResults/TCGA/",folder,"/result_summary/relevantSwitches.tsv"))
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

kk <- kk[rowSums(kk) > 10,]

cancerColor<-read.delim("~/Downloads/known_cancer_colour")  #I am reading here the seperate file that I sent you for cancer colours.
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

geneColors <- rep('white',nrow(kk))

drivers <- rownames(kk) %in% all_drivers$Input.symbol | rownames(kk) %in% all_drivers$Approves.symbol
druggable <- rownames(kk) %in% dgidb_export_all_drivers_bygene_results$Search.Term

geneColors[drivers] <- "red"
geneColors[druggable] <- "black"

png("~/Downloads/relevantGenesPerPatient.png", width=1500, height=1500)
heatmap.2(as.matrix(kk),col=colorpanel(2,"white", "blue"), trace="none", Rowv=FALSE,Colv=FALSE,na.rm=TRUE,
          labRow = rownames(kk), labCol = FALSE,
          xlab = NULL, ylab = NULL,dendrogram='none',margins = c(5, 13),
          ColSideColors = color, RowSideColors=geneColors, cexRow=2,key = FALSE,
          main= "Relevant genes affection")
legend(x="topleft",legend=cancerColor$cancer,fill=cancerColor$color, border=FALSE)
graphics.off()

png("~/Downloads/relevantGenesPerPatient_wDendrogram.png", width=1500, height=1500)
heatmap.2(as.matrix(kk),col=colorpanel(2,"white", "blue"), trace="none", na.rm=TRUE,
          labRow = rownames(kk), labCol = FALSE,
          distfun = function(x) dist(x, method = 'euclidean'),
          hclustfun = function(x) hclust(x, method = 'ward'),
          xlab = NULL, ylab = NULL,dendrogram='both',
          ColSideColors = color,RowSideColors=geneColors, cexRow=2,margins = c(5, 13), 
          key = FALSE,main= "Relevant genes affection")
legend(x="topleft",legend=c("Driver","Druggable"),fill=c("red","black"), border=FALSE)
graphics.off()