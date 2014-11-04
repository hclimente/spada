#!/soft/R/R-3.0.0/bin/Rscript

colLeafs <- function(n,withSwitch,withoutSwitch){
  if (is.leaf(n)){
    a <- attributes(n)
    if (a$label %in% withSwitch){
      attr(n, "nodePar") <- c(a$nodePar, lab.col = "#fe901d")
    } else if (a$label %in% withoutSwitch){
      attr(n, "nodePar") <- c(a$nodePar, lab.col = "#3762D2")
    }
  }
  n
}

countCols <- function(tree,count){
  count <- character()
  for (i in seq(1,length(tree))){
    if (!is.leaf(tree[[i]])){
      count <- c(count,countCols(tree[[i]],count))
    } else {
      count <- c(count,attributes(tree[[i]])$label)
    }
  }
  
  return(count)
  
}

createSubpopulations <- function(x){

	withSwitch <- strsplit(as.character(x[4]),",")[[1]]
	withoutSwitch <- inputData$Replicates[!inputData$Replicates %in% withSwitch]
  
  p <- -log10(as.numeric(x[5]))
  
	relevantDf <- deltaPsis[deltaPsis$Gene==x[1],-1]
	relevantDf <- as.data.frame(t(relevantDf))
	colnames(relevantDf) <- unlist(relevantDf[1,])
	relevantDf <- relevantDf[-1, ]
  

  
	# Ward Hierarchical Clustering
	d <- dist(relevantDf, method="euclidean") # distance matrix
	fit <- tryCatch({
	                  fit <- hclust(d, method="ward")
	                  fit <- as.dendrogram(fit)
	                  
	                  refit <- dendrapply(fit,colLeafs,withSwitch,withoutSwitch)
	                  oneBranch <- countCols(refit[[1]])
	                  otherBranch <- countCols(refit[[2]])
	                  
	                  oneBranchPositives = sum(withSwitch %in% oneBranch)/length(oneBranch)
	                  otherBranchPositives = sum(withSwitch %in% otherBranch)/length(otherBranch)
	                  
	                  if (oneBranchPositives > otherBranchPositives){
                      TP <- sum(withSwitch %in% oneBranch)
                      FP <- sum(withSwitch %in% otherBranch)
                      TN <- sum(withoutSwitch %in% otherBranch)
	                    FN <- sum(withoutSwitch %in% oneBranch)
	                    
	                  } else {
	                    TP <- sum(withSwitch %in% otherBranch)
	                    FP <- sum(withSwitch %in% oneBranch)
	                    TN <- sum(withoutSwitch %in% otherBranch)
	                    FN <- sum(withoutSwitch %in% oneBranch)
	                  }
	                  precision = TP/(TP+FP)
	                  sensitivity = TP/(TP+FN)
	                  
	                  #plot(refit,main=paste(x[1],x[2],x[3],sep="_")) # display dendogram
	                  #cat(paste(x[1],x[2],x[3],sep="_"),"\tSpecificity:",precision,"\tsensitivity:",sensitivity,"\n") # display dendogram
                  },error = function(e){},
                  finally={
                    if (!exists("precision")){
                      precision <- NA
                      sensitivity <- NA
                    }
                  })
	return(c(precision,sensitivity,p))
}

createSubpopulationsKmeans <- function(x){
  
  withSwitch <- strsplit(as.character(x[4]),",")[[1]]
  withoutSwitch <- inputData$Replicates[!inputData$Replicates %in% withSwitch]
  
  p <- -log10(as.numeric(x[5]))
  
  relevantDf <- deltaPsis[deltaPsis$Gene==x[1],-1]
  relevantDf <- as.data.frame(t(relevantDf))
  colnames(relevantDf) <- unlist(relevantDf[1,])
  relevantDf <- relevantDf[-1, ]
  
  psiData <- as.matrix(relevantDf)
  class(psiData) <- "numeric"
  psiData <- na.omit(psiData)
  #psiData <- scale(psiData)
  
  fit <- tryCatch({
    fit <- kmeans(psiData,2)
    aggregate(psiData,by=list(fit$cluster),FUN=mean)
    
    mydata <- data.frame(psiData, fit$cluster) 
    
    oneBranch <- (row.names(mydata))[mydata$fit.cluster == "1"]
    otherBranch <- (row.names(mydata))[mydata$fit.cluster == "2"]
    
    oneBranchPositives = sum(withSwitch %in% oneBranch)/length(oneBranch)
    otherBranchPositives = sum(withSwitch %in% otherBranch)/length(otherBranch)
    
    if (oneBranchPositives > otherBranchPositives){
      TP <- sum(withSwitch %in% oneBranch)
      FP <- sum(withSwitch %in% otherBranch)
      TN <- sum(withoutSwitch %in% otherBranch)
      FN <- sum(withoutSwitch %in% oneBranch)
      
    } else {
      TP <- sum(withSwitch %in% otherBranch)
      FP <- sum(withSwitch %in% oneBranch)
      TN <- sum(withoutSwitch %in% otherBranch)
      FN <- sum(withoutSwitch %in% oneBranch)
    }
    precision = TP/(TP+FP)
    sensitivity = TP/(TP+FN)
    
    #library(cluster)
    #clusplot(mydata, fit$cluster, color=TRUE, shade=TRUE,labels=2, lines=0)
    
  },error = function(e){},
  finally={
    if (!exists("precision")){
      precision <- NA
      sensitivity <- NA
    }
  })
  

  return(c(precision,sensitivity,p))
}

getDeltaPSI <- function(){
  genes = data.frame(Gene=as.character(intraReplicate[[inputData$Replicates[1]]]$Gene),Tx=intraReplicate[[inputData$Replicates[1]]]$Transcript)
  
  psiList <- list()
  for (patient in inputData$Replicates){
    psiList[[patient]] <- intraReplicate[[patient]]$deltaPSI
    #thisDeltaPSI <- intraReplicate[[patient]][,c("Gene","deltaPSI")]
    #result <- merge(result,thisDeltaPSI,by="Gene",suffixes = c("",paste0("_",patient)))
  }  
  psis <- do.call("cbind",psiList)
  result <- cbind(genes,psis)
  return(result)
}

args <- commandArgs(trailingOnly = TRUE)
load(paste0(args[1],"RWorkspaces/2_GetCandidates.RData"))
tag <- args[2]

deltaPsis <- getDeltaPSI()

precisionKmeans <- numeric()
sensitivityKmeans <- numeric()
pvalKmeans <- numeric()

precisionHclust <- numeric()
sensitivityHclust <- numeric()
pvalHclust <- numeric()

for (gene in unique(candidateList[candidateList$pval <= 0.001,"Gene"])){
	switches <- candidateList[candidateList$Gene==gene & candidateList$pval <= 0.001,c("Gene","maxdPSI","mindPSI","Patients","pval")]
  
	kmeans <- apply(switches,1,createSubpopulations)
  
	precisionKmeans <- c(precisionKmeans,kmeans[1])
	sensitivityKmeans <- c(sensitivityKmeans,kmeans[2])
	pvalKmeans <- c(pvalKmeans,kmeans[3])
  
	Hclust <- apply(switches,1,createSubpopulationsKmeans)
	
  precisionHclust <- c(precisionHclust,Hclust[1])
	sensitivityHclust <- c(sensitivityHclust,Hclust[2])
	pvalHclust <- c(pvalHclust,Hclust[3])
}

corrPPKmeans <- cor(precisionKmeans,pvalKmeans,use="complete.obs",method="spearman")
corrSPKmeans <- cor(sensitivityKmeans,pvalKmeans,use="complete.obs",method="spearman")
meanPKmeans <- mean(precisionKmeans,na.rm=T)
meanSKmeans <- mean(sensitivityKmeans,na.rm=T)

corrPPHclust <- cor(precisionHclust,pvalHclust,use="complete.obs",method="spearman")
corrSPHclust <- cor(sensitivityHclust,pvalHclust,use="complete.obs",method="spearman")
meanPHclust <- mean(precisionHclust,na.rm=T)
meanSHclust <- mean(sensitivityHclust,na.rm=T)

#plot(precision,pval)
#plot(sensitivity,pval)

df <- data.frame(cancer=c(tag,tag),precision=c(meanPKmeans,meanPHclust),sensitivity=c(meanSKmeans,meanSHclust),precision_pval=c(corrPPKmeans,corrPPHclust),sensitivity_pval=c(corrSPKmeans,corrSPHclust))

rownames(df) <- c("kmeans","hclust")

write.table(df, paste0(args[1],"/result_summary/accuracy.tsv"), sep="\t", quote=F)