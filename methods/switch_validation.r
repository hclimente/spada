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

withHclust <- function(x){

	withSwitch <- strsplit(as.character(x[4]),",")[[1]]
	withoutSwitch <- inputData$Replicates[!inputData$Replicates %in% withSwitch]
   
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
	return(c(precision,sensitivity))
}

withKmeans <- function(x){
  
  withSwitch <- strsplit(as.character(x[4]),",")[[1]]
  withoutSwitch <- inputData$Replicates[!inputData$Replicates %in% withSwitch]
   
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
       
  },error = function(e){},
  finally={
    if (!exists("precision")){
      precision <- NA
      sensitivity <- NA
    }
  })
  
  return(c(precision,sensitivity))
}

getSensitivityAndPrecision <- function(x){
  switches <- x[c("Gene","maxdPSI","mindPSI","Patients","pval")]

  kmeansV = withKmeans(switches)
  hclustV = withHclust(switches)
  
  precisionKmeans = kmeansV[1]
  sensitivityKmeans = kmeansV[2]
  precisionHclust = hclustV[1]
  sensitivityHclust = hclustV[2]

  return(data.frame(precisionKmeans=precisionKmeans,sensitivityKmeans=sensitivityKmeans,
                    precisionHclust=precisionHclust,sensitivityHclust=sensitivityHclust))
}

suppressMessages( library(logging) )

args <- commandArgs(trailingOnly = TRUE)
load(paste0(args[1],"RWorkspaces/2_GetCandidates.RData"))

logger <- getLogger(name="filter_switches", level=10) #Level debug

addHandler(writeToConsole, logger="filter_switches", level='INFO')
addHandler(writeToFile, logger="filter_switches", file=paste0(out, "rSmartAS.log"), level='DEBUG')

transcripts = data.frame(Gene=as.character(intraReplicate[[inputData$Replicates[1]]]$Gene),Tx=intraReplicate[[inputData$Replicates[1]]]$Transcript)
  
psiList <- list()
for (patient in inputData$Replicates){
  psiList[[patient]] <- intraReplicate[[patient]]$deltaPSI
}  
psis <- do.call("cbind",psiList)
deltaPsis <- cbind(transcripts,psis)

clustVals = apply(candidateList,1,getSensitivityAndPrecision)
clustVals = do.call('rbind',clustVals)

candidateList$precision = apply(cbind(clustVals$precisionKmeans,clustVals$precisionHclust),1,mean,na.rm=T)
candidateList$sensitivity = apply(cbind(clustVals$sensitivityKmeans,clustVals$sensitivityHclust),1,mean,na.rm=T)

candidateList$precisionKmeans = clustVals$precisionKmeans
candidateList$precisionHclust = clustVals$precisionHclust
candidateList$sensitivityKmeans = clustVals$sensitivityKmeans
candidateList$sensitivityHclust = clustVals$sensitivityHclust

write.table(candidateList, file=paste0(out, "candidateList_v2.tsv"), sep="\t", row.names=F, col.names=F, quote=F)

save(isoformExpression, intraReplicate, interReplicate, candidates, candidateList, inputData, allExpressedTranscripts, wd, out, file=paste0(out, "RWorkspaces/3_ClusteringFilter.RData"))