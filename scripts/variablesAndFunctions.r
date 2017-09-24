#!/soft/R/R-3.2.1/bin/Rscript

# 0 - ENVIRONMENT ---------------------
library(plyr)
library(ggplot2)
library(reshape2)
#library(directlabels)
library(gridExtra)
library(magrittr)
library(readr)
library(dplyr)
#library(wordcloud)

cancerTypes <- c("brca","coad","hnsc","kich","kirc","kirp","lihc","luad","lusc","prad","thca")
workingDir <- "/genomics/users/hector/smartas/results"
project <- "/projects_rg/TCGA/users/hector/SmartAS/"

colorPalette <- c("#CC79A7","#636363","#F0E442","#006D2C","#31A354","#74C476","#FC8D59","#08519C","#3182BD","#D55E00","#5E3C99","#000000","#969696","#EDF8FB","#B3CDE3","#8C96C6","#88419D","#D94701","#2171B5","#FD8D3C","#6BAED6","#FD8D3C","#33A02C","#E31A1C","#FF7F00","#6A3D9A","#B15928","#377EB8","#E41A1C")
names(colorPalette) <- c("brca","coad","hnsc","kich","kirc","kirp","lihc","luad","lusc","prad","thca","coad-hyper","coad-hypo","brca-basal","brca-her2","brca-luminala","brca-luminalb","expression-up","expression-down","psi-up","psi-down","odds","a3","a5","mx","ri","se","normal","tumor")

nPatients <- list()
nPatients[["brca"]] <- 1036
nPatients[["coad"]] <- 262
nPatients[["hnsc"]] <- 422
nPatients[["kich"]] <- 62
nPatients[["kirc"]] <- 505
nPatients[["kirp"]] <- 195
nPatients[["lihc"]] <- 197
nPatients[["luad"]] <- 488
nPatients[["lusc"]] <- 483
nPatients[["prad"]] <- 295
nPatients[["thca"]] <- 497
nPatients[["total"]] <- sum(1036,262,422,62,505,195,197,488,483,295,497)

nPatientsDf <- as.data.frame(do.call("rbind",nPatients))
nPatientsDf$Cancer <- rownames(nPatientsDf)
colnames(nPatientsDf) <- c("TotalPatients","Cancer")

# 0.1 - plotting functions ====
getBarplotAsterisks <- function(stat.tests,ranges,categories=cancerTypes,barsize=2, pairAxis=0.5){

  minStep <- max(ranges$y)/40

  r <- 0.5
  t <- seq(0, 180, by = 1) * pi / 180
  x <- r * cos(t)
  y <- 1 * sin(t)
  y[which(y>minStep)] <- minStep
  modelArc <- data.frame(x = x, y = y)

  arcs <- data.frame()
  ast <- data.frame()
  for (kns in stat.tests$Cancer[stat.tests$p < 0.05]){
    i <- which(categories==kns)
    d <- ddply(ranges[ranges$Cancer==kns,],.(variable),summarise,max=max(y))
    j <- max(d$max)
    k <- min(d$max)
    whichBar <- which(levels(ranges$variable)==d$variable[d$max==k]) - barsize
    tmp_arc <- data.frame(Cancer=i, x_arc=modelArc$x+(barsize*i-pairAxis), y_arc=modelArc$y+j)
    #     tmp_arc <- rbind(tmp_arc,
    #       data.frame(Cancer=i, x_arc=2*i+whichBar, y_arc=k))
    if (which(d$max==k)==1){
      tmp_arc <- rbind(tmp_arc,data.frame(Cancer=i, x_arc=barsize*i-1, y_arc=k))
      arcs <- rbind(tmp_arc,arcs)
    } else if (which(d$max==k)==2){
      tmp_arc <- rbind(tmp_arc,data.frame(Cancer=i, x_arc=barsize*i, y_arc=k))
      arcs <- rbind(arcs,tmp_arc)
    }

    tmp_ast <- data.frame(Cancer=i,x_ast=barsize*i-pairAxis ,y_ast=j + 1.2*minStep)
    ast <- rbind(ast,tmp_ast)
  }

  arcs <- arcs[order(-arcs$x_arc),]

  significanceElements <- list("asterisks" = ast, "arcs" = arcs)

  return(significanceElements)

}

getBoxplotAsterisks <- function(stat.tests,ranges,rm.outliers=TRUE,categories=cancerTypes){

  if (rm.outliers){
    plotRanges <- ddply(ranges, .(Cancer,Categories), summarise, Q1=quantile(y, 1/4), Q3=quantile(y, 3/4), IQR=Q3-Q1, upper.limit=Q3+1.5*IQR)
    plotRanges <- ddply(plotRanges, .(Cancer), summarise, y=max(upper.limit))
  } else {
    plotRanges <- ddply(ranges, .(Cancer), summarise, y=max(y))
  }

  minStep <- max(plotRanges$y)/40

  r <- 0.25
  t <- seq(0, 180, by = 1) * pi / 180
  x <- r * cos(t)
  y <- minStep * sin(t)
  y[which(y>minStep/2)] <- minStep/2
  modelArc <- data.frame(x = x, y = y)

  arcs <- data.frame(Cancer=character(),x_arc=numeric() ,y_arc=numeric())
  ast <- data.frame(Cancer=character(),x_ast=numeric() ,y_ast=numeric())
  for (kns in stat.tests$Cancer[stat.tests$p < 0.05]){
    i <- which(categories==kns)
    j <- max(plotRanges$y[plotRanges$Cancer==kns])
    tmp_arc <- data.frame(Cancer=i, x_arc=modelArc$x+(i), y_arc=modelArc$y+j)
    arcs <- rbind(arcs,tmp_arc)

    tmp_ast <- data.frame(Cancer=i,x_ast=i ,y_ast=j + 1.2*minStep)
    ast <- rbind(ast,tmp_ast)
  }

  significanceElements <- list("asterisks" = ast, "arcs" = arcs)

  return(significanceElements)

}

getAnnotationGroup <- function(x){
  if (length(intersect(x,topGroups)) > 0){
    for (y in topGroups){
      if ( y %in% x){
        return(y)
      }
    }
  } else {
    return("Other")
  }
}

splitTextInLines <- function(x){
  y <- gsub("_"," ",as.character(x))
  splittedFeatures <- strsplit(y," ")
  formattedFeatures <- character()
  for (i in splittedFeatures){
    feat <- i[1]
    for (j in 2:length(i)){
      if (j%%5 == 0){
        feat <- paste(feat,i[j],sep="\n")
      } else {
        feat <- paste(feat,i[j],sep=" ")
      }
    }
    formattedFeatures <- c(formattedFeatures,feat)
  }
  return(formattedFeatures)
}

studyGroups <- function(x,y,switchesTable){

  fisherTest <- function(cTable,analysisName){

    if (all(dim(cTable) == c(2,2))){
      f <- fisher.test(cTable)
      #print(cTable)

      if (f$p.value > 0.05){
        groupSignif <- "Non significant"
      } else {
        groupSignif <- "Significant"
      }

      if (f$estimate > 1){
        groupEnrich <- paste0(groupSignif," enrichment of ",rownames(cTable)[1]," in Group 1 and of ",rownames(cTable)[2]," in group 2")
      } else if (f$estimate < 1){
        groupEnrich <- paste0(groupSignif," enrichment of ",rownames(cTable)[2]," in Group 1 and of ",rownames(cTable)[1]," in group 2")
      } else {
        groupEnrich <- "No difference between the groups"
      }

      cat(paste0("*",analysisName," enrichment analysis\n",groupEnrich," (p = ",round(f$p.value,2),"; odds ratio = ",round(f$estimate,2),")\n"))
    } else {
      cat(paste0("*",analysisName," enrichment analysis\ncannot be performed because there are not enough different features.\n"))
    }
  }

  tTest <- function(x,y,analysisName){

    t <- t.test(x,y)

    if (t$p.value > 0.05){
      groupSignif <- "Non significant difference."
    } else {
      groupSignif <- "Significant difference."
    }

    if (f$estimate > 1){
      groupEnrich <- paste0(groupSignif,"Group 1 higher than group 2")
    } else if (f$estimate < 1){
      groupEnrich <- paste0(groupSignif,"Group 2 higher than group 1")
    } else {
      groupEnrich <- "No difference between the groups"
    }

    cat(paste0("*",analysisName," mean difference test\n",groupEnrich," (p = ",round(t$p.value,2),", " , round(mean(x),2)," vs. ", round(mean(y),2),")\n"))
  }

  wilcoxTest <- function(x,y,analysisName){
    tryCatch({
      w <- wilcox.test(x,y)

      if (w$p.value > 0.05){
        groupSignif <- "Non significant difference. "
      } else {
        groupSignif <- "Significant difference. "
      }

      if (median(x) > median(y)){
        groupEnrich <- paste0(groupSignif,"Group 1 higher than group 2")
      } else if (median(x) < median(y)){
        groupEnrich <- paste0(groupSignif,"Group 2 higher than group 1")
      } else {
        groupEnrich <- "No difference between the groups"
      }

      cat(paste0("*",analysisName," median difference test\n",groupEnrich," (",round(median(x),2)," vs. ", round(median(y),2),", p = ",round(w$p.value,2),")\n"))
      }, error = function(err){
        if (err$message=="not enough 'y' observations" | err$message=="not enough 'x' observations"){
          cat(paste0("*",analysisName," median difference test could not be performed\n"))
          w <- "nope"
        } else {
          stop(err)
        }
      }
    )
  }

  cat("Size: group 1 - ",nrow(x),"; group 2 - ",nrow(y),"\n")

  set1 <- merge(x,switchesTable)
  set2 <- merge(y,switchesTable)

  set1$Group <- 1
  set2$Group <- 2

  sets <- rbind(set1,set2)

  # relevance
  cTable <- table(sets[,c("IsFunctional","Group")])
  fisherTest(cTable,"relevance")

  # CDS_change
  cTable <- table(sets[,c("CDS_change","Group")])
  fisherTest(cTable,"CDS change")

  # UTR_change
  cTable <- table(sets[,c("UTR_change","Group")])
  fisherTest(cTable,"UTR change")

  # number of patients
  if ("PatientNumber" %in% colnames(sets) ){
    wilcoxTest(sets$PatientNumber[sets$Group==1],sets$PatientNumber[sets$Group==2],"number of patients")
  }

  if ("NumPatients" %in% colnames(sets) ){
    wilcoxTest(sets$NumPatients[sets$Group==1],sets$NumPatients[sets$Group==2],"number of patients")
  }

  # driver
  cTable <- table(sets[,c("Driver","Group")])
  fisherTest(cTable,"drivers")

  # enrichment in driver types
  cTable <- table(sets[,c("DriverType","Group")])

  for (i in 2:nrow(cTable)){
    thisCTable <- cTable[c(1,i),]
    fisherTest(thisCTable,rownames(cTable)[i])
  }

  # enrichment in driver types
  cTable <- table(sets[,c("DriverAnnotation","Group")])

  for (i in 1:2){
    thisCTable <- cTable[c(3,i),]
    fisherTest(thisCTable,rownames(cTable)[i])
  }

  # hallmarks
  sets$Hallmark <- "No"
  sets$Hallmark[grepl("HALLMARK",sets$Annotation)] <- "Hallmark"

  cTable <- table(sets[,c("Hallmark","Group")])
  fisherTest(cTable,"hallmark")

  # cancer unbalances
  if ("CancerAffected" %in% colnames(sets) ){
    cancerList <- lapply(strsplit(as.character(sets$CancerAffected),","), function(x) { y <- cancerTypes; y[!(y %in% x)] <- "Other"; y })
    cancerList <- do.call("rbind",cancerList)

    # add number of tumors affected
    numTum <- apply(cancerList,1,function(x) { sum(x!="Other",na.rm=T) })
    sets$numberOfTumors <- numTum

    wilcoxTest(sets$numberOfTumors[sets$Group==1],sets$numberOfTumors[sets$Group==2],"number of tumors")
  }

  if ("InteractionsAltered" %in% colnames(sets) ){
    wilcoxTest(sets$InteractionsAltered[sets$Group==1 & sets$InteractionsAltered!=0],sets$InteractionsAltered[sets$Group==2 & sets$InteractionsAltered!=0],"interactions altered")
  }

  if ("InteractionsKept" %in% colnames(sets) ){
    wilcoxTest(sets$InteractionsKept[sets$Group==1 & sets$InteractionsKept!=0 ],sets$InteractionsKept[sets$Group==2 & sets$InteractionsKept!=0],"interactions kept")
  }

  if ("Score" %in% colnames(sets) ){
    wilcoxTest(sets$Score[sets$Group==1],sets$Score[sets$Group==2],"meScore")
  }

}

smartas_theme <- function() {

  library(RColorBrewer)

  # Generate the colors for the chart procedurally with RColorBrewer
  palette <- brewer.pal("Greys", n=9)
  color.background = "white"
  color.grid.major = palette[3]
  color.axis.text = palette[6]
  color.axis.title = palette[7]
  color.title = palette[9]

  # Begin construction of chart
  theme_bw(base_size=20) +

    # Set the entire chart region to a light gray color
    theme(panel.background=element_rect(fill=color.background, color=color.background)) +
    theme(plot.background=element_rect(fill=color.background, color=color.background)) +
    theme(panel.border=element_rect(color=color.background)) +

    # Format the grid
    theme(panel.grid.major=element_line(color=color.grid.major,size=.25)) +
    theme(panel.grid.minor=element_blank()) +
    theme(axis.ticks=element_blank()) +

    # Format the legend, but hide by default
    theme(legend.position="none") +
    theme(legend.background = element_rect(fill=color.background)) +
    theme(legend.text = element_text(size=20,color=color.axis.title)) +

    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(color=color.title, size=25, vjust=1.25)) +
    theme(axis.text.x=element_text(size=15,color=color.axis.text)) +
    theme(axis.text.y=element_text(size=15,color=color.axis.text)) +
    theme(axis.title.x=element_text(size=20,color=color.axis.title, vjust=0)) +
    theme(axis.title.y=element_text(size=20,color=color.axis.title, vjust=1.25)) +

    # Plot margins
    theme(plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm"))
}
