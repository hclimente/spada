library(ComplexHeatmap)
library(GetoptLong)

# arguments
setwd("/genomics/users/hector/TCGA_analysis/hallmark_info/")

# args <- commandArgs(trailingOnly = TRUE)

# geneset_path <- args[1]
geneset_path <- "SWI_SNF_COMPEX_functional_mutations.tsv"
hallmark_filename <- unlist(strsplit(x = geneset_path,split = "\\."))[1]
hallmark_name <- gsub("_", " ", hallmark_filename)
hallmark_name <- paste("OncoPrint",hallmark_name)

filter_out_unfrequent <- TRUE

# get matrix
unmeltHallmarks <- function(path){
  patient_info <- read.delim(path)
  patients <- unique(patient_info$Patient)
  genes <- unique(patient_info$Gene)
  unmelt_table <- as.data.frame(matrix(data="",nrow=length(patients),ncol=length(genes)),stringsAsFactors=FALSE)
  unmelt_table_numeric <- as.data.frame(matrix(data=0,nrow=length(patients),ncol=length(genes)),stringsAsFactors=FALSE)
  patientCancerDf <- unique(patient_info[,c("Patient","Cancer")])
  
  patientCancerCorrespondance <- as.character(patientCancerDf$Cancer)
  names(patientCancerCorrespondance) <- as.character(patientCancerDf$Patient)
  
  colnames(unmelt_table) <- genes
  rownames(unmelt_table) <- patients
  
  colnames(unmelt_table_numeric) <- genes
  rownames(unmelt_table_numeric) <- patients
  
  for (i in 1:nrow(patient_info)){
    thisPatient <- patient_info$Patient[i]
    thisGene <- patient_info$Gene[i]
  
    selectCol <- colnames(unmelt_table) == thisGene
    selectRow <- rownames(unmelt_table) == thisPatient
    
    thisAlteration <- paste0(patient_info$State[i],";")
    currentValue <- unmelt_table[selectRow,selectCol]
    thisAlteration <- paste0(currentValue,thisAlteration)
    
    unmelt_table[selectRow,selectCol] <- as.character(thisAlteration)
    if (patient_info$State[i]=="SWITCH"){
      alterationScore <- 1
    } else if (patient_info$State[i]=="MUT"){
      alterationScore <- 1
    }
    unmelt_table_numeric[selectRow,selectCol] <- unmelt_table_numeric[selectRow,selectCol] + alterationScore
      
  }
  unmelt_table <- t(unmelt_table)
  unmelt_table_numeric <- t(unmelt_table_numeric)
  return(list(unmelt_table,unmelt_table_numeric,patientCancerCorrespondance))
}

# sort columns
memoSort <- function(M, sortGenes=TRUE) {
  if(sortGenes){
    geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix;
  } else {
    geneOrder <- 1:nrow(M)
  }
  scoreCol <- function(x) {
    score <- 0;
    for(i in 1:length(x)) {
      if(x[i]) {
        score <- score + 2^(20+length(x)-i);
        break
      }
    }
    score <- score + sum(x*(length(x):1))
    return(score);
  }
  scores <- apply(M[geneOrder, ], 2, scoreCol);
  sampleOrder <- order(scores, decreasing = TRUE)# sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
  return(M[geneOrder, sampleOrder]);
}

M <- unmeltHallmarks(geneset_path)

oncomatrix_origin <- M[[1]]
oncomatrix_origin_numeric <- M[[2]]
patient_cancer_correspondance <- M[[3]]

# filter out unfrequently switched/mutated genes
if (filter_out_unfrequent){
  pct = apply(oncomatrix_origin, 1, function(x) sum(!grepl("^$", x))/length(x))*100
  oncomatrix_origin <- oncomatrix_origin[pct >=5,]
  oncomatrix_origin_numeric <- oncomatrix_origin_numeric[pct >=5,]
}

oncomatrix_origin_numeric <- memoSort(oncomatrix_origin_numeric)
oncomatrix_origin <- oncomatrix_origin[rownames(oncomatrix_origin_numeric),colnames(oncomatrix_origin_numeric)]
oncomatrix <- oncomatrix_origin[,as.logical(colSums(oncomatrix_origin_numeric))]

altered = ncol(oncomatrix)/ncol(oncomatrix_origin)

type_col = c("SWITCH" = "#2E80B2", "MUT" = "#B2741D")
type_name = c("SWITCH" = "Isoform switch", "MUT" = "Mutation")

add_oncoprint = function(type, x, y, width, height) {
  if(any(type %in% "")) {
    grid.rect(x, y, width - unit(0.5, "mm"), height - unit(1, "mm"), gp = gpar(col = NA, fill = "#CCCCCC"))
  }
  if(any(type %in% "SWITCH")) {
    grid.rect(x, y, width - unit(0.5, "mm"), height - unit(1, "mm"), gp = gpar(col = NA, fill = type_col["SWITCH"]))
  }
  if(any(type %in% "MUT")) {
    grid.rect(x, y, width - unit(0.5, "mm"), height*(1/3), gp = gpar(col = NA, fill = type_col["MUT"]))
  }
}

#####################################################################
# row annotation which shows percent of mutations in all samples
pct = apply(oncomatrix_origin, 1, function(x) sum(!grepl("^$", x))/length(x))*100
pct = paste0(round(pct),"%")
ha_pct = rowAnnotation(pct = anno_text(pct, which = "row"), width = grobWidth(textGrob("100%", gp = gpar(fontsize = 10))))

#####################################################################
# row annotation which is a barplot
anno_row_bar = function(index) {
  n = length(index)
  tb = apply(oncomatrix[index, ], 1, function(x) {
    x = unlist(strsplit(x, ";"))
    x = x[!grepl("^$", x)]
    x = sort(x)
    table(x)
  })
  max_count = max(sapply(tb, sum))
  pushViewport(viewport(xscale = c(0, max_count*1.1), yscale = c(0.5, n + 0.5)))
  for(i in seq_along(tb)) {
    if(length(tb[[i]])) {
      x = cumsum(tb[[i]])
        
      # row order is from top to end while coordinate of y is from bottom to top
      # so here we need to use n-i+1
      grid.rect(x, n-i+1, width = tb[[i]], height = 0.8, default.units = "native", just = "right", gp = gpar(col = NA, fill = type_col[names(tb[[i]])] ))
    }
  }
  breaks = grid.pretty(c(0, max_count))
  grid.xaxis(at = breaks, label = breaks, main = FALSE, gp = gpar(fontsize = 10))
  upViewport()
}

# row annotation which is a barplot
my_anno_row_bar = function(index) {
  n = length(index)
  tb = apply(oncomatrix[index, ], 1, function(x) {
    x = unlist(strsplit(x, ";"))
    x = x[!grepl("^$", x)]
    x = sort(x)
    table(x)
  })
  max_count = max(sapply(tb, sum))
  pushViewport(viewport(xscale = c(0, max_count*1.1), yscale = c(0.5, n + 0.5)))
  for(i in ncol(tb)) {
    if(length(tb[,i])) {
      x = sum(tb[,i])
      
      # row order is from top to end while coordinate of y is from bottom to top
      # so here we need to use n-i+1
      grid.rect(x, n-i+1, width = tb[,i], height = 0.8, default.units = "native", just = "right", gp = gpar(col = NA, fill = type_col[names(tb[,i])] ))
    }
  }
  breaks = grid.pretty(c(0, max_count))
  grid.xaxis(at = breaks, label = breaks, main = FALSE, gp = gpar(fontsize = 10))
  upViewport()
}

ha_row_bar = rowAnnotation(row_bar = my_anno_row_bar, width = unit(4, "cm"))
###################################################################
# column annotation which is also a barplot
anno_column_bar = function(index) {
  n = length(index)
  tb = apply(oncomatrix[, index], 2, function(x) {
    x = unlist(strsplit(x, ";"))
    x = x[!grepl("^$", x)]
    x = sort(x)
    table(x)
  })
  max_count = max(sapply(tb, sum))
  pushViewport(viewport(yscale = c(0, max_count*1.1), xscale = c(0.5, n + 0.5)))
  for(i in seq_along(tb)) {
    if(length(tb[[i]])) {
      y = cumsum(tb[[i]])
      grid.rect(i, y, height = tb[[i]], width = 0.8, default.units = "native", just = "top", gp = gpar(col = NA, fill = type_col[names(tb[[i]])]))
    }
  }
  breaks = grid.pretty(c(0, max_count))
  grid.yaxis(at = breaks, label = breaks, gp = gpar(fontsize = 10))
  upViewport()
}

# heatmap annotation with cancer type
colorPalette <- list( Cancer=c("brca"="#CC79A7","coad"="#636363","hnsc"="#F0E442","kich"="#006D2C","kirc"="#31A354","kirp"="#74C476","lihc"="#FC8D59","luad"="#08519C","lusc"="#3182BD","prad"="#D55E00","thca"="#5E3C99","coad-hyper"="#000000","coad-hypo"="#969696","brca-basal"="#EDF8FB","brca-her2"="#B3CDE3","brca-luminala"="#8C96C6","brca-luminalb"="#88419D","expression-up"="#D94701","expression-down"="#2171B5","psi-up"="#FD8D3C","psi-down"="#6BAED6","odds"="#FD8D3C","a3"="#33A02C","a5"="#E31A1C","mx"="#FF7F00","ri"="#6A3D9A","se"="#B15928","normal"="#377EB8","tumor"="#E41A1C"))

cancer_df = data.frame(Cancer = as.character(patient_cancer_correspondance[colnames(oncomatrix)]))
colorPalette$Cancer <- colorPalette$Cancer[unique(cancer_df$Cancer)]

ha_column_bar = HeatmapAnnotation(column_bar = anno_column_bar, which = "column",
                                  df = cancer_df, col=colorPalette)

#####################################################################
# the main matrix
ht = Heatmap(oncomatrix, rect_gp = gpar(type = "none"), cell_fun = function(j, i, x, y, width, height, fill) {
  type = unique(strsplit(oncomatrix[i, j], ";")[[1]])
  add_oncoprint(type, x, y, width, height)
}, row_names_gp = gpar(fontsize = 10), show_column_names = FALSE, show_heatmap_legend = FALSE,
top_annotation = ha_column_bar, top_annotation_height = unit(2, "cm"))
#ht_list = ha_pct + ht + ha_row_bar
ht_list = ha_pct + ht
#########################################################
# legend
legend = legendGrob(labels = type_name[names(type_col)], pch = 15, gp = gpar(col = type_col), nrow = 1)
pdf(paste0(hallmark_filename,".pdf"), width = 10, height = 10)
draw(ht_list, newpage = FALSE, annotation_legend_side = "bottom", annotation_legend_list = list(legend), column_title = qq(hallmark_name), column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()
