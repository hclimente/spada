library(ComplexHeatmap)
library(GetoptLong)

# get matrix
unmeltHallmarks <- function(path){
  patient_info <- read.delim(path)
  patients <- unique(patient_info$patient)
  genes <- unique(patient_info$gene)
  unmelt_table <- as.data.frame(matrix(data="",nrow=length(patients),ncol=length(genes)),stringsAsFactors=FALSE)
  unmelt_table_numeric <- as.data.frame(matrix(data=0,nrow=length(patients),ncol=length(genes)),stringsAsFactors=FALSE)
  
  colnames(unmelt_table) <- genes
  rownames(unmelt_table) <- patients
  
  colnames(unmelt_table_numeric) <- genes
  rownames(unmelt_table_numeric) <- patients
  
  
  for (i in 1:nrow(patient_info)){
    thisPatient <- patient_info$patient[i]
    thisGene <- patient_info$gene[i]
  
    selectCol <- colnames(unmelt_table) == thisGene
    selectRow <- rownames(unmelt_table) == thisPatient
    
    thisAlteration <- paste0(patient_info$alteration[i],";")
    currentValue <- unmelt_table[selectRow,selectCol]
    thisAlteration <- paste0(currentValue,thisAlteration)
    
    unmelt_table[selectRow,selectCol] <- as.character(thisAlteration)
    unmelt_table_numeric[selectRow,selectCol] <- 1
      
  }
  unmelt_table <- t(unmelt_table)
  unmelt_table_numeric <- t(unmelt_table_numeric)
  return(list(unmelt_table,unmelt_table_numeric))
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

M <- unmeltHallmarks("HALLMARK_TEST.txt")

oncomatrix_origin <- M[[1]]
oncomatrix_origin_numeric <- M[[2]]

oncomatrix_origin_numeric <- memoSort(oncomatrix_origin_numeric)
oncomatrix_origin <- oncomatrix_origin[rownames(oncomatrix_origin_numeric),colnames(oncomatrix_origin_numeric)]

oncomatrix <- oncomatrix_origin[,as.logical(colSums(oncomatrix_origin_numeric))]

altered = ncol(oncomatrix)/ncol(oncomatrix_origin)

type_col = c("SWITCH" = "red", "MUT" = "blue")
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
      grid.rect(x, n-i+1, width = tb[[i]], height = 0.8, default.units = "native", just = "right", gp = gpar(col = NA, fill = type_col[names(tb[[i]])]))
    }
  }
  breaks = grid.pretty(c(0, max_count))
  grid.xaxis(at = breaks, label = breaks, main = FALSE, gp = gpar(fontsize = 10))
  upViewport()
}
ha_row_bar = rowAnnotation(row_bar = anno_row_bar, width = unit(4, "cm"))
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
ha_column_bar = HeatmapAnnotation(column_bar = anno_column_bar, which = "column")
#####################################################################
# the main matrix
ht = Heatmap(oncomatrix, rect_gp = gpar(type = "none"), cell_fun = function(j, i, x, y, width, height, fill) {
  type = unique(strsplit(oncomatrix[i, j], ";")[[1]])
  add_oncoprint(type, x, y, width, height)
}, row_names_gp = gpar(fontsize = 10), show_column_names = FALSE, show_heatmap_legend = FALSE,
top_annotation = ha_column_bar, top_annotation_height = unit(2, "cm"))
ht_list = ha_pct + ht + ha_row_bar
#########################################################
# legend
legend = legendGrob(labels = type_name[names(type_col)], pch = 15, gp = gpar(col = type_col), nrow = 1)
#pdf("oncoprint.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, annotation_legend_side = "bottom", annotation_legend_list = list(legend), column_title = qq("OncoPrint"), column_title_gp = gpar(fontsize = 12, fontface = "bold"))
#dev.off()