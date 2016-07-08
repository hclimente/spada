source("~/smartas/pipeline/scripts/variablesAndFunctions.r")

pancancer.files <- "~/smartas/analyses/pancancer/"

all.proteome.muts <- list()
all.proteome.fts <- list()
all.switch.prosite <- list()
all.switch.pfam <- list()
for (cancer in cancerTypes){
  cancer.files <- paste0("~/smartas/analyses/",cancer,"/")
  proteome.muts.file <- paste0(cancer.files,"mutations/proteome_mutations.txt")
  proteome.fts.file <- paste0(cancer.files,"mutations/proteome_features.txt")
  switch.prosite.file <- paste0(cancer.files,"structural_analysis/prosite_analysis.tsv")
  switch.pfam.file <- paste0(cancer.files,"structural_analysis/interpro_analysis.tsv")
  
  # read mutations
  all.proteome.muts[[cancer]] <- read_tsv(proteome.muts.file)
  all.proteome.fts[[cancer]] <- read_tsv(proteome.fts.file)
  all.switch.prosite[[cancer]] <- read_tsv(switch.prosite.file)
  all.switch.pfam[[cancer]] <- read_tsv(switch.pfam.file)

}

all.proteome.muts <- do.call("rbind",all.proteome.muts)
proteome.muts.pan.file <- paste0(pancancer.files,"mutations/proteome_mutations.txt")
write_tsv(all.proteome.muts,proteome.muts.pan.file)

all.proteome.fts <- do.call("rbind",all.proteome.fts)
proteome.fts.pan.file <- paste0(pancancer.files,"mutations/proteome_features.txt")
write_tsv(all.proteome.fts,proteome.fts.pan.file)

all.switch.prosite <- do.call("rbind",all.switch.prosite)
switch.prosite.pan.file <- paste0(pancancer.files,"structural_analysis/prosite_analysis.tsv")
write_tsv(all.switch.prosite,switch.prosite.pan.file)

all.switch.pfam <- do.call("rbind",all.switch.pfam)
switch.pfam.pan.file <- paste0(pancancer.files,"structural_analysis/interpro_analysis.tsv")
write_tsv(all.switch.pfam,switch.pfam.pan.file)

system2('Rscript', c('~/smartas/pipeline/methods/mutated_features_analysis.R',pancancer.files))