###############
# RNA-seq module genetics analysis
###############
library(stringr)

setwd("~/rimod/RNAseq/humanbase_genetics/")

# Load FTD-asscoiated mutations
ftd <- read.csv("frontotemporal_dementia_omim.csv", header=T, sep=";")
genes <- ftd$Approved.Symbol
genes <- as.character(genes[!duplicated(genes)])

# Load modules
mods <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_mapt_filtered_down_modules.txt", header=T)

test = mods$CLUSTER_GENES[1]
test <- str_split(test, pattern="[.]")[[1]]

for (i in 1:nrow(mods)){
  print(as.character(mods$CLUSTER_NAME[i]))
  mod.genes <- mods$CLUSTER_GENES[i]
  mod.genes <- str_split(mod.genes, pattern=",")[[1]]
  
  ovl <- intersect(mod.genes, genes)
  print(length(mod.genes))
  print(ovl)
}

# Perform enrichment analysis
ovl.len <- length(ovl)
n_genes = length(genes)
mod.genes.len <- length(mod.genes)
background = 20000


phyper(q=ovl.len -1, m=n_genes, n=)