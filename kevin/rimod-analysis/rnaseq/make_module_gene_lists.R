#####################
# Make gene lists for every module for further pathway enrichment analysis
#####################
library(stringr)
setwd("~/rimod/RNAseq/analysis/human_base/")

# grn up
grn_up <- read.table("rnaseq_grn_filtered_up_modules.txt", sep="\t", header=T)
for (m in grn_up$CLUSTER_NAME) {
  print(m)
  mod <- as.character(grn_up[grn_up$CLUSTER_NAME == m,]$CLUSTER_GENES)
  mod <- str_split(mod, pattern=",")[[1]]
  write.table(mod, paste("module_genes/grn_up_", m, ".txt", sep=""), quote=F, col.names = F, row.names = F)
}


# grn down
grn_up <- read.table("rnaseq_grn_filtered_down_modules.txt", sep="\t", header=T)
for (m in grn_up$CLUSTER_NAME) {
  print(m)
  mod <- as.character(grn_up[grn_up$CLUSTER_NAME == m,]$CLUSTER_GENES)
  mod <- str_split(mod, pattern=",")[[1]]
  write.table(mod, paste("module_genes/grn_down_", m, ".txt", sep=""), quote=F, col.names = F, row.names = F)
}


# mapt up
mapt <- read.table("rnaseq_mapt_filtered_up_modules.txt", sep="\t", header=T)
for (m in mapt$CLUSTER_NAME) {
  print(m)
  mod <- as.character(mapt[mapt$CLUSTER_NAME == m,]$CLUSTER_GENES)
  mod <- str_split(mod, pattern=",")[[1]]
  write.table(mod, paste("module_genes/mapt_up_", m, ".txt", sep=""), quote=F, col.names = F, row.names = F)
}


# mapt down
mapt <- read.table("rnaseq_mapt_filtered_down_modules.txt", sep="\t", header=T)
for (m in mapt$CLUSTER_NAME) {
  print(m)
  mod <- as.character(mapt[mapt$CLUSTER_NAME == m,]$CLUSTER_GENES)
  mod <- str_split(mod, pattern=",")[[1]]
  write.table(mod, paste("module_genes/mapt_down_", m, ".txt", sep=""), quote=F, col.names = F, row.names = F)
}


