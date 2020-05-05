setwd("~/rimod/RNAseq/analysis/human_base/module_genes/")

library(stringr)

# grn up
grn <- read.table("../rnaseq_grn_filtered_up_modules.txt", stringsAsFactors = F, header=T)
for (i in 1:nrow(grn)){
  m <- grn$CLUSTER_NAME[i]
  mod <- grn$CLUSTER_GENES[i]
  genes <- str_split(mod, pattern=",")[[1]]
  write.table(genes, paste0("GRN_",m,"_up_genes.txt"), quote=F, row.names=F, col.names=F)
}

# grn down
grn <- read.table("../rnaseq_grn_filtered_down_modules.txt", stringsAsFactors = F, header=T)
for (i in 1:nrow(grn)){
  m <- grn$CLUSTER_NAME[i]
  mod <- grn$CLUSTER_GENES[i]
  genes <- str_split(mod, pattern=",")[[1]]
  write.table(genes, paste0("GRN_",m,"_down_genes.txt"), quote=F, row.names=F, col.names=F)
}



# mapt up
mapt <- read.table("../rnaseq_mapt_filtered_up_modules.txt", stringsAsFactors = F, header=T)
for (i in 1:nrow(grn)){
  m <- mapt$CLUSTER_NAME[i]
  mod <- mapt$CLUSTER_GENES[i]
  genes <- str_split(mod, pattern=",")[[1]]
  write.table(genes, paste0("MAPT_",m,"_up_genes.txt"), quote=F, row.names=F, col.names=F)
}


# mapt up
mapt <- read.table("../rnaseq_mapt_filtered_down_modules.txt", stringsAsFactors = F, header=T)
for (i in 1:nrow(grn)){
  m <- mapt$CLUSTER_NAME[i]
  mod <- mapt$CLUSTER_GENES[i]
  genes <- str_split(mod, pattern=",")[[1]]
  write.table(genes, paste0("MAPT_",m,"_down_genes.txt"), quote=F, row.names=F, col.names=F)
}
