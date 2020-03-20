###############################
# Figure 3: only TF-module networks
################################
library(stringr)
library(igraph)
library(ggplot2)

setwd("~/rimod/paper/figures/figure3/")

####
# MAPT
####

# TF
tf.up <- read.table("~/rimod/integrative_analysis/tf_cage_rnaseq/MAPT_common_TFs_up.txt", sep="\t", header=T, stringsAsFactors = F)
tf.down <- read.table("~/rimod/integrative_analysis/tf_cage_rnaseq/MAPT_common_TFs_down.txt", sep="\t", header=T, stringsAsFactors = F)
deg <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/MAPT_cell_composition_filterd_DEGs.txt", sep="\t", header=T, stringsAsFactors = F)

# Modules
mod.up <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_mapt_filtered_up_modules.txt", sep="\t", header=T)
mod.down <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_mapt_filtered_down_modules.txt", sep="\t", header=T)




####
# Down-modules
down.mods <- list()
for (m in as.character(mod.down$CLUSTER_NAME)) {
  tmp <- mod.down[mod.down$CLUSTER_NAME == m,]
  genes <- as.character(tmp$CLUSTER_GENES)
  genes <- str_split(genes, pattern=",")[[1]]
  down.mods$tmp <- genes
  names(down.mods)[length(down.mods)] <- m
}


# Format Down-TF
tf.down <- tf.down[, c(3, 6)]
tf.down <- tf.down[tf.down$TF %in% deg$hgnc_symbol,]

####
# Up-modules
## Up-modules
up.mods <- list()
for (m in as.character(mod.up$CLUSTER_NAME)) {
  tmp <- mod.up[mod.up$CLUSTER_NAME == m,]
  genes <- as.character(tmp$CLUSTER_GENES)
  genes <- str_split(genes, pattern=",")[[1]]
  up.mods$tmp <- genes
  names(up.mods)[length(up.mods)] <- m
}

tf.up <- tf.up[, c(3, 6)]
tf.up <- tf.up[tf.up$TF %in% deg$hgnc_symbol,]

# 
tf.up <- deg[deg$hgnc_symbol %in% tf.up$TF,]
tf.down <- deg[deg$hgnc_symbol %in% tf.down$TF,]
tfs <- rbind(tf.up, tf.down)

df <- data.frame(TF = factor(tfs$hgnc_symbol, levels=as.character(tfs$hgnc_symbol)), LFC = tfs$log2FoldChange, Pvalue = tfs$padj)

p <- ggplot(df, aes(x=TF, y=LFC, fill = Pvalue)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_minimal()
p
ggsave("mapt_TFs_barplot.png", width=3, height=3)

####
# GRN
####

#####
# FTD-GRN
#####

# TF
tf.up <- read.table("~/rimod/integrative_analysis/tf_cage_rnaseq/GRN_common_TFs_up.txt", sep="\t", header=T, stringsAsFactors = F)
tf.down <- read.table("~/rimod/integrative_analysis/tf_cage_rnaseq/GRN_common_TFs_down.txt", sep="\t", header=T, stringsAsFactors = F)
deg <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/GRN_cell_composition_filterd_DEGs.txt", sep="\t", header=T, stringsAsFactors = F)

# Modules
mod.up <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_grn_filtered_up_modules.txt", sep="\t", header=T)
mod.down <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_grn_filtered_down_modules.txt", sep="\t", header=T)




####
# Down-modules
down.mods <- list()
for (m in as.character(mod.down$CLUSTER_NAME)) {
  tmp <- mod.down[mod.down$CLUSTER_NAME == m,]
  genes <- as.character(tmp$CLUSTER_GENES)
  genes <- str_split(genes, pattern=",")[[1]]
  down.mods$tmp <- genes
  names(down.mods)[length(down.mods)] <- m
}

tf.down <- tf.down[, c(3, 6)]
tf.down <- tf.down[tf.down$TF %in% deg$hgnc_symbol,]
## Up-modules
up.mods <- list()
for (m in as.character(mod.up$CLUSTER_NAME)) {
  tmp <- mod.up[mod.up$CLUSTER_NAME == m,]
  genes <- as.character(tmp$CLUSTER_GENES)
  genes <- str_split(genes, pattern=",")[[1]]
  up.mods$tmp <- genes
  names(up.mods)[length(up.mods)] <- m
}

tf.up <- tf.up[, c(3, 6)]
tf.up <- tf.up[tf.up$TF %in% deg$hgnc_symbol,]

# PLOTTING
# 
tf.up <- deg[deg$hgnc_symbol %in% tf.up$TF,]
tf.down <- deg[deg$hgnc_symbol %in% tf.down$TF,]
tfs <- rbind(tf.up, tf.down)

df <- data.frame(TF = factor(tfs$hgnc_symbol, levels=as.character(tfs$hgnc_symbol)), LFC = tfs$log2FoldChange, Pvalue = tfs$padj)

p <- ggplot(df, aes(x=TF, y=LFC, fill = Pvalue)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_minimal()
p
ggsave("grn_TFs_barplot.png", width=3, height=3)


