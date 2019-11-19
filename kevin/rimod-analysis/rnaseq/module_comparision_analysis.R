##############################
# Module Comparison
##############################
library(stringr)
library(biomaRt)
library(igraph)
setwd("~/rimod/RNAseq/analysis/module_integration/")

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


getModule <- function(modules, mod){
  genes <- str_split(modules[modules$CLUSTER_NAME == mod,]$CLUSTER_GENES, pattern=",")[[1]]
  return(genes)
}

getSymbols <- function(module, mart){
  bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=module$V1, mart=mart)
  genes <- bm$hgnc_symbol
  genes <- genes[!genes == ""]
  return(genes)
}

# Extrac gene accessions
getGenes <- function(x){
  genes <- as.character(x$GencodeCompV12_NAME)
  genes <- genes[!genes == ""]
  genes <- as.character(sapply(genes, function(y){strsplit(y, split="[;]")[[1]][[1]]}))
  genes <- genes[!duplicated(genes)]
  return(genes)
}


###
# GRN modules
###
grn.up.modules <- read.table("../human_base/rnaseq_grn_filtered_up_modules.txt", header=T, stringsAsFactors = F)
grn.down.modules <- read.table("../human_base/rnaseq_grn_filtered_down_modules.txt", header=T, stringsAsFactors = F)
grn.up.cnc <- read.table("../human_base/rnaseq_grn_filtered_up_connections.txt", header=T, stringsAsFactors = F)
grn.down.cnc <- read.table("../human_base/rnaseq_grn_filtered_down_connections.txt", header=T, stringsAsFactors = F)

##
# Down-regulated modules, miRNAs
grn.mir <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/GRN_miRNA_target_edge_table.txt", sep="\t", header=T)

for (m in grn.down.modules$CLUSTER_NAME) {
  print(m)
  genes <- getModule(grn.down.modules, m)
  print(length(genes))
  tmp <- grn.mir[grn.mir$targets %in% genes,]
  no_targets <- nrow(tmp[!duplicated(tmp$targets),])
  print(no_targets)
  print(table(tmp$mirna))
  print(no_targets/ length(genes))
}

# make graph for module M4
grn.m4 <- getModule(grn.down.modules, "M4")
grn.m4.cnc <- grn.down.cnc[grn.down.cnc$SOURCE %in% grn.m4,]
grn.m4.cnc <- grn.m4.cnc[grn.m4.cnc$TARGET %in% grn.m4,]
# create edges
edges <- c()
for (i in 1:nrow(grn.m4.cnc)) {
  e <- c(as.character(grn.m4.cnc[i,1]), as.character(grn.m4.cnc[i,2]))
  edges <- c(edges, e)
}
# add miRNA-target edges
grn.mir.m4 <- grn.mir[grn.mir$targets %in% grn.m4,]
for (i in 1:nrow(grn.mir.m4)) {
  e <- c(as.character(grn.mir.m4[i,1]), as.character(grn.mir.m4[i,2]))
  edges <- c(edges, e)
}

g <- graph(edges = edges)
# add weights to the graph
weights <- c(grn.m4.cnc$WEIGHT, rep(0.5, nrow(grn.mir.m4)))
E(g)$weight = weights

# Add type to the graph
node_type = rep("Gene", length(V(g)$name))
node_type[grepl("miR", V(g)$name)] <- "miRNA"
V(g)$type = node_type

write_graph(g, file = "GRN_M4_network.gml", format = "gml")

#===========================================#

###
# MAPT modules
###
mapt.up.modules <- read.table("../human_base/rnaseq_mapt_filtered_up_modules.txt", header=T, stringsAsFactors = F)
mapt.down.modules <- read.table("../human_base/rnaseq_mapt_filtered_down_modules.txt", header=T, stringsAsFactors = F)
mapt.up.cnc <- read.table("../human_base/rnaseq_mapt_filtered_up_connections.txt", header=T, stringsAsFactors = F)
mapt.down.cnc <- read.table("../human_base/rnaseq_mapt_filtered_down_connections.txt", header=T, stringsAsFactors = F)

##
# Down-regulated modules, miRNAs
mapt.mir <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/MAPT_miRNA_target_edge_table.txt", sep="\t", header=T)

for (m in mapt.down.modules$CLUSTER_NAME) {
  print(m)
  genes <- getModule(mapt.down.modules, m)
  print(length(genes))
  tmp <- mapt.mir[mapt.mir$targets %in% genes,]
  no_targets <- nrow(tmp[!duplicated(tmp$targets),])
  print(no_targets)
  print(table(tmp$mirna))
  print(no_targets/ length(genes))
}

# make graph for module M4
mapt.m4 <- getModule(mapt.down.modules, "M4")
mapt.m4.cnc <- mapt.down.cnc[mapt.down.cnc$SOURCE %in% mapt.m4,]
mapt.m4.cnc <- mapt.m4.cnc[mapt.m4.cnc$TARGET %in% mapt.m4,]
# create edges
edges <- c()
for (i in 1:nrow(mapt.m4.cnc)) {
  e <- c(as.character(mapt.m4.cnc[i,1]), as.character(mapt.m4.cnc[i,2]))
  edges <- c(edges, e)
}
# add miRNA-target edges
mapt.mir.m4 <- mapt.mir[mapt.mir$targets %in% mapt.m4,]
for (i in 1:nrow(mapt.mir.m4)) {
  e <- c(as.character(mapt.mir.m4[i,1]), as.character(mapt.mir.m4[i,2]))
  edges <- c(edges, e)
}

g <- graph(edges = edges)
# add weights to the graph
weights <- c(mapt.m4.cnc$WEIGHT, rep(0.5, nrow(mapt.mir.m4)))
E(g)$weight = weights

# Add type to the graph
node_type = rep("Gene", length(V(g)$name))
node_type[grepl("miR", V(g)$name)] <- "miRNA"
V(g)$type = node_type

write_graph(g, file = "MAPT_M4_network.gml", format = "gml")

###
# Check for epigenetic enrichment in modules
###
# MAPT
mapt.met <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_mapt.ndc_quant.txt", sep="\t")
mapt.met <- mapt.met[mapt.met$adj.P.Val <= 0.01,]
mapt.met.genes <- getGenes(mapt.met)
# iterate over modules
for (m in mapt.down.modules$CLUSTER_NAME) {
  print(m)
  genes <- getModule(mapt.down.modules, m)
  print(length(genes))
  tmp <- mapt.met.genes[mapt.met.genes %in% genes]
  print(length(tmp))
  print(length(tmp) / length(genes))
}


####
# Check for TFs
####
# mapt up
tf.mapt <- read.table("~/rimod/CAGE/cage_analysis/chea_tf_regulators/frontal/Integrated_meanRank_MAPT_up.tsv", sep="\t", header=T, stringsAsFactors = F)
tf.mapt <- tf.mapt[1:50,]

for (m in mapt.up.modules$CLUSTER_NAME) {
  print(m)
  genes <- getModule(mapt.up.modules, m)
  
  gene_overlaps <- c()
  for (i in 1:nrow(tf.mapt)) {
    ovl <- as.character(tf.mapt$Overlapping_Genes[i])
    ovl <- str_split(ovl, pattern=",")[[1]]
    ovl <- ovl[ovl %in% genes]
    gene_overlaps <- c(gene_overlaps, length(ovl))
  }
  tmp <- tf.mapt
  tmp$overlap <- gene_overlaps
  tmp <- tmp[order(tmp$overlap, decreasing = T),]
  tmp <- tmp[1:10,]
  print(tmp$TF)
  print(tmp$overlap)
}

# mapt down
tf.mapt <- read.table("~/rimod/CAGE/cage_analysis/chea_tf_regulators/frontal/Integrated_meanRank_MAPT_down.tsv", sep="\t", header=T, stringsAsFactors = F)
tf.mapt <- tf.mapt[1:50,]

for (m in mapt.down.modules$CLUSTER_NAME) {
  print(m)
  genes <- getModule(mapt.down.modules, m)
  print(length(genes))
  gene_overlaps <- c()
  for (i in 1:nrow(tf.mapt)) {
    ovl <- as.character(tf.mapt$Overlapping_Genes[i])
    ovl <- str_split(ovl, pattern=",")[[1]]
    ovl <- ovl[ovl %in% genes]
    gene_overlaps <- c(gene_overlaps, length(ovl))
  }
  tmp <- tf.mapt
  tmp$overlap <- gene_overlaps
  tmp <- tmp[order(tmp$overlap, decreasing = T),]
  tmp <- tmp[1:10,]
  print(tmp$TF)
  print(tmp$overlap)
}

# Chechk for overlap with promotor shifting
ps <- read.table("~/rimod/CAGE/cage_analysis/promotor_shifting/frontal/mapt_promotor_shifting_genes_fro.txt", sep="\t", header=T)
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=ps$x, mart=ensembl)
ps <- bm$hgnc_symbol

for (m in mapt.down.modules$CLUSTER_NAME) {
  print(m)
  genes <- getModule(mapt.down.modules, m)
  print(length(genes))
  ovl = intersect(genes, ps)
  print(length(ovl))
  print(length(ovl) / length(genes))
}

# Check for alternative splicing
as <- read.table("~/rimod/RNAseq/as_analysis/majiq/mapt_AS_genes_dPSI_0.2.txt", sep="\t", header=T)
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=as$x, mart=ensembl)
as <- bm$hgnc_symbol
for (m in mapt.down.modules$CLUSTER_NAME) {
  print(m)
  genes <- getModule(mapt.down.modules, m)
  print(length(genes))
  ovl = intersect(genes, as)
  print(length(ovl))
  print(length(ovl) / length(genes))
}





