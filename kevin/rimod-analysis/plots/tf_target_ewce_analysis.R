################
# EWCE analysis of TF targets
################

#####
# Perform EWCE enrichment analysis with candidate TF targets
#####
library(EWCE)
library(stringr)
library(biomaRt)
library(pheatmap)
library(viridis)
library(extrafont)
loadfonts()

setwd("/Users/kevin/dzne/rimod_package/integrative_analysis/tf_cage_rnaseq/")

# mart
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")


# Generate Celltype Data from Darmanis dataset
mat <- read.table("/Users/kevin/dzne/rimod_package/analysis/deconvolution/ewce_analysis/GSE67835_norm_counts_all.txt",
                  sep="\t", header=T, row.names = 1)
ct <- read.table("/Users/kevin/dzne/rimod_package/analysis/deconvolution/ewce_analysis/GSE67835_celltypes.txt",
                 sep="\t", header=T, row.names = 1)
mat <- t(mat)

annotLevel <- list(l1=ct$Celltype)
ct_data <- generate.celltype.data(exp=mat, annotLevels = annotLevel, groupName = "Darmanis", no_cores=1)
load(ct_data)

# Define Background set
gene.mat <- read.table("/Users/kevin/dzne/rimod_package/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_vst_values_2019-10-23_13.33.11.txt", sep="\t", header=T, row.names = 1)
genes <- row.names(gene.mat)
genes <- str_split(genes, pattern="[.]", simplify = T)[,1]
rownames(gene.mat) <- genes

# get HGNC
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=genes, mart=ensembl)
gene.mat <- merge(gene.mat, bm, by.x="row.names", by.y="ensembl_gene_id")
gene.mat <- gene.mat[!duplicated(gene.mat$hgnc_symbol),]
gene.mat <- gene.mat[!gene.mat$hgnc_symbol == "",]
gene.mat <- na.omit(gene.mat)
genes <- gene.mat$hgnc_symbol
rownames(gene.mat) <- gene.mat$hgnc_symbol
gene.mat <- gene.mat[, c(-1, -ncol(gene.mat))]




# Get enrichment results for all TFs
perform_tf_ewce <- function(tfs, tf.df){
  enrichment.list <- list()
  for (tf in tfs){
    tmp <- tf.df[tf.df$TF == tf,]
    targets <- str_split(tmp$Overlapping_Genes, pattern=',')[[1]]
    targets <- targets[!duplicated(targets)]
    targets <- targets[targets %in% genes]
    res <- bootstrap.enrichment.test(ctd, hits=targets, bg = genes, genelistSpecies = "human", sctSpecies = "human", reps=1000)$result
    
    enrichment.list$tmp <- res
    names(enrichment.list)[length(enrichment.list)] <- tf
  }
  return(enrichment.list)
}



###
# MAPT
###
# UP TFs
tf.df <- read.table("MAPT_common_TFs_up.txt", sep="\t", header=T, stringsAsFactors = F)
tfs <- tf.df$TF
mapt.up.enrichment <- perform_tf_ewce(tfs = tfs, tf.df = tf.df)

# Down TFs
tf.df <- read.table("MAPT_common_TFs_down.txt", sep="\t", header=T, stringsAsFactors = F)
tfs <- tf.df$TF
mapt.down.enrichment <- perform_tf_ewce(tfs = tfs, tf.df = tf.df)



###
# GRN
###
# UP TFs
tf.df <- read.table("GRN_common_TFs_up.txt", sep="\t", header=T, stringsAsFactors = F)
tfs <- tf.df$TF
grn.up.enrichment <- perform_tf_ewce(tfs = tfs, tf.df = tf.df)

# Down TFs
tf.df <- read.table("GRN_common_TFs_down.txt", sep="\t", header=T, stringsAsFactors = F)
tfs <- tf.df$TF
grn.down.enrichment <- perform_tf_ewce(tfs = tfs, tf.df = tf.df)


###
# C9orf72
###
# UP TFs
tf.df <- read.table("C9orf72_common_TFs_up.txt", sep="\t", header=T, stringsAsFactors = F)
tfs <- tf.df$TF
c9.up.enrichment <- perform_tf_ewce(tfs = tfs, tf.df = tf.df)

# Down TFs
tf.df <- read.table("C9orf72_common_TFs_down.txt", sep="\t", header=T, stringsAsFactors = F)
tfs <- tf.df$TF
c9.down.enrichment <- perform_tf_ewce(tfs = tfs, tf.df = tf.df)


#====================#



####
# make heatmaps
####
library(pheatmap)
library(viridis)

makeDataFrame <- function(enr){
  df <- data.frame(enr[[1]]$p)
  for (i in 2:length(enr)) {
    tmp <- data.frame(enr[[i]]$p)
    df <- cbind(df, tmp)
  }
  colnames(df) <- names(enr)
  rownames(df) <- rownames(enr[[1]])
  return(df)
}

# GRN
# up
grn.df <- makeDataFrame(grn.up.enrichment)
grn.df <- grn.df[-nrow(grn.df),]
grn.df <- t(grn.df)
pheatmap(grn.df, color = viridis(200, option="D"), cluster_rows = F, cluster_cols = F, angle_col = "90",
         height = 5, width = 4, filename = "ewce_heatmap_grn_up.png")

# down
grn.df <- makeDataFrame(grn.down.enrichment)
grn.df <- grn.df[-nrow(grn.df),]
grn.df <- t(grn.df)
pheatmap(grn.df, color = viridis(200, option="D"), cluster_rows = F, cluster_cols = F, angle_col = "90",
         height = 5, width = 4, filename = "ewce_heatmap_grn_down.png")

# MAPT
# up
mapt.df <- makeDataFrame(mapt.up.enrichment)
mapt.df <- mapt.df[-nrow(mapt.df),]
mapt.df <- t(mapt.df)
pheatmap(mapt.df, color = viridis(200, option="D"), cluster_rows = F, cluster_cols = F, angle_col = "90",
         height = 5, width = 4, filename = "ewce_heatmap_mapt_up.png")

# down
mapt.df <- makeDataFrame(mapt.down.enrichment)
mapt.df <- mapt.df[-nrow(mapt.df),]
mapt.df <- t(mapt.df)
pheatmap(mapt.df, color = viridis(200, option="D"), cluster_rows = F, cluster_cols = F, angle_col = "90",
         height = 5, width = 4, filename = "ewce_heatmap_mapt_down.png")