################
# EWCE analysis of miRNA targets
################

#####
# Perform EWCE enrichment analysis with human snRNA-seq data
#####
library(EWCE)
library(stringr)
library(biomaRt)
library(pheatmap)
library(viridis)
library(extrafont)
loadfonts()

setwd("/Users/kevin/dzne/rimod_package/analysis/deconvolution/ewce_analysis/")

# Get genes from a module
getModule <- function(modules, mod){
  genes <- str_split(modules[modules$CLUSTER_NAME == mod,]$CLUSTER_GENES, pattern=",")[[1]]
  return(genes)
}
# mart
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")



# Generate Celltype Data from Darmanis dataset
mat <- read.table("GSE67835_norm_counts_all.txt",
                  sep="\t", header=T, row.names = 1)
ct <- read.table("GSE67835_celltypes.txt",
                 sep="\t", header=T, row.names = 1)
mat <- t(mat)
# keep only genes that are expressed commonly
#rs <- apply(mat, 1, sum)
#keep <- rs > 100
#mat <- mat[keep,]

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

# Get enrichment results for all miRNAs
perform_mirna_ewce <- function(mir.genes, mirs){
  enrichment.list <- list()
  for (mir in mir.genes){
    tmp <- mirs[mirs$V1 == mir,]
    targets <- tmp$V2
    bm <- getBM(attributes = c("hgnc_symbol", "refseq_mrna"), filters="refseq_mrna", values=targets, mart=ensembl)
    targets <- bm$hgnc_symbol
    targets <- targets[!duplicated(targets)]
    targets <- targets[targets %in% genes]
    res <- bootstrap.enrichment.test(ctd, hits=targets, bg = genes, genelistSpecies = "human", sctSpecies = "human", reps=1000)$result
    
    enrichment.list$tmp <- res
    names(enrichment.list)[length(enrichment.list)] <- mir
  }
  return(enrichment.list)
}



###
# MAPT
###
up.targets <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/target_mrna_correlation_analysis_0819/MAPT_downMir_correlated_targets_Refseq.txt", sep="\t", header=T)
bm <- getBM(attributes = c("hgnc_symbol", "refseq_mrna"), filters="refseq_mrna", values=up.targets$x, mart=ensembl)
up.targets <- bm$hgnc_symbol
up.targets <- up.targets[!duplicated(up.targets)]

down.targets <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/target_mrna_correlation_analysis_0819/MAPT_upMir_correlated_targets_Refseq.txt", sep="\t", header=T)
bm <- getBM(attributes = c("hgnc_symbol", "refseq_mrna"), filters="refseq_mrna", values=down.targets$x, mart=ensembl)
down.targets <- bm$hgnc_symbol
down.targets <- down.targets[!duplicated(down.targets)]
down.targets <- down.targets[down.targets %in% genes]


res.mapt.up <- bootstrap.enrichment.test(ctd, hits=up.targets, bg = genes, genelistSpecies = "human", sctSpecies = "human", reps=1000)
res.mapt.down <- bootstrap.enrichment.test(ctd, hits=down.targets, bg = genes, genelistSpecies = "human", sctSpecies = "human", reps=1000)


# Perform enrichment for sinlge miRNAs
mirs <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/mirna_target_analysis_0719/MAPT_DEG_targets.txt", sep="\t", header=T, stringsAsFactors = F)
mir.genes <- as.character(unique(mirs$V1))

mapt.mir.enrichment <- perform_mirna_ewce(mir.genes=mir.genes, mirs=mirs)



#====================#


###
# GRN
###
up.targets <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/target_mrna_correlation_analysis_0819/GRN_downMir_correlated_targets_Refseq.txt", sep="\t", header=T)
bm <- getBM(attributes = c("hgnc_symbol", "refseq_mrna"), filters="refseq_mrna", values=up.targets$x, mart=ensembl)
up.targets <- bm$hgnc_symbol
up.targets <- up.targets[!duplicated(up.targets)]

down.targets <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/target_mrna_correlation_analysis_0819/GRN_upMir_correlated_targets_Refseq.txt", sep="\t", header=T)
bm <- getBM(attributes = c("hgnc_symbol", "refseq_mrna"), filters="refseq_mrna", values=down.targets$x, mart=ensembl)
down.targets <- bm$hgnc_symbol
down.targets <- down.targets[!duplicated(down.targets)]


res.grn.up <- bootstrap.enrichment.test(ctd, hits=up.targets, bg = genes, genelistSpecies = "human", sctSpecies = "human", reps=1000)
res.grn.down <- bootstrap.enrichment.test(ctd, hits=down.targets, bg = genes, genelistSpecies = "human", sctSpecies = "human", reps=1000)


# Perform enrichment for sinlge miRNAs
mirs <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/mirna_target_analysis_0719/GRN_DEG_targets.txt", sep="\t", header=T, stringsAsFactors = F)
mir.genes <- as.character(unique(mirs$V1))

grn.mir.enrichment <- perform_mirna_ewce(mir.genes=mir.genes, mirs=mirs)

#====================#


###
# C9orf72
###
up.targets <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/target_mrna_correlation_analysis_0819/C9_downMir_correlated_targets_Refseq.txt", sep="\t", header=T)
bm <- getBM(attributes = c("hgnc_symbol", "refseq_mrna"), filters="refseq_mrna", values=up.targets$x, mart=ensembl)
up.targets <- bm$hgnc_symbol
up.targets <- up.targets[!duplicated(up.targets)]

down.targets <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/target_mrna_correlation_analysis_0819/C9_upMir_correlated_targets_Refseq.txt", sep="\t", header=T)
bm <- getBM(attributes = c("hgnc_symbol", "refseq_mrna"), filters="refseq_mrna", values=down.targets$x, mart=ensembl)
down.targets <- bm$hgnc_symbol
down.targets <- down.targets[!duplicated(down.targets)]


res.c9.up <- bootstrap.enrichment.test(ctd, hits=up.targets, bg = genes, genelistSpecies = "human", sctSpecies = "human", reps=1000)
res.c9.down <- bootstrap.enrichment.test(ctd, hits=down.targets, bg = genes, genelistSpecies = "human", sctSpecies = "human", reps=1000)

# Perform enrichment for sinlge miRNAs
mirs <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/mirna_target_analysis_0719/C9_DEG_targets.txt", sep="\t", header=T, stringsAsFactors = F)
mir.genes <- as.character(unique(mirs$V1))

C9.mir.enrichment <- perform_mirna_ewce(mir.genes=mir.genes, mirs=mirs)
#====================#



####
# make heatmaps
####
library(pheatmap)
library(viridis)
setwd("/Users/kevin/dzne/rimod_analysis/figure4/")

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
grn.df <- makeDataFrame(grn.mir.enrichment)
grn.df <- grn.df[-nrow(grn.df),]
grn.df <- t(grn.df)

pheatmap(grn.df, color = viridis(200, option="D"), cluster_rows = F, cluster_cols = F, angle_col = "90",
         height = 5, width = 4, filename = "ewce_heatmap_grn.png")

# MAPT
mapt.df <- makeDataFrame(mapt.mir.enrichment)
mapt.df <- mapt.df[-nrow(mapt.df),]
mapt.df <- t(mapt.df)

pheatmap(mapt.df, color = viridis(200, option="D"), cluster_rows = F, cluster_cols = F, angle_col = "90",
         height = 5, width = 4, filename = "ewce_heatmap_mapt.png")
# C9orf72
c9.df <- makeDataFrame(C9.mir.enrichment)
c9.df <- c9.df[-nrow(c9.df),]
c9.df <- t(c9.df)

pheatmap(c9.df, color = viridis(200, option="D"), cluster_rows = F, cluster_cols = F, angle_col = "90",
         height = 5, width = 4, filename = "ewce_heatmap_c9orf72.png")


# Only certain number of genes
n_mirs = 10
# mapt
deg <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/analysis_0719/deseq_result_mapt.ndc_frontal_smRNAseq.txt", sep="\t", header=T, row.names=1)
deg <- deg[deg$padj <= 0.05,]
deg <- deg[deg$log2FoldChange > 0.6,]
deg <- deg[order(deg$padj),]
mapt.df <- mapt.df[rownames(mapt.df) %in% rownames(deg),]
deg <- deg[rownames(deg) %in% rownames(mapt.df),]
mapt.df <- mapt.df[match(rownames(deg), rownames(mapt.df)),]
mapt.df <- mapt.df[1:n_mirs, ]
pheatmap(mapt.df, color = viridis(200, option="D"), cluster_rows = F, cluster_cols = F, angle_col = "90",
         height = 5, width = 4, filename = "upDegs20_ewce_heatmap_mapt.png")

# grn
deg <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/analysis_0719/deseq_result_grn.ndc_frontal_smRNAseq.txt", sep="\t", header=T, row.names=1)
deg <- deg[deg$padj <= 0.05,]
deg <- deg[deg$log2FoldChange > 0.6,]
deg <- deg[order(deg$padj),]
grn.df <- grn.df[rownames(grn.df) %in% rownames(deg),]
deg <- deg[rownames(deg) %in% rownames(grn.df),]
grn.df <- grn.df[match(rownames(deg), rownames(grn.df)),]
grn.df <- grn.df[1:n_mirs,]


pheatmap(grn.df, color = viridis(200, option="D"), cluster_rows = F, cluster_cols = F, angle_col = "90",
         height = 5, width = 4, filename = "upDegs20_ewce_heatmap_grn.png")

# c9orf72
deg <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/analysis_0719/deseq_result_c9.ndc_frontal_smRNAseq.txt", sep="\t", header=T, row.names=1)
deg <- deg[deg$padj <= 0.05,]
deg <- deg[deg$log2FoldChange > 0.6,]
deg <- deg[order(deg$padj),]
c9.df <- c9.df[rownames(c9.df) %in% rownames(deg),]
deg <- deg[rownames(deg) %in% rownames(c9.df),]
c9.df <- c9.df[match(rownames(deg), rownames(c9.df)),]
c9.df <- c9.df[1:n_mirs,]

pheatmap(c9.df, color = viridis(200, option="D"), cluster_rows = F, cluster_cols = F, angle_col = "90",
         height = 5, width = 4, filename = "upDegs20_ewce_heatmap_c9orf72.png")

