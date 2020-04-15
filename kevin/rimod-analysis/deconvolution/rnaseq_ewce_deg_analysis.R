#####
# Perform EWCE enrichment analysis with human snRNA-seq data
#####
library(EWCE)
library(stringr)
library(biomaRt)

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

###
# MAPT analysis
###
mapt.up <- read.table("/Users/kevin/dzne/rimod_package/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/DEGs_UP_mapt.ndc.txt", sep="\t", header=T)
mapt.down <- read.table("/Users/kevin/dzne/rimod_package/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/DEGs_DOWN_mapt.ndc.txt", sep="\t", header=T)
# transform to HGNC
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=mapt.up$x, mart=ensembl)
mapt.up <- bm$hgnc_symbol
mapt.up <- mapt.up[mapt.up %in% genes]
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=mapt.down$x, mart=ensembl)
mapt.down <- bm$hgnc_symbol
mapt.down <- mapt.down[mapt.down %in% genes]

res.mapt.up <- bootstrap.enrichment.test(ctd, hits=mapt.up, bg = genes, genelistSpecies = "human", sctSpecies = "human", reps=1000)
res.mapt.down <- bootstrap.enrichment.test(ctd, hits=mapt.down, bg = genes, genelistSpecies = "human", sctSpecies = "human", reps=1000)

#=======================#


###
# GRN analysis
###
mapt.up <- read.table("/Users/kevin/dzne/rimod_package/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/DEGs_UP_grn.ndc.txt", sep="\t", header=T)
mapt.down <- read.table("/Users/kevin/dzne/rimod_package/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/DEGs_DOWN_grn.ndc.txt", sep="\t", header=T)
# transform to HGNC
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=mapt.up$x, mart=ensembl)
mapt.up <- bm$hgnc_symbol
mapt.up <- mapt.up[mapt.up %in% genes]
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=mapt.down$x, mart=ensembl)
mapt.down <- bm$hgnc_symbol
mapt.down <- mapt.down[mapt.down %in% genes]

res.grn.up <- bootstrap.enrichment.test(ctd, hits=mapt.up, bg = genes, genelistSpecies = "human", sctSpecies = "human", reps=1000)
res.grn.down <- bootstrap.enrichment.test(ctd, hits=mapt.down, bg = genes, genelistSpecies = "human", sctSpecies = "human", reps=1000)


###
# C9orf72 analysis
###
mapt.up <- read.table("/Users/kevin/dzne/rimod_package/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/DEGs_UP_c9.ndc.txt", sep="\t", header=T)
mapt.down <- read.table("/Users/kevin/dzne/rimod_package/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/DEGs_Down_c9.ndc.txt", sep="\t", header=T)
# transform to HGNC
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=mapt.up$x, mart=ensembl)
mapt.up <- bm$hgnc_symbol
mapt.up <- mapt.up[mapt.up %in% genes]
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=mapt.down$x, mart=ensembl)
mapt.down <- bm$hgnc_symbol
mapt.down <- mapt.down[mapt.down %in% genes]

res.c9.up <- bootstrap.enrichment.test(ctd, hits=mapt.up, bg = genes, genelistSpecies = "human", sctSpecies = "human", reps=1000)
res.c9.down <- bootstrap.enrichment.test(ctd, hits=mapt.down, bg = genes, genelistSpecies = "human", sctSpecies = "human", reps=1000)

