#########################
# Exploration of sRNA-seq data analyzed with nf-core rnaseq pipeline
########################
library(edgeR)
library(DESeq2)
library(factoextra)

setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/sRNA_frontal/results/")
counts <- read.table("featureCounts/merged_gene_counts.txt", sep="\t", header = TRUE)
rownames(counts) <- counts$Geneid
gene_map <- counts[, c(1,2)] # save mapping for later
counts <- counts[, c(-1, -2)] # remove names and IDs
colnames(counts) <- as.character(sapply(colnames(counts), function(x){strsplit(x, split='Aligned')[[1]][[1]]}))

# Filter file
csum <- rowSums(counts)
counts <- counts[csum > 1,]

counts <- cpm(counts, log=T)


# PCA analysis
pca = prcomp(t(counts), center = TRUE)

fviz_eig(pca)

fviz_pca_ind(pca, col.ind = 'cos2', gradient.cols=c("#00AFBB", "#E7B800", "#FC4E07"), repel=T)
