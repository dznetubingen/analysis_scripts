############
# Look enrichment of activate TF targets in module genes
############
library(biomaRt)

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

setwd("~/rimod/RNAseq/analysis/human_base/")


tfs <- read.table("~/rimod/CAGE/cage_analysis/tf_enrichment_analysis_050420/GRN_up_targets.txt", sep="\t", header=T)
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=tfs$TargetGene, mart=ensembl)
tfs <- merge(tfs, bm, by.x="TargetGene", by.y="ensembl_gene_id")


# load them modules
m1 <- read.table("module_genes/GRN_M1_up_genes.txt")
m3 <- read.table("module_genes/GRN_M3_up_genes.txt")
m4 <- read.table("module_genes/GRN_M4_up_genes.txt")
m6 <- read.table("module_genes/GRN_M6_up_genes.txt")

m1.tf <- tfs[tfs$hgnc_symbol %in% m1$V1,]
m3.tf <- tfs[tfs$hgnc_symbol %in% m3$V1,]
m4.tf <- tfs[tfs$hgnc_symbol %in% m4$V1,]
m6.tf <- tfs[tfs$hgnc_symbol %in% m6$V1,]

tf.m1 <- tfs[tfs$TF %in% m1$V1,]
tf.m1 <- tf.m1[!duplicated(tf.m1$TF),]
tf.m3 <- tfs[tfs$TF %in% m3$V1,]
tf.m3 <- tf.m3[!duplicated(tf.m3$TF),]
tf.m4 <- tfs[tfs$TF %in% m4$V1,]
tf.m4 <- tf.m4[!duplicated(tf.m4$TF),]
tf.m6 <- tfs[tfs$TF %in% m6$V1,]
tf.m6 <- tf.m6[!duplicated(tf.m6$TF),]
