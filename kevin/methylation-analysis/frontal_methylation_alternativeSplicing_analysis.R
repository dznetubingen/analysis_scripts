###################
# Methylation-Alternative splicing analysis
# check for correlation of methylation with alternatively spliced genes
###
library(biomaRt)
library(stringr)

setwd("~/rimod/Methylation/frontal_methylation_0818/")


# load mart
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Load DMP tables
mapt <- read.table("DMPs_mapt.ndc_quant.txt", sep = "\t", header=T)
grn <- read.table("DMPs_grn.ndc_quant.txt", sep="\t", header=T)
cn <- read.table("DMPs_c9orf72.ndc_quant.txt", sep="\t", header=T)

# Load alternatively spliced genes
as.mapt <- read.table("~/rimod/RNAseq/as_analysis/majiq/grn_AS_genes_dPSI_0.2.txt", sep="\t", header=T)
as.mapt <- as.character(as.mapt$x)
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=as.mapt, mart=ensembl)
as.mapt <- bm$hgnc_symbol

mapt <- grn[grn$adj.P.Val <= 0.05,]
genes <- as.character(mapt$GencodeBasicV12_NAME)
genes <- genes[!genes == ""]
genes <- unlist(str_split(genes, pattern=";"))
genes <- genes[!duplicated(genes)]

