###########################
## ESCRT pathway analysis
###########################
library(stringr)
library(biomaRt)
setwd("~/rimod/integrative_analysis/specifc_analyses/escrt_pathway_analysis/")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


# load ESCRT pathway
escrt <- read.table("~/rimod/integrative_analysis/specifc_analyses/escrt_pathway_analysis/reactome_escrt_pathway_molecules.tsv", sep="\t", header=T)
esc.genes <- as.character(escrt$MoleculeName)
esc.genes <- str_split(esc.genes, pattern=" ", simplify = T)[,2]

###
# MAPT analysis
###
mapt <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2019-08-15_14.49.08_frontal/deseq_result_mapt.ndc_2019-08-15_14.49.08.txt",
                   sep="\t", header=T)
mapt <- mapt[mapt$padj <= 0.05,]
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=mapt$X, mart=ensembl)
mapt <- merge(mapt, bm, by.x="X", by.y="ensembl_gene_id")

table(esc.genes %in% mapt$hgnc_symbol)

mapt.esc <- mapt[mapt$hgnc_symbol %in% esc.genes,]

###
# GRN analysis
###
grn <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2019-08-15_14.49.08_frontal/deseq_result_grn.ndc_2019-08-15_14.49.08.txt",
                   sep="\t", header=T)
grn <- grn[grn$padj <= 0.05,]
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=grn$X, mart=ensembl)
grn <- merge(grn, bm, by.x="X", by.y="ensembl_gene_id")

table(esc.genes %in% grn$hgnc_symbol)

grn.esc <- grn[grn$hgnc_symbol %in% esc.genes,]
