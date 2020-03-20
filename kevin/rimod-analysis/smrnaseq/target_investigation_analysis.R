library(biomaRt)

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

ats <- listAttributes(ensembl)
ats <- ats[,1]
ats[grepl("refs", ats)]


setwd("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819")

# MAPT upMir targets
mapt <- read.table("MAPT_upMir_correlated_targets_Refseq.txt", sep="\t", header=T)
mapt <- as.character(mapt$x)
mapt <- mapt[!duplicated(mapt)]

bm <- getBM(attributes = c("hgnc_symbol", "refseq_mrna"), filters = "refseq_mrna", values=mapt, mart=ensembl)
mapt <- bm$hgnc_symbol
mapt <- mapt[!duplicated(mapt)]


deg <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/MAPT_cell_composition_filterd_DEGs.txt", sep="\t", header=T)

mapt <- mapt[mapt %in% deg$hgnc_symbol]
write.table(mapt, "~/tmp_stuff/mapt_test.txt", row.names = F, col.names = F, quote=F)

# grn up mir targets
grn <- read.table("C9_upMir_correlated_targets_Refseq.txt", sep="\t", header=T)
grn <- as.character(grn$x)
grn <- grn[!duplicated(grn)]
bm <- getBM(attributes = c("hgnc_symbol", "refseq_mrna"), filters = "refseq_mrna", values=grn, mart=ensembl)
grn <- bm$hgnc_symbol
grn <- grn[!duplicated(grn)]

#deg <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/C", sep="\t", header=T)

#grn <- grn[grn %in% deg$hgnc_symbol]
write.table(grn, "~/tmp_stuff/c9orf72_test.txt", row.names = F, col.names = F, quote=F)






