setwd("~/rimod/paper_v2/figures/figure5_MMP/")


grn <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/GRN_cell_composition_filterd_DEGs.txt", sep="\t", header=T)

grn <- grn[grn$log2FoldChange > 0,]
grn <- grn[grn$log2FoldChange > 0.5,]

write.table(grn$Row.names, "grn_up_genes_lfc>0.5.txt", quote=F, row.names = F, col.names = F)
