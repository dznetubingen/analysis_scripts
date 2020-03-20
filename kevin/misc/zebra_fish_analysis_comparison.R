################
# Comparison of Zebra-fish results from Tina to our Rimod-Data results form humans
###############


setwd("~/rimod/public_data/Re__vascularization_phenotype/")

# zebra fish stuff
zeb.up <- read.csv("Suppl. Table1_Up_0.01_0.5_sorted_FC.csv", sep=";")
zeb.down <- read.csv("Suppl. Table2_Down_0.01_0.5_sorted_FC.csv", sep=";")

# GRN genes
grn <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/GRN_cell_composition_filterd_DEGs.txt", sep="\t", header=T)
grn.up <- grn[grn$log2FoldChange > 0,]
grn.down <- grn[grn$log2FoldChange < 0,]

up.ovl <- intersect(zeb.up$GENE, grn.up$hgnc_symbol)
down.ovl <- intersect(zeb.down$GENE, grn.down$hgnc_symbol)
write.table(up.ovl, "up_genes_overlap_GRN.txt", sep="\t", quote=F, row.names = F, col.names = F)
write.table(down.ovl, "down_genes_overlap_GRN.txt", sep="\t", quote=F, row.names = F, col.names = F)


# test for significant overlap of DEGs
# up genes
m <- matrix(c(length(up.ovl), nrow(zeb.up), nrow(grn.up), 20000), nrow=2)
fisher.test(m)

# down genes
m <- matrix(c(length(down.ovl), nrow(zeb.down), nrow(grn.down), 20000), nrow=2)
fisher.test(m)


