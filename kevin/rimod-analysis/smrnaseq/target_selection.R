###
# Help to select interesting target genes

setwd("~/rimod/smallRNA/frontal/target_selection/")

mirs <- c("hsa-miR-150-5p", "hsa-miR-193a-3p")

# load target edges
mapt <- read.table("../analysis/target_mrna_correlation_analysis_0819/MAPT_miRNA_target_edge_table.txt", sep="\t", header=T)
grn <- read.table("../analysis/target_mrna_correlation_analysis_0819/GRN_miRNA_target_edge_table.txt", sep="\t", header=T)
c9 <- read.table("../analysis/target_mrna_correlation_analysis_0819/C9_miRNA_target_edge_table.txt", sep="\t", header=T)

# mir1
mapt.m1 <- mapt[mapt$mirna == mirs[1],]
grn.m1 <- grn[grn$mirna == mirs[1],]
c9.m1 <- c9[c9$mirna == mirs[1],]
mir1.targets <- intersect(mapt.m1$targets, intersect(grn.m1$targets, c9.m1$targets))


# mir2
mapt.m2 <- mapt[mapt$mirna == mirs[2],]
grn.m2 <- grn[grn$mirna == mirs[2],]
c9.m2 <- c9[c9$mirna == mirs[2],]
mir2.targets <- intersect(mapt.m2$targets, intersect(grn.m2$targets, c9.m2$targets))

# save targets
write.table(mir1.targets, "targets_miR150.txt", sep="\t", quote=F, row.names = F, col.names = F)
write.table(mir2.targets, "targets_miR193.txt", sep="\t", quote=F, row.names = F, col.names = F)


## load DE analysis results
deg.mapt <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/MAPT_cell_composition_filterd_DEGs.txt", sep="\t", header=T)
deg.grn <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/GRN_cell_composition_filterd_DEGs.txt", sep="\t", header=T)
deg.c9 <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/C9orf72_cell_composition_filterd_DEGs.txt", sep="\t", header=T)

table(mir1.targets %in% deg.mapt$hgnc_symbol)
table(mir1.targets %in% deg.grn$hgnc_symbol)

mir2.targets[mir2.targets %in% deg.grn$hgnc_symbol]

test <- mir1.targets[mir1.targets %in% deg.mapt$hgnc_symbol]
write.table(test, "mir1_mapt_degs.txt", sep="\t", quote=F, col.names = F, row.names=F)
