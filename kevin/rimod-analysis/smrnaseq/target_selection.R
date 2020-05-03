###
# Help to select interesting target genes

setwd("/Users/kevin/dzne/rimod_analysis/target_selection/")

mirs <- c("hsa-miR-150-5p", "hsa-miR-193a-3p", "hsa-miR-19b-3p")

# load target edges
mapt <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/target_mrna_correlation_pval0.1_0420//MAPT_miRNA_target_edge_table.txt", sep="\t", header=T)
grn <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/target_mrna_correlation_pval0.1_0420//GRN_miRNA_target_edge_table.txt", sep="\t", header=T)
c9 <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/target_mrna_correlation_pval0.1_0420//C9_miRNA_target_edge_table.txt", sep="\t", header=T)

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

# mir3
mapt.m3 <- mapt[mapt$mirna == mirs[3],]

# save targets
write.table(mir1.targets, "targets_miR150.txt", sep="\t", quote=F, row.names = F, col.names = F)
write.table(mir2.targets, "targets_miR193.txt", sep="\t", quote=F, row.names = F, col.names = F)
write.table(mapt.m3$targets, "targets_mir19b.txt", sep="\t", quote=F, row.names= F, col.names = F)


## load DE analysis results
deg.mapt <- read.table("/Users/kevin/dzne/rimod_package/analysis/deconvolution/cell_type_specificity/MAPT_cell_composition_filterd_DEGs.txt", sep="\t", header=T)
deg.grn <- read.table("/Users/kevin/dzne/rimod_package/analysis/deconvolution/cell_type_specificity/GRN_cell_composition_filterd_DEGs.txt", sep="\t", header=T)
deg.c9 <- read.table("/Users/kevin/dzne/rimod_package/analysis/deconvolution/cell_type_specificity/C9orf72_cell_composition_filterd_DEGs.txt", sep="\t", header=T)

cmn_degs <- intersect(deg.mapt$hgnc_symbol, deg.grn$hgnc_symbol)
# Get Cell type specificity values
spec <- read.table("/Users/kevin/dzne/rimod_package/analysis/deconvolution/cell_type_specificity/darmanis_cell_type_specificity.txt", sep="\t", header=T, row.names = 1)


# Select targets for mir-150
m.spec <- spec[, colnames(spec) %in% mir1.targets]
m.spec <- as.data.frame(t(m.spec))
m.spec <- m.spec[m.spec$Neurons > 0.2,]
neural_targets <- rownames(m.spec)
neural_targets <- neural_targets[neural_targets %in% cmn_degs]
tmp <- deg.mapt[deg.mapt$hgnc_symbol %in% neural_targets,]
write.table(tmp, "mir150_neura_DE_targets.txt", sep="\t", quote=F, row.names = F)
write.table(tmp$hgnc_symbol, "mir150_only_targets.txt", quote=F, col.names = F, row.names = F)

# Select targets for mir-193
m.spec <- spec[, colnames(spec) %in% mir2.targets]
m.spec <- as.data.frame(t(m.spec))
m.spec <- m.spec[m.spec$Neurons > 0.2,]
neural_targets <- rownames(m.spec)
neural_targets <- neural_targets[neural_targets %in% cmn_degs]
tmp <- deg.mapt[deg.mapt$hgnc_symbol %in% neural_targets,]
write.table(tmp, "mir193_neura_DE_targets.txt", sep="\t", quote=F, row.names = F)
write.table(tmp$hgnc_symbol, "mir193_only_targets.txt", quote=F, col.names = F, row.names = F)

# Select targets for mir-19b
m.spec <- spec[, colnames(spec) %in% mapt.m3$targets]
m.spec <- as.data.frame(t(m.spec))
m.spec <- m.spec[m.spec$Neurons > 0.2,]
neural_targets <- rownames(m.spec)
neural_targets <- neural_targets[neural_targets %in% cmn_degs]
tmp <- deg.mapt[deg.mapt$hgnc_symbol %in% neural_targets,]
write.table(tmp, "mir19b_neura_DE_targets.txt", sep="\t", quote=F, row.names = F)
write.table(tmp$hgnc_symbol, "mir19b_only_targets.txt", quote=F, col.names = F, row.names = F)

