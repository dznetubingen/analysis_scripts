setwd("~/rimod/smallRNA/frontal/")

mapt <- read.table("analysis/target_mrna_correlation_analysis_0819/MAPT_miRNA_target_edge_table.txt", sep="\t", header=T)
grn <- read.table("analysis/target_mrna_correlation_analysis_0819/GRN_miRNA_target_edge_table.txt", sep="\t", header=T)
cn <- read.table("analysis/target_mrna_correlation_analysis_0819/C9_miRNA_target_edge_table.txt", sep="\t", header=T)

mirs <- intersect(mapt$mirna, intersect(grn$mirna, cn$mirna))
# use only up-regulated miRNAs

mir1 <- mirs[1]
mir2 <- mirs[2]
mir3 <- "hsa-miR-19b-3p"


mir <- mir1
tmp1 <- mapt[mapt$mirna == mir,]
tmp2 <- grn[grn$mirna == mir,]
tmp3 <- cn[cn$mirna == mir,]
targets <- union(tmp1$targets, intersect(tmp2$targets, tmp3$targets))
write.table(targets, "targets_miR-150-5p.txt", col.names = F, row.names=F)

mir <- mir2
tmp1 <- mapt[mapt$mirna == mir,]
tmp2 <- grn[grn$mirna == mir,]
tmp3 <- cn[cn$mirna == mir,]
targets <- union(tmp1$targets, intersect(tmp2$targets, tmp3$targets))
write.table(targets, "targets_miR-193a-3p.txt", col.names = F, row.names=F)


tmp1 <- mapt[mapt$mirna == mir3,]
write.table(tmp1$targets, "targets_miR-19b-3p.txt", col.names = F, row.names=F)

