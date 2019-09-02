####################################
# analysis of hsa-mir-19 effects
####################################
library(stringr)
library(biomaRt)

setwd("~/rimod/integrative_analysis/specifc_analyses/mir19_effects/")

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

####
# MAPT analysis
####

mir.targets <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/MAPT_miRNA_target_edge_table.txt", sep="\t", header=T)
mir19a <- mir.targets[grepl("hsa-miR-19a", mir.targets$mirna),]
mir19b <- mir.targets[grepl("hsa-miR-19b", mir.targets$mirna),]
mir.targets <- rbind(mir19a, mir19b)

write.table(mir.targets$targets, "MAPT_mir19_targets.txt", row.names = F, quote=F)

# Now make table for STRING enrichment analysis
deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-08-12_07.58.35/deseq_result_mapt.ndc_fro_2019-08-12_07.58.35.txt", sep="\t", header=T)
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=deg$X, mart=ensembl)
deg <- merge(deg, bm, by.x="X", by.y="ensembl_gene_id")

deg <- deg[deg$hgnc_symbol %in% mir.targets$targets,]
string <- deg[, c("hgnc_symbol", "log2FoldChange")]
write.table(string, "MAPT_mir19_string_network_file.txt", sep="\t", row.names = F, quote=F)

#==========================================================#

####
# GRN analysis
####

mir.targets <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/GRN_miRNA_target_edge_table.txt", sep="\t", header=T)
mir19a <- mir.targets[grepl("hsa-miR-19a", mir.targets$mirna),]
mir19b <- mir.targets[grepl("hsa-miR-19b", mir.targets$mirna),]
mir.targets <- rbind(mir19a, mir19b)

write.table(mir.targets$targets, "GRN_mir19_targets.txt", row.names = F, quote=F)

#==========================================================#
