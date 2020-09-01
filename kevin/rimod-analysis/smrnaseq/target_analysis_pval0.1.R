##################################
# Analysis of miRNAs in combination with their predicted targets
# (TargetScan)
###################################
setwd("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/mirna_target_analysis_0719/")


## Load targets
targets <- read.table("miRDB_v6.0_prediction_result.txt", sep="\t")
targets <- targets[grepl("hsa-", targets$V1),]
## Load annotations
anno <- read.table("miR_Family_Info.txt", sep="\t", header=T)

# Load Mapt DEGs
#mapt <- read.table("../analysis_0719/DEGs_P0.05_LFC0.6_result_mapt.ndc_frontal_smRNAseq.txt", sep="\t", header=T)
mapt <- read.table("../analysis_0719/deseq_result_mapt.ndc_frontal_smRNAseq.txt", sep="\t", header=T)
mapt <- mapt[mapt$padj <= 0.1,]
# filter for miRNAs with available targets
mapt <- mapt[mapt$X %in% targets$V1,]
mapt.targets <- targets[targets$V1 %in% mapt$X,]

# Load GRN DEGs
#grn <- read.table("../analysis_0719/DEGs_P0.05_LFC0.6_result_grn.ndc_frontal_smRNAseq.txt", sep = "\t", header=T)
grn <- read.table("../analysis_0719/deseq_result_grn.ndc_frontal_smRNAseq.txt", sep="\t", header=T)
grn <- grn[grn$padj <= 0.1,]
grn <- grn[grn$X %in% targets$V1,]
grn.targets <- targets[targets$V1 %in% grn$X,]

# Load C9 DEGs
#cn <- read.table("../analysis_0719/DEGs_P0.05_LFC0.6result_c9.ndc_frontal_smRNAseq.txt", sep="\t", header=T)
cn <- read.table("../analysis_0719/deseq_result_c9.ndc_frontal_smRNAseq.txt", sep="\t", header=T)
cn <- cn[cn$padj <= 0.1,]
cn <- cn[cn$X %in% targets$V1,]
cn.target <- targets[targets$V1 %in% cn$X,]

# Save results
write.table(mapt.targets, "MAPT_DEG_targets_pval0.1.txt", sep="\t", quote=F)
write.table(grn.targets, "GRN_DEG_targets_pval0.1.txt", sep="\t", quote=F)
write.table(cn.target, "C9_DEG_targets_pval0.1.txt", sep="\t", quote=F)
