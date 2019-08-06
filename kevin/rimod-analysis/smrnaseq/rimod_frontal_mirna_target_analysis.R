##################################
# Analysis of miRNAs in combination with their predicted targets
# (TargetScan)
###################################
setwd("~/rimod/smallRNA/frontal/analysis/mirna_target_analysis_0719/")


## Load targets
targets <- read.table("miRDB_v6.0_prediction_result.txt", sep="\t")
targets <- targets[grepl("hsa-", targets$V1),]
## Load annotations
anno <- read.table("miR_Family_Info.txt", sep="\t", header=T)

# Load Mapt DEGs
mapt <- read.table("../analysis_0719/DEGs_P0.05_LFC0.8_result_mapt.ndc_frontal_smRNAseq.txt", sep="\t", header=T)
# filter for miRNAs with available targets
mapt <- mapt[mapt$X %in% targets$V1,]
mapt.targets <- targets[targets$V1 %in% mapt$X,]

# Load GRN DEGs
grn <- read.table("../analysis_0719/DEGs_P0.05_LFC0.8_result_grn.ndc_frontal_smRNAseq.txt", sep = "\t", header=T)
grn <- grn[grn$X %in% targets$V1,]
grn.targets <- targets[targets$V1 %in% grn$X,]

# Load C9 DEGs
cn <- read.table("../analysis_0719/DEGs_P0.05_LFC0.8result_c9.ndc_frontal_smRNAseq.txt", sep="\t", header=T)
cn <- cn[cn$X %in% targets$V1,]
cn.target <- targets[targets$V1 %in% cn$X,]

# Save results
write.table(mapt.targets, "MAPT_DEG_targets.txt", sep="\t", quote=F)
write.table(grn.targets, "GRN_DEG_targets.txt", sep="\t", quote=F)
write.table(cn.target, "C9_DEG_targets.txt", sep="\t", quote=F)
