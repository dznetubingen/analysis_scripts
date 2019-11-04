#############################
# Create subsets of exprssion data for YETI
#############################
library(stringr)
setwd("~/rimod/RNAseq/analysis/pathways_analysis/yeti/")

mat <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_vst_values_2019-10-23_13.33.11.txt", sep="\t", header=T, row.names = 1)
rownames(mat) <- str_split(rownames(mat), pattern="[.]", simplify = T)[,1]

# C9orf72
deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_result_c9.ndc_fro_2019-10-23_13.33.11.txt", sep="\t", header=T)
deg <- deg[deg$padj <= 0.05,]
c9.degs <- deg$X
mat.sub <- mat[deg$X,]
write.table(mat.sub, "C9orf72_subMatrix_Yeti.txt", quote=F, col.names = NA, sep="\t")

# MAPT
deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_result_mapt.ndc_fro_2019-10-23_13.33.11.txt", sep="\t", header=T)
deg <- deg[deg$padj <= 0.05,]
mapt.degs <- deg$X
mat.sub <- mat[deg$X,]
write.table(mat.sub, "MAPT_subMatrix_Yeti.txt", quote=F, col.names = NA, sep="\t")

# GRN
deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_result_grn.ndc_fro_2019-10-23_13.33.11.txt", sep="\t", header=T)
deg <- deg[deg$padj <= 0.05,]
mat.sub <- mat[deg$X,]
grn.degs <- deg$X
write.table(mat.sub, "GRN_subMatrix_Yeti.txt", quote=F, col.names = NA, sep="\t")


all.degs <- union(c9.degs, union(mapt.degs, grn.degs))
mat.sub <- mat[all.degs,]
write.table(mat.sub, "FTD_subMatrix_Yeti.txt", quote=F, col.names=NA, sep="\t")
