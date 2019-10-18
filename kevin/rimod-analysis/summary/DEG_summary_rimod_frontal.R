###########################
# RiMod DE summary (frontal)
# Collect all DEGs from RNA-seq and CAGE-seq, count them, and compare between the two datatypes
###########################
library(UpSetR)
library(VennDiagram)
setwd("~/rimod/summary/RiMod_DEGs/")

# Cutoffs
pval <- 0.05
lfc <- 0.6
filterLFC <- TRUE
###
# MAPT
###
mapt.rna <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-08-12_07.58.35/deseq_result_mapt.ndc_fro_2019-08-12_07.58.35.txt",
                       sep="\t", header=T, row.names=1)
mapt.cage <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2019-08-15_14.49.08_frontal/deseq_result_mapt.ndc_2019-08-15_14.49.08.txt",
                        sep="\t", header=T, row.names=1)
mapt.rna.ct <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/MAPT_cell_composition_filterd_DEGs.txt", sep="\t", header=T, row.names=2)
mapt.cage.ct <- read.table("~/rimod/CAGE/cage_analysis/deconvolution_frontal/cell_specificity_filtering/MAPT_cell_composition_filterd_DEGs.txt",
                           sep="\t", header=T, row.names=2)
# filter with pvalue
mapt.rna <- mapt.rna[mapt.rna$padj <= pval,]
mapt.cage <- mapt.cage[mapt.cage$padj <= pval,]
mapt.rna.ct <- mapt.rna.ct[mapt.rna.ct$padj <= pval,]
mapt.cage.ct <- mapt.cage.ct[mapt.cage.ct$padj <= pval,]
# filter with lfc
if (filterLFC){
  mapt.rna <- mapt.rna[abs(mapt.rna$log2FoldChange) >= lfc,]
  mapt.cage <- mapt.cage[abs(mapt.cage$log2FoldChange) >= lfc,]
  mapt.rna.ct <- mapt.rna.ct[abs(mapt.rna.ct$log2FoldChange) >= lfc,]
  mapt.cage.ct <- mapt.cage.ct[abs(mapt.cage.ct$log2FoldChange) >= lfc,]
}


mapt.list <- list(rownames(mapt.rna),
                  rownames(mapt.cage),
                  rownames(mapt.rna.ct),
                  rownames(mapt.cage.ct))
names(mapt.list) <- c("MAPT.RNA", "MAPT.CAGE", "MAPT.RNA.FLT", "MAPT.CAGE.FLT")

upset(fromList(mapt.list), order.by = "freq")


###
# GRN
###
grn.rna <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-08-12_07.58.35/deseq_result_grn.ndc_fro_2019-08-12_07.58.35.txt",
                       sep="\t", header=T, row.names=1)
grn.cage <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2019-08-15_14.49.08_frontal/deseq_result_grn.ndc_2019-08-15_14.49.08.txt",
                        sep="\t", header=T, row.names=1)
grn.rna.ct <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/GRN_cell_composition_filterd_DEGs.txt", sep="\t", header=T, row.names=1)
grn.rna.ct <- grn.rna.ct[!duplicated(grn.rna.ct$Row.names),]
rownames(grn.rna.ct) <- grn.rna.ct$Row.names
grn.cage.ct <- read.table("~/rimod/CAGE/cage_analysis/deconvolution_frontal/cell_specificity_filtering/GRN_cell_composition_filterd_DEGs.txt",
                           sep="\t", header=T, row.names=2)

# filter with pvalue
grn.rna <- grn.rna[grn.rna$padj <= pval,]
grn.cage <- grn.cage[grn.cage$padj <= pval,]
grn.rna.ct <- grn.rna.ct[grn.rna.ct$padj <= pval,]
grn.cage.ct <- grn.cage.ct[grn.cage.ct$padj <= pval,]

# filter with lfc
if (filterLFC){
  grn.rna <- grn.rna[abs(grn.rna$log2FoldChange) >= lfc,]
  grn.cage <- grn.cage[abs(grn.cage$log2FoldChange) >= lfc,]
  grn.rna.ct <- grn.rna.ct[abs(grn.rna.ct$log2FoldChange) >= lfc,]
  grn.cage.ct <- grn.cage.ct[abs(grn.cage.ct$log2FoldChange) >= lfc,]
}


grn.list <- list(rownames(grn.rna),
                  rownames(grn.cage),
                  rownames(grn.rna.ct),
                  rownames(grn.cage.ct))
names(grn.list) <- c("GRN.RNA", "GRN.CAGE", "GRN.RNA.FLT", "GRN.CAGE.FLT")

upset(fromList(grn.list), order.by = "freq")


###
# C9orf72
###
c9.rna <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-08-12_07.58.35/deseq_result_c9.ndc_fro_2019-08-12_07.58.35.txt",
                      sep="\t", header=T, row.names=1)
c9.cage <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2019-08-15_14.49.08_frontal/deseq_result_c9.ndc_2019-08-15_14.49.08.txt",
                       sep="\t", header=T, row.names=1)
c9.rna.ct <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/C9orf72_cell_composition_filterd_DEGs.txt", sep="\t", header=T, row.names=1)
c9.cage.ct <- read.table("~/rimod/CAGE/cage_analysis/deconvolution_frontal/cell_specificity_filtering/C9orf72_cell_composition_filterd_DEGs.txt",
                          sep="\t", header=T, row.names=2)

# filter with pvalue
c9.rna <- c9.rna[c9.rna$padj <= pval,]
c9.cage <- c9.cage[c9.cage$padj <= pval,]
c9.rna.ct <- c9.rna.ct[c9.rna.ct$padj <= pval,]
c9.cage.ct <- c9.cage.ct[c9.cage.ct$padj <= pval,]
# filter with lfc
if (filterLFC){
  c9.rna <- c9.rna[abs(c9.rna$log2FoldChange) >= lfc,]
  c9.cage <- c9.cage[abs(c9.cage$log2FoldChange) >= lfc,]
  c9.rna.ct <- c9.rna.ct[abs(c9.rna.ct$log2FoldChange) >= lfc,]
  c9.cage.ct <- c9.cage.ct[abs(c9.cage.ct$log2FoldChange) >= lfc,]
}

c9.list <- list(rownames(c9.rna),
                 rownames(c9.cage),
                 rownames(c9.rna.ct),
                 rownames(c9.cage.ct))
names(c9.list) <- c("C9orf72.RNA", "C9orf72.CAGE", "C9orf72.RNA.FLT", "C9orf72.CAGE.FLT")

upset(fromList(c9.list), order.by = "freq")





###
# Check overlap between comparisons

# RNA-seq
rna.list <- list(rownames(mapt.rna),
                 rownames(grn.rna),
                 rownames(c9.rna))
names(rna.list) <- c("MAPT", "GRN", "C9orf72")
upset(fromList(rna.list), order.by = "freq")

# CAGE-seq
cage.list <- list("MAPT"=rownames(mapt.cage),
                  "GRN"=rownames(grn.cage),
                  "C9orf72"=rownames(c9.cage))
upset(fromList(cage.list), order.by = "freq")










##############
# miRNA DEGs comparison
##############

mir.mapt <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_result_mapt.ndc_frontal_smRNAseq.txt",
                       sep = "\t", header=T, row.names = 1)
mir.grn <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_result_grn.ndc_frontal_smRNAseq.txt",
                      sep="\t", header=T, row.names=1)
mir.c9 <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_result_c9.ndc_frontal_smRNAseq.txt",
                     sep="\t", header=T, row.names=1)

# Filtering
mir.mapt <- mir.mapt[mir.mapt$padj <= pval,]
mir.grn <- mir.grn[mir.grn$padj <= pval,]
mir.c9 <- mir.c9[mir.c9$padj <= pval,]
if (filterLFC){
  mir.mapt <- mir.mapt[abs(mir.mapt$log2FoldChange) >= lfc,]
  mir.grn <- mir.grn[abs(mir.grn$log2FoldChange) >= lfc,]
  mir.c9 <- mir.c9[abs(mir.c9$log2FoldChange) >= lfc,]
}

# upset plot
mir.list <- list("MAPT"=rownames(mir.mapt),
                 "GRN"=rownames(mir.grn),
                 "C9orf72"=rownames(mir.c9))
upset(fromList(mir.list), order.by="freq")
