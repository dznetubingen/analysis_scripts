#################
# RiMod miRNA experiment mimic/inhibitor
##############

library(ggplot2)
library(DESeq2)

setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/mirna_mimics/")

cts <- read.csv("results/salmon/salmon_merged_gene_counts.csv", check.names = F, row.names = 1)
md <- read.table("rimod_mirna_neurons_md.txt", sep="\t", header=T)
cts <- round(cts)
dds <- DESeqDataSetFromMatrix(cts,
                              colData = md,
                              design = ~ group)

dds <- estimateSizeFactors(dds)
#keep <- rowSums(counts(dds, normalized=TRUE) >= 10) >= 3
keep <- rowSums(counts(dds, normalized=TRUE)) > 10
dds <- dds[keep,]


# run DESeq2
dds <- DESeq(dds)


#== Extract results ==#
# 19B Mimic
res.19b <- results(dds, c("group", "Mimic19B", "NegMimic"))
res.19b <- na.omit(res.19b)
deg.19b <- res.19b[res.19b$padj <= 0.05,]

deg.19b.up <- deg.19b[deg.19b$log2FoldChange > 0,]
deg.19b.down <- deg.19b[deg.19b$log2FoldChange < 0,]
write.table(rownames(deg.19b.down), "mimic_19B_down.txt", row.names = F, quote=F)
write.table(rownames(deg.19b.up), "mimic_19B_up.txt", row.names = F, quote=F)


# 19B Inhibitor
res.inhib.19b <- results(dds, c("group", "Inhibitor19B", "NegInhibitor"))
res.inhib.19b <- na.omit(res.inhib.19b)
deg.inhib.19b <- res.inhib.19b[res.inhib.19b$padj <= 0.05,]

deg.inhib.19b.up <- deg.inhib.19b[deg.inhib.19b$log2FoldChange > 0,]
deg.inhib.19b.down <- deg.inhib.19b[deg.inhib.19b$log2FoldChange < 0,]
write.table(rownames(deg.inhib.19b.down), "inhib_19B_down.txt", row.names = F, quote=F)
write.table(rownames(deg.inhib.19b.up), "inhib_19B_up.txt", row.names = F, quote=F)


# 150 Mimic
res.150 <- results(dds, c("group", "Mimic150", "NegMimic"))
res.150 <- na.omit(res.150)
deg.150 <- res.150[res.150$padj <= 0.05,]

deg.150.up <- deg.150[deg.150$log2FoldChange > 0,]
deg.150.down <- deg.150[deg.150$log2FoldChange < 0,]
write.table(rownames(deg.150.down), "mimic_150_down.txt", row.names = F, quote=F)
write.table(rownames(deg.150.up), "mimic_150_up.txt", row.names = F, quote=F)


# 150 Inhibitor
res.inhib.150 <- results(dds, c("group", "Inhibitor150", "NegInhibitor"))
res.inhib.150 <- na.omit(res.inhib.150)
deg.inhib.150 <- res.inhib.150[res.inhib.150$padj <= 0.05,]

deg.inhib.150.up <- deg.inhib.150[deg.inhib.150$log2FoldChange > 0,]
deg.inhib.150.down <- deg.inhib.150[deg.inhib.150$log2FoldChange < 0,]
write.table(rownames(deg.inhib.150.down), "inhib_150_down.txt", row.names = F, quote=F)
write.table(rownames(deg.inhib.150.up), "inhib_150_up.txt", row.names = F, quote=F)

vst.mat <- varianceStabilizingTransformation(dds)

plotPCA(vst.mat, intgroup="group")


# testing of targets
targets <- read.table("~/rimod/smallRNA/frontal/target_selection/targets_miR-19b-3p.txt")

test1 <- intersect(targets$V1, rownames(deg.19b.down))
test2 <- intersect(targets$V1, rownames(deg.inhib.19b.up))


mapt <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/MAPT_cell_composition_filterd_DEGs.txt", sep="\t")
mapt.inter <- intersect(mapt$hgnc_symbol, rownames(deg.19b.down))
mapt.inter2 <- intersect(mapt$hgnc_symbol, rownames(deg.inhib.19b.up))

