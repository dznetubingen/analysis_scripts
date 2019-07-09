########################################
# RiMod ChIP-seq investigative analysis
########################################
library(stringr)
library(DESeq2)
setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/chipseq/")


# GRN
grn <- read.table("results_grn_narrow//bwa/mergedLibrary/macs/narrowPeak/consensus/H3K4me3/deseq2/H3K4me3.consensus_peaks.featureCounts.txt",
                  sep = "\t", header=T)


grn <- grn[, c(-2,-3,-4,-5,-6)]
rownames(grn) <- grn$Geneid
grn <- grn[, -1]
group <- colnames(grn)
group <- gsub("Brain_", "", group)
group <- str_split(group, pattern="_", simplify = T)[,1]
md <- data.frame(group)
rownames(md) <- colnames(grn)

dds <- DESeqDataSetFromMatrix(grn,
                              colData = md,
                              design = ~ group)

dds  <- DESeq(dds)
resnames = resultsNames(dds)

res = results(dds, c("group", "NDC" ,"GRN"))
res <- na.omit(res)
res <- res[res$pvalue <= 0.05,]
res <- res[order(res$log2FoldChange),]

grn <- read.table("results_grn_narrow//bwa/mergedLibrary/macs/narrowPeak/consensus/H3K4me3/deseq2/H3K4me3.consensus_peaks.featureCounts.txt",
sep = "\t", header=T)
grn <- grn[grn$Geneid %in% rownames(res),]
