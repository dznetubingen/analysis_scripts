#################
# Rat experiment differential expression analysis
#################
library(DESeq2)
library(stringr)
library(edgeR)
library(ggplot2)
setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rat_data/")

# load count table
counts <- read.table("results/featureCounts/merged_gene_counts.txt", sep="\t", header=T)
rownames(counts) <- counts$Geneid
gene_mapping <- counts[, c(1,2)]
counts <- counts[, c(-1,-2)]

# modify colnames
samples <- colnames(counts)
samples <- gsub("X", "", samples)
samples <- str_split(samples, pattern="_", simplify = T)[,1]
samples <- paste0("sample_", samples)
colnames(counts) <- samples


# load metadata
md <- read.csv("metadata.csv", sep=",", header = F)
md$sample <- paste0("sample_", md$V1)
md <- md[match(samples, md$sample),]
colnames(md) <- c("nr", "group", "gender", "sample")

# make proper group column
groups <- as.character(md$group)
groups <- str_sub(groups, start = 8)
groups <- gsub("_", "", groups)
md$group <- factor(groups)

# Perform DESeq2 Analysis
dds <- DESeqDataSetFromMatrix(counts,
                              colData = md,
                              design = ~ group)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)
resnames <- resultsNames(dds)

# extract results
comp1 <- results(dds, c("group", "KOCerebellumTP1", "WTCerebellumTP1"))
comp2 <- results(dds, c("group", "KOCerebellumTP2", "WTCerebellumTP2"))
comp3 <- results(dds, c("group", "KOBrainStemTP1", "WTBrainStemTP1"))
comp4 <- results(dds, c("group", "KOBrainStemTP2", "WTBrainStemTP2"))

comp1 <- as.data.frame(na.omit(comp1))
comp2 <- as.data.frame(na.omit(comp2))
comp3 <- as.data.frame(na.omit(comp3))
comp4 <- as.data.frame(na.omit(comp4))

write.table(comp1, "Pclo_cerebellum_KOvsWT_TP1.txt", sep = "\t", quote=F, col.names = NA)
write.table(comp2, "Pclo_cerebellum_KOvsWT_TP2.txt", sep = "\t", quote=F, col.names = NA)
write.table(comp3, "Pclo_brainstem_KOvsWT_TP1.txt", sep = "\t", quote=F, col.names = NA)
write.table(comp4, "Pclo_brainstem_KOvsWT_TP2.txt", sep = "\t", quote=F, col.names = NA)

# Get DEGs for profiling
pval_cut = 0.05
deg1 <- comp1[comp1$padj <= pval_cut,]
deg2 <- comp2[comp2$padj <= pval_cut,]
deg3 <- comp3[comp3$padj <= pval_cut,]
deg4 <- comp4[comp4$padj <= pval_cut,]

write.table(rownames(deg1), "DEGs_Pclo_cerebellum_KOvsWT_TP1.txt", col.names = F, row.names = F, quote=F)
write.table(rownames(deg2), "DEGs_Pclo_cerebellum_KOvsWT_TP2.txt", col.names = F, row.names = F, quote=F)
write.table(rownames(deg3), "DEGs_Pclo_brainstem_KOvsWT_TP1.txt", col.names = F, row.names = F, quote=F)
write.table(rownames(deg4), "DEGs_Pclo_brainstem_KOvsWT_TP2.txt", col.names = F, row.names = F, quote=F)

# Divide up and down
# deg2
deg2.up <- deg2[deg2$log2FoldChange > 0,]
deg2.down <- deg2[deg2$log2FoldChange < 0,]
write.table(rownames(deg2.up), "Up_DEGs_Pclo_cerebellum_KOvsWT_TP2.txt", col.names = F, row.names = F, quote=F)
write.table(rownames(deg2.down), "Down_DEGs_Pclo_cerebellum_KOvsWT_TP2.txt", col.names = F, row.names = F, quote=F)

# deg4
deg4.up <- deg4[deg4$log2FoldChange > 0,]
deg4.down <- deg4[deg4$log2FoldChange < 0,]
write.table(rownames(deg4.up), "Up_DEGs_Pclo_brainstem_KOvsWT_TP2.txt", col.names = F, row.names = F, quote=F)
write.table(rownames(deg4.down), "Down_DEGs_Pclo_brainstem_KOvsWT_TP2.txtt", col.names = F, row.names = F, quote=F)




# vst transform
vst.mat <- vst(dds)
plotPCA(ntop=1000, vst.mat, intgroup="group")

plotPCA(ntop=1000, vst.mat, intgroup="gender")


# Test for expression of PCLA
pclo <- gene_mapping[gene_mapping$gene_name == "Pclo",]
cpms <- cpm(counts)
pclo.cpm <- cpms[as.character(pclo$Geneid),]

df <- data.frame(cpm = pclo.cpm)
df$sample <- colnames(counts)
df$group <- groups
df$ko <- str_sub(groups, 1, 2)
df$region <- str_sub(groups, 3, 12)

p <- ggplot(df, aes(x=sample, y=cpm, color=ko)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle=45))
p


# controls only
ctrl <- df[df$ko == "WT",]
p <- ggplot(ctrl, aes(x=region, y=cpm, color=region)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=45))
p

####
# Analysis of TP2 only with Gender as covariate
####
# keep only TP2
keep <- grepl("TP2", md$group)
md.tp2 <- md[keep,]
counts.tp2 <- counts[,keep]

# Perform DESeq2 Analysis
# Include Gender in model design
dds <- DESeqDataSetFromMatrix(counts.tp2,
                              colData = md.tp2,
                              design = ~ gender + group)


keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)
resnames <- resultsNames(dds)

# extract results
comp2 <- results(dds, c("group", "KOCerebellumTP2", "WTCerebellumTP2"))
comp4 <- results(dds, c("group", "KOBrainStemTP2", "WTBrainStemTP2"))

comp2 <- as.data.frame(na.omit(comp2))
comp4 <- as.data.frame(na.omit(comp4))

write.table(comp2, "Pclo_cerebellum_KOvsWT_TP2_GenderCorrected.txt", sep = "\t", quote=F, col.names = NA)
write.table(comp4, "Pclo_brainstem_KOvsWT_TP2_GenderCorrected.txt", sep = "\t", quote=F, col.names = NA)

# Get DEGs for profiling
pval_cut = 0.05
deg2 <- comp2[comp2$padj <= pval_cut,]
deg4 <- comp4[comp4$padj <= pval_cut,]


write.table(rownames(deg2), "DEGs_Pclo_cerebellum_KOvsWT_TP2_GenderCorrected.txt", col.names = F, row.names = F, quote=F)
write.table(rownames(deg4), "DEGs_Pclo_brainstem_KOvsWT_TP2_GenderCorrected.txt", col.names = F, row.names = F, quote=F)

vst.mat <- vst(dds)
plotPCA(vst.mat, intgroup="gender")
plotPCA(vst.mat, intgroup="group")

mat <- assay(vst.mat)
mat <- limma::removeBatchEffect(mat, md.tp2$gender)
assay(vst.mat) <- mat
plotPCA(vst.mat, intgroup = "gender")
