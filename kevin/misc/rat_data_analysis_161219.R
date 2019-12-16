#################
# Rat experiment differential expression analysis
#################
library(DESeq2)
library(stringr)
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
md <- read.csv("metadata.csv", sep=";", header = F)
md$sample <- paste0("sample_", md$V1)
md <- md[match(samples, md$sample),]
colnames(md) <- c("nr", "group", "sample")

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


# vst transform
vst.mat <- vst(dds)
plotPCA(vst.mat, intgroup="group")

# Test for expression of PCLA
pclo <- gene_mapping[gene_mapping$gene_name == "Pclo",]
vst.mat <- assay(vst.mat)
cpms <- cpm(counts)
pclo.cpm <- cpms[as.character(pclo$Geneid),]

library(ggplot2)
df <- data.frame(cpm = pclo.cpm)
df$sample <- colnames(counts)
df$group <- groups
df$ko <- str_sub(groups, 1, 2)

p <- ggplot(df, aes(x=sample, y=cpm, color=ko)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle=45))
p



