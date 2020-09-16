#########
# DE analsysis for Rimod Mouse MAPT data
########
library(DESeq2)

setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/mouse/mapt_mice/")

# Load data
counts <- read.table("RiMod_MouseMAPT_aggrGeneCounts_CAGEseq_all.txt", sep="\t", header=T)
md <- read.csv("rimod_mouse_mapt_meteadata.csv")
md$sample <- paste0("sample_", md$sample_id)

md <- md[match(colnames(counts), md$sample),]
md$age_bin <- rep("old", nrow(md))
md$age_bin[md$age == 1.5] <- "young"
md$age_bin[md$age == 3] <- "middle"

colnames(md)[4] <- "groupType"
md$group <- paste(md$groupType, md$age_bin, sep="_")

# remove outliers
outliers <- c("sample_16055", "sample_17452", "sample_17450", "sample_15111")
keep <- !md$sample %in% outliers
md <- md[keep,]
counts <- counts[,keep]

# Make DESeq2 object
dds <- DESeqDataSetFromMatrix(counts,
                              colData = md,
                              design = ~ sex  + group)

# Remove unexpressed genes
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)
resnames <- resultsNames(dds)

# Extract results
setwd("deseq_results/")
# young mice
res.young <- na.omit(results(dds, c("group", "Transgenic_young", "Control_young")))
deg.young <- res.young[res.young$padj <= 0.05,]
write.table(rownames(deg.young), "MAPT_mice_young_DEGs.txt", quote=F, row.names=F, col.names=F)

# middle mice
res.middle <- na.omit(results(dds, c("group", "Transgenic_middle", "Control_middle")))
deg.middle <- res.middle[res.middle$padj <= 0.05,]
write.table(rownames(deg.middle), "MAPT_mice_middle_DEGs.txt", quote=F, row.names=F, col.names=F)

# old mide
res.old <- na.omit(results(dds, c("group", "Transgenic_old", "Control_old")))
deg.old <- res.old[res.old$padj <= 0.05,]
write.table(rownames(deg.old), "MAPT_mice_old_DEGs.txt", quote=F, row.names=F, col.names=F)


# PCA
vst.vals <- varianceStabilizingTransformation(dds, blind=FALSE)
plotPCA(vst.vals, intgroup="group")

# outlier identification
library(ggplot2)
pca = as.data.frame(prcomp(t(assay(vst.vals)))$x)
pca$sample <- md$sample
ggplot(pca, aes(x=PC1, y=PC2)) + 
  geom_text(aes(label=sample))
