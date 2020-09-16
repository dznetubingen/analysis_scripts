#########
# DE analsysis for Rimod Mouse GRN data
########
library(DESeq2)

setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/mouse/grn_mice/")

# Load data
counts <- read.table("RiMod_MouseGRN_aggrGeneCounts_CAGEseq_all.txt", sep="\t", header=T)
md <- read.csv("rimod_grn_mouse_metadata.csv", sep=";")
md$sample <- paste0("sample_GRN_", md$sample_id)

setwd("deseq_analysis/")
# for now, change HET sample to HOM
md[4,3] <- "HO"
md$group <- paste(md$genotype, md$age_bin, sep="_")

# bring in matching order
md <- md[match(colnames(counts), md$sample),]

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

# young mice
res.young <- na.omit(results(dds, c("group", "HO_young", "NTG_young")))
deg.young <- res.young[res.young$padj <= 0.05,]
write.table(rownames(deg.young), "GRN_mice_young_DEGs.txt", quote=F, row.names=F, col.names=F)

# middle mice
res.middle <- na.omit(results(dds, c("group", "HO_middle", "NTG_middle")))
deg.middle <- res.middle[res.middle$padj <= 0.05,]
write.table(rownames(deg.middle), "GRN_mice_middle_DEGs.txt", quote=F, row.names=F, col.names=F)

# old mide
res.old <- na.omit(results(dds, c("group", "HO_old", "NTG_old")))
deg.old <- res.old[res.old$padj <= 0.05,]
write.table(rownames(deg.old), "GRN_mice_old_DEGs.txt", quote=F, row.names=F, col.names=F)


# Save up and down genes
deg.old.up <- deg.old[deg.old$log2FoldChange > 0,]
deg.old.down <- deg.old[deg.old$log2FoldChange < 0,]
write.table(rownames(deg.old.up), "GRN_mice_old_UP.DEGS.txt", quote=F, row.names=F, col.names = F)
write.table(rownames(deg.old.down), "GRN_mice_old_DOWN.DEGS.txt", quote=F, row.names=F, col.names = F)

# Save for STRING-DB
lfc_cutoff <- 0.1
young <- res.young[abs(res.young$log2FoldChange) > lfc_cutoff,]
middle <- res.middle[abs(res.middle$log2FoldChange) > lfc_cutoff,]
old <- res.old[abs(res.old$log2FoldChange) > lfc_cutoff,]
young <- data.frame(gene = rownames(young), lfc = young$log2FoldChange)
middle <- data.frame(gene = rownames(middle), lfc = middle$log2FoldChange)
old <- data.frame(gene = rownames(old), lfc = old$log2FoldChange)
write.table(young, "GRN_young_forString.txt", sep="\t", row.names = F, quote=F)
write.table(middle, "GRN_middle_forString.txt", sep="\t", row.names = F, quote=F)
write.table(old, "GRN_old_forString.txt", sep="\t", row.names = F, quote=F)


# PCA
vst.vals <- varianceStabilizingTransformation(dds, blind=FALSE)
plotPCA(vst.vals, intgroup="group", ntop=500)
