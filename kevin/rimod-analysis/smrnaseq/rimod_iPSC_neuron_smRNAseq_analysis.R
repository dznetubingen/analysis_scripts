#######################################################
# Analysis of fronal smRNA-seq data of RiMOD project
#######################################################
library(DESeq2)
library(stringr)
library(pheatmap)
library(viridis)
library(limma)

# Parameters
row_sum_cutoff = 5
pval_cutoff = 0.05
lfc_cutoff = 0.6

setwd("~/rimod/smallRNA/iPSC/analysis/")


# Load Count data
counts <- read.table("~/rimod/smallRNA/iPSC/iPSCNeurons_smRNAseq_counts.txt", sep="\t", header=T, row.names = 1, check.names = F)

# Load Metadata
md <- read.table("~/rimod/smallRNA/iPSC/iPSCNeurons_smRNAseq_metadata_formatted_2.txt", sep="\t", header=T)
md$group <- factor(md$group)

# format md and counts
md <- md[md$sample %in% colnames(counts),]
counts <- counts[, colnames(counts) %in% md$sample]
md <- md[match(colnames(counts), md$sample),]
rownames(md) <- colnames(counts)

# remove outliers
#keep <- !md$sample %in% c("F_FD_07Cl23", "LZ_FD_13Cl10")
#md <- md[keep,]
#counts <- counts[,keep]

# Make DESeq2 object
dds <- DESeqDataSetFromMatrix(counts,
                              colData = md,
                              design = ~ gender + batchN + group)



# Specify control group
dds$dc <- relevel(dds$group, ref = "control")
keep <- rowSums(counts(dds)) >= row_sum_cutoff
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)
resnames <- resultsNames(dds)

#== Extract results ==#
# NOTE: use independent hypothesis weighting for filter
# Igantiadis et al., Nature Methods, 2016
### MAPT - control
res.mapt <- results(dds, c("group", "MAPT", "control"))
res.mapt <- na.omit(res.mapt)
deg.mapt <- res.mapt[res.mapt$padj <= pval_cutoff,]
#deg.mapt <- deg.mapt[abs(deg.mapt$log2FoldChange) >= lfc_cutoff,]
print(deg.mapt)

### GRN - control
res.grn <- results(dds, c("group", "GRN", "control"))
res.grn <- na.omit(res.grn)
deg.grn <- res.grn[res.grn$padj <= pval_cutoff,]
#deg.grn <- deg.grn[abs(deg.grn$log2FoldChange) >= lfc_cutoff,]
print(deg.grn)

### C9orf72 - control
res.c9 <- results(dds, c("group", "C9", "control"))
res.c9 <- na.omit(res.c9)
deg.c9 <- res.c9[res.c9$padj <= pval_cutoff,]
#deg.c9 <- deg.c9[abs(deg.c9$log2FoldChange) >= lfc_cutoff,]
print(deg.c9)

# Save all results
write.table(res.mapt, "deseq_result_mapt.ndc_iPSCNeurons_smRNAseq.txt", sep="\t", quote=F, col.names = NA)
write.table(res.grn, "deseq_result_grn.ndc_iPSCNeurons_smRNAseq.txt", sep="\t", quote=F, col.names = NA)
write.table(res.c9, "deseq_result_c9.ndc_iPSCNeurons_smRNAseq.txt", sep="\t", quote=F, col.names = NA)

# Save differentially expressed miRNAs according to specified cutoff
write.table(deg.mapt, paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_result_mapt.ndc_iPSCNeurons_smRNAseq.txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(deg.grn, paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_result_grn.ndc_iPSCNeurons_smRNAseq.txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(deg.c9, paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"result_c9.ndc_iPSCNeurons_smRNAseq.txt", sep=""), sep="\t", quote=F, col.names = NA)

# Save only DEGs (without ohter info) for use in Pathway tools
write.table(rownames(deg.mapt), paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_mapt.ndc_iPSCNeurons_smRNAseq_miRNAs.txt", sep=""), sep="\t", quote=F, row.names = FALSE)
write.table(rownames(deg.grn), paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_grn.ndc_iPSCNeurons_smRNAseq_miRNAs.txt", sep=""), sep="\t", quote=F, row.names = FALSE)
write.table(rownames(deg.c9), paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_c9.ndc_iPSCNeurons_smRNAseq_miRNAs.txt", sep=""), sep="\t", quote=F, row.names = FALSE)



########################################
## Generate count table and rLog table
########################################

# normalized count values
norm.counts <- counts(dds, normalized=TRUE)
write.table(norm.counts, "deseq_normalized_counts_iPSCNeurons_smRNA.txt", sep="\t", quote=F, col.names = NA)

# reg log transformed values
vst.vals <- varianceStabilizingTransformation(dds, blind=FALSE)
vst.mat <- assay(vst.vals)
write.table(vst.mat, "deseq_vst_values_iPSCNeurons_smRNA.txt", sep="\t", quote=F, col.names = NA)


## PCA
pca <- plotPCA(vst.vals, intgroup = "group")
png("PCA_rimod_frontal_VST_group.png", width=800, height=600)
pca
dev.off()

plotPCA(vst.vals, intgroup = "batchN")
plotPCA(vst.vals, intgroup = "batchS")


library(ggplot2)

pca <- prcomp(t(assay(vst.vals)))
pca.df <- as.data.frame(pca$x)
pca.df$sample <- rownames(pca.df)

ggplot(pca.df, aes(x=PC1, y=PC2)) + 
  geom_text(label=pca.df$sample)


mat <- assay(vst.vals)
library(limma)
design <- model.matrix(~md$group)
test <- removeBatchEffect(mat, batch=md$batchN, design=design)

pca.df <- as.data.frame(prcomp(t(test))$x)
pca.df$sample <- rownames(pca.df)
pca.df$batchN <- md$batchN
pca.df$batchS <- md$batchS
pca.df$group <- md$group

ggplot(pca.df, aes(x=PC1, y=PC2, color=group)) + 
  geom_text(label=pca.df$sample)


# plot expression
gene = "hsa-miR-142-3p"
g <- melt(vst.mat[gene,])
g <- merge(g, md, by.x="row.names", by.y="sample")
ggplot(g, aes(x=group, y=value, fill=group)) + 
  geom_boxplot() + geom_point()
