#######################################################
# Analysis of temporal smRNA-seq data of RiMOD project
#######################################################
library(DESeq2)
library(stringr)
library(pheatmap)
library(viridis)
library(limma)

# Parameters
row_sum_cutoff = 10
pval_cutoff = 0.05
lfc_cutoff = 0.2

setwd("~/rimod/smallRNA/temporal/analysis/analysis_0819/")


# Load Count data
counts <- read.table("~/rimod/smallRNA/temporal/rimod_human_temporal_smRNAseq_counts_150120.txt", sep="\t", header=T, row.names = 1, check.names = F)
colnames(counts) <- as.character(str_split(colnames(counts), pattern = "temporal", simplify = T)[,1])

# Load Metadata
md <- read.table("~/rimod/smallRNA/temporal/rimod_human_temporal_smRNAseq_metadata_150120.txt", sep="\t", header=T, check.names=F, row.names = 1)
rownames(md) <- str_pad(md$sample_id, side='left', width = 5, pad = "0")
md$sample_id <- as.factor(md$sample_id)
md$dc <- make.names(md$dc)


# Bring in correct order
rownames(md) <- rownames(md)[match(colnames(counts), rownames(md))]


dds <- DESeqDataSetFromMatrix(counts,
                              colData = md,
                              design = ~ gender + dc)



# Specify control group
dds$dc <- relevel(dds$dc, ref = "control")
keep <- rowSums(counts(dds)) >= row_sum_cutoff
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)
resnames <- resultsNames(dds)

#== Extract results ==#
### MAPT - control
res.mapt <- results(dds, c("dc", "FTD.MAPT", "control"))
res.mapt <- na.omit(res.mapt)
deg.mapt <- res.mapt[res.mapt$padj <= pval_cutoff,]
deg.mapt <- deg.mapt[abs(deg.mapt$log2FoldChange) >= lfc_cutoff,]

### GRN - control
res.grn <- results(dds, c("dc", "FTD.GRN", "control"))
res.grn <- na.omit(res.grn)
deg.grn <- res.grn[res.grn$padj <= pval_cutoff,]
deg.grn <- deg.grn[abs(deg.grn$log2FoldChange) >= lfc_cutoff,]

### C9orf72 - control
res.c9 <- results(dds, c("dc", "FTD.C9", "control"))
res.c9 <- na.omit(res.c9)
deg.c9 <- res.c9[res.c9$padj <= pval_cutoff,]
deg.c9 <- deg.c9[abs(deg.c9$log2FoldChange) >= lfc_cutoff,]

# Save results

# Save all results
write.table(res.mapt, "deseq_result_mapt.ndc_temporal_smRNAseq.txt", sep="\t", quote=F, col.names = NA)
write.table(res.grn, "deseq_result_grn.ndc_temporal_smRNAseq.txt", sep="\t", quote=F, col.names = NA)
write.table(res.c9, "deseq_result_c9.ndc_temporal_smRNAseq.txt", sep="\t", quote=F, col.names = NA)

# Save differentially expressed miRNAs according to specified cutoff
write.table(deg.mapt, paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_result_mapt.ndc_temporal_smRNAseq.txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(deg.grn, paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_result_grn.ndc_temporal_smRNAseq.txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(deg.c9, paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"result_c9.ndc_temporal_smRNAseq.txt", sep=""), sep="\t", quote=F, col.names = NA)

# Save only DEGs (without ohter info) for use in Pathway tools
write.table(rownames(deg.mapt), paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_mapt.ndc_temporal_smRNAseq_miRNAs.txt", sep=""), sep="\t", quote=F, row.names = FALSE)
write.table(rownames(deg.grn), paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_grn.ndc_temporal_smRNAseq_miRNAs.txt", sep=""), sep="\t", quote=F, row.names = FALSE)
write.table(rownames(deg.c9), paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_c9.ndc_temporal_smRNAseq_miRNAs.txt", sep=""), sep="\t", quote=F, row.names = FALSE)



########################################
## Generate count table and rLog table
########################################

# normalized count values
norm.counts <- counts(dds, normalized=TRUE)
write.table(norm.counts, "deseq_normalized_counts_temporal_smRNA.txt", sep="\t", quote=F, col.names = NA)

# reg log transformed values
vst.vals <- varianceStabilizingTransformation(dds, blind=FALSE)
vst.mat <- assay(vst.vals)
write.table(vst.mat, "deseq_rLog_values_temporal_smRNA.txt", sep="\t", quote=F, col.names = NA)


## PCA
pca <- plotPCA(vst.vals, intgroup = "dc")
png("PCA_rimod_frontal_VST_group.png", width=800, height=600)
pca
dev.off()
pca
plotPCA(vst.vals, intgroup = "dc")

plotPCA(vst.vals, intgroup = "age")

# compare genes
mapt_genes <- as.character(rownames(deg.mapt))
grn_genes <- as.character(rownames(deg.grn))
c_genes <- as.character(rownames(deg.c9))

common <- intersect(mapt_genes, intersect(grn_genes, c_genes))
write.table(common, "combined_DE_miRNAs.txt", sep="\t", row.names=F, quote=F)

# make pca
pca <- as.data.frame(prcomp(t(vst.mat))$x)
pca$samples <- rownames(pca)
library(ggplot2)
ggplot(pca, aes(x=PC1, y=PC2)) + 
  geom_text(label=pca$samples)
