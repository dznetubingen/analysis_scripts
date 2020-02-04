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

setwd("~/rimod/smallRNA/iPSC/")


# Load Count data
counts <- read.table("~/rimod/smallRNA/iPSC/iPSCNeurons_smRNAseq_counts.txt", sep="\t", header=T, row.names = 1, check.names = F)

# Load Metadata
md <- read.table("iPSCNeurons_smRNAseq_metadata_formatted.txt", sep="\t", header=T)
md$group <- factor(md$group)

# format md and counts
md <- md[md$sample %in% colnames(counts),]
counts <- counts[, colnames(counts) %in% md$sample]
md <- md[match(colnames(counts), md$sample),]
rownames(md) <- colnames(counts)

# Make DESeq2 object
dds <- DESeqDataSetFromMatrix(counts,
                              colData = md,
                              design = ~ group)



# Specify control group
dds$dc <- relevel(dds$group, ref = "control")
keep <- rowSums(counts(dds)) >= row_sum_cutoff
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)
resnames <- resultsNames(dds)

#== Extract results ==#
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

# Save results

# Save all results
write.table(res.mapt, "deseq_result_mapt.ndc_iPSCNeurons_smRNAseq.txt", sep="\t", quote=F, col.names = NA)
write.table(res.grn, "deseq_result_grn.ndc_iPSCNeurons_smRNAseq.txt", sep="\t", quote=F, col.names = NA)
write.table(res.c9, "deseq_result_c9.ndc_iPSCNeurons_smRNAseq.txt", sep="\t", quote=F, col.names = NA)

# Save differentially expressed miRNAs according to specified cutoff
write.table(deg.mapt, paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_result_mapt.ndc_iPSCNeurons_smRNAseq.txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(deg.mapt, paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_result_mapt.ndc_iPSCNeurons_smRNAseq.txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(deg.mapt, paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_result_mapt.ndc_iPSCNeurons_smRNAseq.txt", sep=""), sep="\t", quote=F, col.names = NA)


# VST and PCA
vst <- varianceStabilizingTransformation(dds)
plotPCA(vst, intgroup = "group", ntop=100)

# remove batch effect with limma
design <- model.matrix(~ md$dc)
x_noBatch <- removeBatchEffect(rld.mat, batch = md$batch, design=design)
nb <- rld
assay(nb) <- x_noBatch

plotPCA(nb, intgroup = "dc")
png("PCA_sRNA_rimod_frontal_rlog_batchCorrected.png", width=800, height=600)
plotPCA(nb, intgroup = "batch")
dev.off()

# compare genes
mapt_genes <- as.character(rownames(deg.mapt))
grn_genes <- as.character(rownames(deg.grn))
c_genes <- as.character(rownames(deg.c9))

common <- intersect(mapt_genes, intersect(grn_genes, c_genes))
write.table(common, "combined_DE_miRNAs.txt", sep="\t", row.names=F, quote=F)
