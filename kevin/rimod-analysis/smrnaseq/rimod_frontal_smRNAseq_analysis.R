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

setwd("~/rimod/smallRNA/frontal/analysis/analysis_0719/")


# Load Count data
counts <- read.table("~/rimod/smallRNA/frontal/rimod_human_frontal_smRNAseq_counts.txt", sep="\t", header=T, row.names = 1, check.names = F)
# only keep human miRNAs
keep <- grepl("hsa-", rownames(counts))
counts <- counts[keep,]

# Load Metadata
md <- read.table("~/rimod/smallRNA/frontal/rimod_human_frontal_smRNAseq_metadata.txt", sep="\t", header=T, check.names=F, row.names = 1)
md$id <- as.factor(md$id)

# PH
ph <- read.csv("~/rimod/files/FTD_Brain.csv", stringsAsFactors = F)
ph <- ph[ph$REGION == "frontal",]
ph <- ph[ph$SAMPLEID %in% md$id,]
ph <- ph[match(md$id, ph$SAMPLEID),]
md$ph <- as.numeric(ph$PH)
ph.mean <- mean(na.omit(md$ph))
md$ph[is.na(md$ph)] <- ph.mean


dds <- DESeqDataSetFromMatrix(counts,
                              colData = md,
                              design = ~ ph + batch +  gender + dc)



# Specify control group
dds$dc <- relevel(dds$dc, ref = "NDC")
keep <- rowSums(counts(dds)) >= row_sum_cutoff
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)
resnames <- resultsNames(dds)

#== Extract results ==#
### MAPT - control
res.mapt <- results(dds, c("dc", "FTD.MAPT", "NDC"))
res.mapt <- na.omit(res.mapt)
deg.mapt <- res.mapt[res.mapt$padj <= pval_cutoff,]
deg.mapt <- deg.mapt[abs(deg.mapt$log2FoldChange) >= lfc_cutoff,]

### GRN - control
res.grn <- results(dds, c("dc", "FTD.GRN", "NDC"))
res.grn <- na.omit(res.grn)
deg.grn <- res.grn[res.grn$padj <= pval_cutoff,]
deg.grn <- deg.grn[abs(deg.grn$log2FoldChange) >= lfc_cutoff,]

### C9orf72 - control
res.c9 <- results(dds, c("dc", "FTD.C9", "NDC"))
res.c9 <- na.omit(res.c9)
deg.c9 <- res.c9[res.c9$padj <= pval_cutoff,]
deg.c9 <- deg.c9[abs(deg.c9$log2FoldChange) >= lfc_cutoff,]


# Save all results
write.table(res.mapt, "deseq_result_mapt.ndc_frontal_smRNAseq.txt", sep="\t", quote=F, col.names = NA)
write.table(res.grn, "deseq_result_grn.ndc_frontal_smRNAseq.txt", sep="\t", quote=F, col.names = NA)
write.table(res.c9, "deseq_result_c9.ndc_frontal_smRNAseq.txt", sep="\t", quote=F, col.names = NA)

# Save differentially expressed miRNAs according to specified cutoff
write.table(deg.mapt, paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_result_mapt.ndc_frontal_smRNAseq.txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(deg.grn, paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_result_grn.ndc_frontal_smRNAseq.txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(deg.c9, paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"result_c9.ndc_frontal_smRNAseq.txt", sep=""), sep="\t", quote=F, col.names = NA)

# Save only DEGs (without ohter info) for use in Pathway tools
write.table(rownames(deg.mapt), paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_mapt.ndc_frontal_smRNAseq_miRNAs.txt", sep=""), sep="\t", quote=F, row.names = FALSE)
write.table(rownames(deg.grn), paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_grn.ndc_frontal_smRNAseq_miRNAs.txt", sep=""), sep="\t", quote=F, row.names = FALSE)
write.table(rownames(deg.c9), paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_c9.ndc_frontal_smRNAseq_miRNAs.txt", sep=""), sep="\t", quote=F, row.names = FALSE)



########################################
## Generate count table and rLog table
########################################

# normalized count values
norm.counts <- counts(dds, normalized=TRUE)
write.table(norm.counts, "deseq_normalized_counts_temporal_smRNA.txt", sep="\t", quote=F, col.names = NA)

# reg log transformed values
vst.vals <- varianceStabilizingTransformation(dds, blind=FALSE)
vst.mat <- assay(vst.vals)
write.table(vst.mat, "deseq_vst_values_frontal_smRNA.txt", sep="\t", quote=F, col.names = NA)


## PCA
pca <- plotPCA(vst.vals, intgroup = "dc")
png("PCA_rimod_frontal_VST_group.png", width=800, height=600)
pca
dev.off()
pca
plotPCA(rld, intgroup = "dc")

# remove batch effect with limma
design <- model.matrix(~ md$dc)
x_noBatch <- removeBatchEffect(vst.mat, batch = md$batch, design=design)
nb <- vst.vals
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
