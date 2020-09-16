############
# iPSC-neurons mRNA-seq analysis
############
library(tximport)
library(DESeq2)
library(stringr)
library(pheatmap)
library(biomaRt)
library(GenomicFeatures)
library(reshape2)
library(ggplot2)
library(biomaRt)

setwd("~/rimod/Neurons/rimod_neurons/")


###
# Load data as TXI object
samples <- list.files("salmon_quantification")
files <- file.path("salmon_quantification", samples, "quant.sf")
names(files) <- samples
# create txdb
txdb <- makeTxDbFromGFF("gencode.v34.annotation.gff3")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
# load counts
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar=TRUE, countsFromAbundance = "lengthScaledTPM")

counts <- txi$counts
write.table(counts, "frontal_lengthScaeldTPM_counts.txt", sep="\t", quote=F, col.names = NA)

#====== end generating length scaled TPMs ==============#

rownames(counts) <- str_split(rownames(counts), pattern="[.]", simplify = T)[,1]
###########
# Start with analysis
######

md <- read.table("iPSCNeurons_mRNAseq_metadata.txt", sep="\t", header=T, stringsAsFactors = F)
md <- md[md$sample %in% colnames(counts),]
md <- md[match(colnames(counts), md$sample),]

# for testing
#md[14,3] <- "MAPT2"

cts <- round(counts) # round to integer counts
dds <- DESeqDataSetFromMatrix(cts,
                              colData = md,
                              design = ~ flowcell + gender + group)



# Specify control group
dds$DISEASE.CODE <- relevel(dds$group, ref = "control")

# apply prefiltering
row_sum_cutoff = 10
row_sum_samples_nr = 5
dds <- estimateSizeFactors(dds)
keep <- rowSums((counts(dds, normalized=TRUE) >= row_sum_cutoff)) >= row_sum_samples_nr
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds)
resnames <- resultsNames(dds)

#== Extract results ==#
pval_cut <- 0.05
### MAPT - control
res.mapt <- results(dds, c("group", "MAPT", "control"))
res.mapt <- na.omit(res.mapt)
deg.mapt <- res.mapt[res.mapt$padj <= pval_cut,]
print(dim(deg.mapt))

### GRN - control
res.grn <- results(dds, c("group", "GRN", "control"))
res.grn <- na.omit(res.grn)
deg.grn <- res.grn[res.grn$padj <= pval_cut,]
print(dim(deg.grn))

### C9orf72 - control
res.c9 <- results(dds, c("group", "C9", "control"))
res.c9 <- na.omit(res.c9)
deg.c9 <- res.c9[res.c9$padj <= pval_cut,]
print(dim(deg.c9))

###########
## Save results
# Adjust rownames
write.table(res.mapt, "deseq_result_mapt.ndc.txt", sep="\t", quote=F, col.names = NA)
write.table(res.grn, "deseq_result_grn.ndc.txt", sep="\t", quote=F, col.names = NA)
write.table(res.c9, "deseq_result_c9.ndc.txt", sep="\t", quote=F, col.names = NA)


# Save only significant genes for online tools
write.table(rownames(deg.mapt),"DEGs_mapt.ndc.txt", sep="\t", quote=F, row.names=F)
write.table(rownames(deg.grn), "DEGs_grn.ndc.txt", sep="\t", quote=F, row.names=F)
write.table(rownames(deg.c9), "DEGs_c9.ndc.txt", sep="\t", quote=F, row.names=F)

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=rownames(deg.grn), mart=ensembl)
write.table(bm$hgnc_symbol, "DEGs_grn.ndc_HGNC.txt", sep="\t", quote=F, row.names=F)
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=rownames(deg.mapt), mart=ensembl)
write.table(bm$hgnc_symbol, "DEGs_mapt.ndc_HGNC.txt", sep="\t", quote=F, row.names=F)
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=rownames(deg.c9), mart=ensembl)
write.table(bm$hgnc_symbol, "DEGs_c9orf72.ndc_HGNC.txt", sep="\t", quote=F, row.names=F)

# also save for pvale < 0.1
deg.mapt_2 <- res.mapt[res.mapt$padj <= 0.1,]
deg.grn_2 <- res.grn[res.grn$padj <= 0.1,]
deg.c9_2 <- res.c9[res.c9$padj <= 0.1,]

bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=rownames(deg.grn_2), mart=ensembl)
write.table(bm$hgnc_symbol, "DEGs_grn.ndc_HGNC_pval0.1.txt", sep="\t", quote=F, row.names=F)
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=rownames(deg.mapt_2), mart=ensembl)
write.table(bm$hgnc_symbol, "DEGs_mapt.ndc_HGNC_pval0.1.txt", sep="\t", quote=F, row.names=F)
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=rownames(deg.c9_2), mart=ensembl)
write.table(bm$hgnc_symbol, "DEGs_c9orf72.ndc_HGNC_pval0.1.txt", sep="\t", quote=F, row.names=F)


# Divde in up and down regulated genes
# MAPT
mapt.up <- deg.mapt[deg.mapt$log2FoldChange > 0,]
mapt.down <- deg.mapt[deg.mapt$log2FoldChange < 0,]
write.table(rownames(mapt.up), "DEGs_UP_mapt.ndc.txt", quote=F, row.names=F)
write.table(rownames(mapt.down), "DEGs_Down_mapt.ndc.txt", quote=F, row.names=F)
# GRN
grn.up <- deg.grn[deg.grn$log2FoldChange > 0,]
grn.down <- deg.grn[deg.grn$log2FoldChange < 0,]
write.table(rownames(grn.up), "DEGs_UP_grn.ndc.txt", quote=F, row.names=F)
write.table(rownames(grn.down), "DEGs_Down_grn.ndc.txt", quote=F, row.names=F)
# C9orf72
c9.up <- deg.c9[deg.c9$log2FoldChange > 0,]
c9.down <- deg.c9[deg.c9$log2FoldChange < 0,]
write.table(rownames(c9.up), "DEGs_UP_c9.ndc.txt", quote=F, row.names=F)
write.table(rownames(c9.down), "DEGs_Down_c9.ndc.txt", quote=F, row.names=F)

# Filter on fold-change and make tables for STRING-DB
# MAPT
string <- res.mapt[abs(res.mapt$log2FoldChange) > 0.8,]
print(dim(string))
string <- data.frame(gene = rownames(string), lfc = string$log2FoldChange)
write.table(string, "MAPT_string_lfc0.5.txt", sep="\t", row.names=F, quote=F)

# GRN
string <- res.grn[abs(res.grn$log2FoldChange) > 0.8,]
print(dim(string))
string <- data.frame(gene = rownames(string), lfc = string$log2FoldChange)
write.table(string, "GRN_string_lfc0.5.txt", sep="\t", row.names=F, quote=F)

# C9orf72
string <- res.c9[abs(res.c9$log2FoldChange) > 0.8,]
print(dim(string))
string <- data.frame(gene = rownames(string), lfc = string$log2FoldChange)
write.table(string, "C9orf72_string_lfc0.5.txt", sep="\t", row.names=F, quote=F)

###
# PCA plotting and stuff
###
dds.vst <- varianceStabilizingTransformation(dds)
mat <- assay(dds.vst)

write.table(mat, "RiMod_iPSCNeurons_vst_values.txt", sep="\t", quote=F, col.names=NA)

plotPCA(dds.vst, intgroup="mutation", ntop=100)


ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=rownames(res.grn), mart=ensembl)
res.grn <- merge(as.data.frame(res.grn), bm, by.x="row.names", by.y="ensembl_gene_id")
res.grn <- res.grn[order(abs(res.grn$log2FoldChange), decreasing = T),]

# Plot gene expression
plotGene <- function(gene_hgnc){
  gene = bm[bm$hgnc_symbol == gene_hgnc,]$ensembl_gene_id
  e <- melt(mat[gene,])
  e <- merge(e, md, by.x="row.names", by.y="sample")
  ggplot(e, aes(x=group, y=value, fill=group)) + 
    geom_boxplot() +
    geom_point() +
    ggtitle(gene_hgnc)
}

plotGene("CHCHD10")

####
# Plotting for figure
# Make publication ready plots 



