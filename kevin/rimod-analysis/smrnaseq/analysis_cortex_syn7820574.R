########
# Analysis of public data from syn7820574
# smRNA-seq
# for validation of our own results
########
library(DESeq2)

# Parameters
row_sum_cutoff = 5
pval_cutoff = 0.05
lfc_cutoff = 0.6


setwd("~/rimod/public_data/miRNA_syn7820574/")

cortex <- read.table("RawCounts_miRDeep2_miRNA_HumanFTD_Cortex.txt", sep="\t", header=T)
cortex <- round(cortex)

# Load and format metadata
md <- read.table("metadata_HumanFTD_Cortex.txt", sep="\t", header=T)
md$dc <- make.names(md$Diagnosis)
rownames(md) <- md$SampleID
md$batch <- factor(md$Batch.RNA.Isolation)


# Make DESeq object
dds <- DESeqDataSetFromMatrix(cortex,
                              colData = md,
                              design = ~ batch + Sex + dc)


# Specify control group
dds$dc <- relevel(dds$dc, ref = "Control")
keep <- rowSums(counts(dds)) >= row_sum_cutoff
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)
resnames <- resultsNames(dds)

#== Extract results ==#
### Tau-positive vs control
res.mapt <- results(dds, c("dc", "FTD.Tau", "Control"))
res.mapt <- na.omit(res.mapt)
deg.mapt <- res.mapt[res.mapt$padj <= pval_cutoff,]

### Tau-negative vs control
res.tauneg <- results(dds, c("dc", "FTD.TauNeg", "Control"))
res.tauneg <- na.omit(res.tauneg)
deg.tauneg <- res.tauneg[res.tauneg$padj <= pval_cutoff,]

# reg log transformed values
rld <- varianceStabilizingTransformation(dds)
rld.mat <- assay(rld)


## PCA
plotPCA(rld, intgroup = "dc")

