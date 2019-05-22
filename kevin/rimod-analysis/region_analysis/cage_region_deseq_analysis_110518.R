###################################################
## Comparison of frontal and temporal to other
## brain regions
###################################################
library(DESeq2)
library(stringr)

#===================================================#
# Load and format the data
setwd("~/rimod/CAGE/region_analysis/")
md <- read.csv("~/rimod/files/FTD_Brain.csv")
counts <- read.table("cage_all7regions_3kbgr_aggr.txt", sep="\t", row.names=1, header = T)
# Get sample IDs
cols <- colnames(counts)
samples <- as.character(sapply(cols, function(x){strsplit(x, split="_no_")[[1]][[1]]}))
samples <- str_sub(samples, 8, -5)




# Include only controls in analysis
md <- md[md$CASE.CONTROL == "control",]
md <- md[!is.na(md$SAMPLEID),]
md$SAMPLEID <- str_pad(md$SAMPLEID, 5, "left", "0")

# Filter counts for controls
counts_keep <- samples %in% md$SAMPLEID
samples <- samples[counts_keep]
counts <- counts[,counts_keep]

# Extract regions of samples
regs <- as.character(sapply(colnames(counts), function(x){strsplit(x, split="_no_")[[1]][[1]]}))
regs <- str_sub(regs, -3, -1)
regs[regs == "tem"] <- "frotem"
regs[regs == "fro"] <- "frotem"
regs[regs != "frotem"] <- "rest"

cd <- data.frame(region = regs)

#========================================#
# DESeq2 Analysis
#
dds <- DESeqDataSetFromMatrix(counts,
                              colData = cd,
                              design = ~ region)
dds$region <- relevel(dds$region, ref = "rest")

# prefilter
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds)
resnames <- resultsNames(dds)

# Extract results
res <- results(dds, c("region", "frotem", "rest"))
res <- na.omit(res)
deg <- res[res$padj <= 0.0001,]
deg <- deg[abs(deg$log2FoldChange) >= 1,]
write.table(res, "deseq_result_frotem_rest_controls.txt", sep="\t", quote=F, col.names=NA)

# Get norm counts
norm.counts <- counts(dds, normalized=TRUE)
write.table(norm.counts, "controls_norm_counts.txt", sep="\t", quote=F, col.names=NA)
# Rlog transform
rld <- rlog(dds, blind=FALSE)
rld.mat <- assay(rld)
write.table(rld.mat, "controls_rlog_vals.txt", sep="\t", quote=F, col.names=NA)

# Plot PCA
pca <- plotPCA(rld, intgroup="region")
png("pca_frotem_rest_controls.png", height = 900, width=1200)
pca
dev.off()

