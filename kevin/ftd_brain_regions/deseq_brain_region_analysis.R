######
# DESeq2 analysis of the Brain Region Data 
######
library(DESeq2)
library(pcaExplorer)
library(ggplot2)

setwd("~/../work_dzne/rimod_brain_regions/")

###
# Data loading and formatting
###
# metadata
md <- read.csv("FTD_Brain_corrected.csv")
md$sample <- paste("sample_", md$GIVENSAMPLENAME, sep="")
dc <- md$DISEASE.CODE
dc[grepl("orad", dc)] <- "FTD-sporadic"
md$DISEASE.CODE <- dc

# cpm data
cpm <- read.table("results_annotation/RiMod_aggrGeneCounts_CAGEseq_all.txt", sep="\t", header=T, row.names=1)

# filter on 246 common samples
cmn_samples <- intersect(md$sample, colnames(cpm))
cpm <- cpm[, colnames(cpm) %in% cmn_samples]
md <- md[md$sample %in% cmn_samples,]
md <- md[match(colnames(cpm), md$sample),]

design <- data.frame(group = md$DISEASE.CODE, gender = md$GENDER, region = md$REGION)
rownames(design) <- colnames(cpm)

# Remove NAs
keep <- !is.na(design$group)
design <- design[keep,]
cpm <- cpm[,keep]

# saver names
design$group <- make.names(design$group)

#===== end data loading =====#


###
# DESeq2 analysis
###
dds <- DESeqDataSetFromMatrix(countData = cpm,
                              colData = design,
                              design = ~ gender + region + group)

# pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# specify control
dds$condition <- relevel(dds$group, ref = "control")

# Run DE analysis
dds <- DESeq(dds)

# Perform VST
dds.vst <- varianceStabilizingTransformation(dds)

# plot PCA
p <- plotPCA(dds.vst, intgroup="region", ntop=1000)
p <- p + theme_minimal()
ggsave("vst_brain_region_pca_ntop1000.png", height=8, width=8)
plotPCA(dds.vst, ntop=1000)
ggsave("vst_brain_regionGroup_pca_ntop1000.png", height=8, width=8)

# Plot PCA for each region
for (reg in as.character(levels(design$region))) {
  print(reg)
  keep <- design$region == reg
  tmp.vsg <- dds.vst[,keep]
  p <- plotPCA(tmp.vsg, intgroup="group", ntop=1000)
  p <- p + theme_minimal() + ggtitle(reg)
  ggsave(paste(reg,"pca_ntop1000.png", sep="_"), height=4, width=4)
}




# Fire up pcaExplorer

pcaExplorer(dds = dds, dst = dds.vst)