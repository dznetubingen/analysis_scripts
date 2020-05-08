library(ggplot2)
library(pheatmap)
setwd("~/../work_dzne/rimod_brain_regions/")

# metadata
md <- read.csv("FTD_Brain_corrected.csv")
md$sample <- paste("sample_", md$GIVENSAMPLENAME, sep="")
dc <- md$DISEASE.CODE
dc[grepl("orad", dc)] <- "FTD-sporadic"
md$DISEASE.CODE <- dc

# cpm data
cpm <- read.table("results_annotation/RiMod_aggrGeneCPM_CAGEseq_all.txt", sep="\t", header=T, row.names=1)

# filter on 246 common samples
cmn_samples <- intersect(md$sample, colnames(cpm))
cpm <- cpm[, colnames(cpm) %in% cmn_samples]
md <- md[md$sample %in% cmn_samples,]
md <- md[match(colnames(cpm), md$sample),]

# Make a PCA
pca <- as.data.frame(prcomp(t(cpm))$x)
pca$region <- md$REGION
pca$group <- md$DISEASE.CODE

p <- ggplot(pca, aes(x=PC1, y=PC2, color=region)) +
  geom_point() +
  theme_minimal()
p


# Make separate PCAs for the different regions
regions <- as.character(levels(md$REGION))

for (reg in regions) {
  print(reg)
  keep <- md$REGION == reg
  reg.md <- md[keep,]
  reg.cpm <- cpm[,keep]
  
  
  pca <- as.data.frame(prcomp(t(reg.cpm))$x)
  pca$group <- reg.md$DISEASE.CODE
  
  p <- ggplot(pca, aes(x=PC1, y=PC2, color=group)) +
    geom_point() +
    theme_minimal() +
    ggtitle(reg)
  print(p)
}

# Cluster for outlier detection
v <- apply(cpm, 1, var)
cpm <- cpm[v > 10,]

pheatmap(cpm, scale="row")
