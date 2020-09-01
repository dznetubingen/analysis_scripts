#####
# Nice PCA plot for RiMod paper
#####
library(stringr)
library(ggplot2)
setwd("~/rimod/paper/")

###
# RNA-seq
###
mat <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2020-05-04_15.45.57/deseq_vst_values_2020-05-04_15.45.57.txt", sep="\t", header=T, row.names = 1)
colnames(mat) <- gsub("X", "", colnames(mat))
samples <- str_sub( colnames(mat), 1, 5)
samples <- gsub("_", "", samples)
samples <- str_pad(samples, width=5, pad="0", side="left")

# load metadatasamp
md <- read.table("~/rimod/RNAseq/rimod_frontal_rnaseq_metadata_v2.txt", sep="\t", header=T, stringsAsFactors = F)
md <- md[md$SampleID %in% samples,]
md <- md[match(samples, md$SampleID),]

# perform PCA
pca <- prcomp(t(mat))
pca <- as.data.frame(pca$x)
# add metadat
pca$DC <- md$Disease.Code
  
# plotting
p <- ggplot(pca, aes(x=PC1, y=PC2, color=DC, fill=DC)) + 
  geom_point() +
  theme_minimal()
p
ggsave("rnaseq_pca.png", width=4, height = 3)


# cage frontal
mat <- read.table("~/rimod/CAGE/results_annotation/RiMod_aggrGeneCPM_CAGEseq_fro.txt", sep="\t", header=T, row.names=1)
# perform PCA
pca <- prcomp(t(mat))
pca <- as.data.frame(pca$x)
# add metadat
pca$DC <- md$Disease.Code

# plotting
p <- ggplot(pca, aes(x=PC1, y=PC2)) + 
  geom_point() +  theme_minimal() + geom_text(aes(label=rownames(pca)))
p
ggsave("PCA_frontal_cpm.png")



# cage temporal
mat <- read.table("~/rimod/CAGE/results_annotation/RiMod_aggrGeneCPM_CAGEseq_tem.txt", sep="\t", header=T, row.names=1)
# perform PCA
pca <- prcomp(t(mat))
pca <- as.data.frame(pca$x)
# add metadat
pca$DC <- md$Disease.Code

# plotting
p <- ggplot(pca, aes(x=PC1, y=PC2)) + 
  geom_point() +  theme_minimal()+ geom_text(aes(label=rownames(pca)))
p
ggsave("PCA_temporal_cpm.png")
