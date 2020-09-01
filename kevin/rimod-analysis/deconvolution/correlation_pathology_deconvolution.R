# Correlation analysis of pathological report with the deconvolution data
# for the RiMod RNA-seq data of the frontal lobe
library(stringr)
library(ggplot2)

setwd("~/rimod/")

# load and format pathology data
deg <- read.csv("Paraffin_slides_to_be_checked_MN.csv", sep=";")
deg <- deg[, c(1,2)]
deg <- deg[!is.na(deg$Deg),]
deg$sample <- gsub("S", "", gsub("/", "", deg$sample))


# load deconvolution data
res <- read.table("~/rimod/RNAseq/analysis/deconvolution/cdn_predictions.txt", sep="\t", header=T)
# format deconvolution data
sample <- res$X
sample <- gsub("X", "", sample)
sample <- str_split(sample, pattern="_", simplify = T)[,1]
res$X <- sample

# only use common samples
cmn <- intersect(res$X, deg$sample)

res <- res[res$X %in% cmn,]
deg <- deg[deg$sample %in% cmn,]

deg <- deg[match(res$X, deg$sample),]

# Perform correlation
cor.test(deg$Deg, res$ExNeurons)


plot(deg$Deg, res$ExNeurons)
abline(lm(res$ExNeurons ~ deg$Deg))

df <- data.frame(Score = deg$Deg, ExNeurons = res$ExNeurons)

ggplot(df, aes(x=Score, y=ExNeurons)) + 
  geom_point() + 
  geom_smooth(method = lm, se=F) +
  theme_minimal()
ggsave("~/rimod/paper_v2/figures/figure3/pathology_correlation.png", width=3, height = 3)
