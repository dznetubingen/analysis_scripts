###########################
# Analysis of repeat expansion results
# using STRique and tandem-genotypes
###########################
library(stringr)
library(ggplot2)

setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/c9orf72_enrichment/analysis_061219/nanopore_enrichmen_101219/")

cutoff = 20

# PBMC22181
pbmc22181 <- read.table("pbmc22181_strique_hg38.tsv", sep="\t", header=T)
pbmc22181 <- pbmc22181[order(pbmc22181$count, decreasing = T),]
# remove small expansions
pbmc22181 <- pbmc22181[pbmc22181$count > cutoff,]
pbmc22181.df <- data.frame(expansion = pbmc22181$count, strand = pbmc22181$strand)
p <- ggplot(pbmc22181.df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "PBMC22181_histogram_STRique_hg38.png")


# PBMC27683
pbmc27683 <- read.table("pbmc27683_strique_hg38.tsv", sep="\t", header=T)
pbmc27683 <- pbmc27683[order(pbmc27683$count, decreasing = T),]
# remove small expansions
pbmc27683 <- pbmc27683[pbmc27683$count > cutoff,]
pbmc27683.df <- data.frame(expansion = pbmc27683$count, strand = pbmc27683$strand)
p <- ggplot(pbmc27683.df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "PBMC27683_histogram_STRique_hg38.png")





