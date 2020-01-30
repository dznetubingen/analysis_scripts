###########################
# Analysis of repeat expansion results
# using STRique and tandem-genotypes
###########################
library(stringr)
library(ggplot2)

setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/c9_enrichment/strique_results/")

cutoff = 10

###
# Sample 29254
###
strique <- read.table("c9orf72_29254.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
# remove small expansions
strique <- strique[strique$count > cutoff,]
df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "Sample_29254_histogram_STRique_hg38.png")

#=======================#


###
# Sample 29915
###
strique <- read.table("c9orf72_29915.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
# remove small expansions
strique <- strique[strique$count > cutoff,]
df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "Sample_29915_histogram_STRique_hg38.png")

#=======================#


###
# Sample 30070
###
strique <- read.table("c9orf72_30070.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
# remove small expansions
strique <- strique[strique$count > cutoff,]
df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "Sample_30070_histogram_STRique_hg38.png")

#=======================#


###
# Sample 30929
###
strique <- read.table("c9orf72_30929.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
# remove small expansions
strique <- strique[strique$count > cutoff,]
df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "Sample_30929_histogram_STRique_hg38.png")

#=======================#
