###########################
# Analysis of repeat expansion results
# using STRique and tandem-genotypes
###########################
library(stringr)
library(ggplot2)

setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/c9_enrichment/analysis_06022020/")

cutoff = 10

###
# Sample 73
###
strique <- read.table("c9orf72_rotterdam73.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
# remove small expansions
strique <- strique[strique$count > cutoff,]
df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "Rotterdam_73_histogram_STRique_hg38.png")

# no cutoff
strique <- read.table("c9orf72_rotterdam73.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "Rotterdam_73_histogram_STRique_hg38_noCutoff.png")


#=======================#


###
# Sample 74
###
strique <- read.table("c9orf72_rotterdam74.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
# remove small expansions
strique <- strique[strique$count > cutoff,]
df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "Rotterdam_74_histogram_STRique_hg38.png")

# no cutoff
strique <- read.table("c9orf72_rotterdam74.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "Rotterdam_74_histogram_STRique_hg38_noCutoff.png")


#=======================#
