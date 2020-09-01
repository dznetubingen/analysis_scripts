###########################
# Analysis of repeat expansion results
# using STRique and tandem-genotypes
###########################
library(stringr)
library(ggplot2)

setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/c9_enrichment/nanopre_cas9nanoprotocol_110520")

cutoff = 10

###
# DN19 C1
###
strique <- read.table("c9orf72_DN19_C1.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
# remove small expansions
strique <- strique[strique$count > cutoff,]
df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "DN19_C1_histogram_STRique_hg38.png")

# no cutoff
strique <- read.table("c9orf72_DN19_C1.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "DN19_C1_histogram_STRique_hg38_noCutoff.png")


#=======================#


###
# DN19 C2
###
strique <- read.table("c9orf72_DN19_C2.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
# remove small expansions
strique <- strique[strique$count > cutoff,]
df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "DN19_C2_histogram_STRique_hg38.png")

# no cutoff
strique <- read.table("c9orf72_DN19_C2.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "DN19_C2_histogram_STRique_hg38_noCutoff.png")


#=======================#


###
# STCl49
###
strique <- read.table("c9orf72_STCl49.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
# remove small expansions
strique <- strique[strique$count > cutoff,]
df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "STCl49_histogram_STRique_hg38.png")

# no cutoff
strique <- read.table("c9orf72_STCl49.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "STCl49_histogram_STRique_hg38_noCutoff.png")


#=======================#
