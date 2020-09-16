###########################
# Analysis of repeat expansion results
# using STRique and tandem-genotypes
###########################
library(stringr)
library(ggplot2)

setwd("/home/kevin/dzne/nanopore_results_300320/")

cutoff = 10

###
# Sample DN19Cl2
###
strique <- read.table("dn19/c9orf72_dn19cl2.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
# remove small expansions
strique <- strique[strique$count > cutoff,]
df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "DN19Cl2_multiplex_histogram_STRique_hg38.png")

# no cutoff
strique <- read.table("dn19/c9orf72_dn19cl2.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "DN19Cl2_multiplex_histogram_STRique_hg38_noCutoff.png")


#=======================#


###
# Sample STCL49
###
strique <- read.table("stcl/c9orf72_stcl49.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
# remove small expansions
strique <- strique[strique$count > cutoff,]
df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "STCL49_multiplex_histogram_STRique_hg38.png")

# no cutoff
strique <- read.table("stcl/c9orf72_stcl49.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "STCL49_multiplex_histogram_STRique_hg38_noCutoff.png")


#=======================#
