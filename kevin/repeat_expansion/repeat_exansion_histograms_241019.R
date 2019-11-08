##################
# Analysis of second run of repeat expansions
# using tandem-genotypes only
###############################
library(stringr)
library(ggplot2)
setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/c9orf72_enrichment/")

cutoff = 10
# DN19CTL1
dn <- read.table("results_DN19CL1_1017/dn19_tg_disease.txt", sep="\t")
dn <- dn[1,]
# Visualize tandemo-genotypes output
dn.fwd <- as.numeric(unlist(str_split(dn$V7, pattern=",")))
dn.back <- as.numeric(unlist(str_split(dn$V8, pattern=",")))
# remove small expansions
dn.fwd <- dn.fwd[dn.fwd > cutoff]
dn.back <- dn.back[dn.back > cutoff]
dn.df <- data.frame(expansion = c(dn.fwd, dn.back), strand = c(rep("forwad", length(dn.fwd)), rep("reverse", length(dn.back))))
p <- ggplot(dn.df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "DN19CTl1_histogram_TG.png")

# Strique results
strique <- read.table("results_DN19CL1_1017/dn19cl1_hg19.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
# remove small expansions
strique <- strique[strique$count > cutoff,]
st.df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(st.df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "DN19_histogram_STRique.png")




# STCL49
dn <- read.table("results_STCL49_1017/stcl49_tg_disease.txt", sep="\t")
dn <- dn[1,]
# Visualize tandemo-genotypes output
dn.fwd <- as.numeric(unlist(str_split(dn$V7, pattern=",")))
dn.back <- as.numeric(unlist(str_split(dn$V8, pattern=",")))
# remove small expansions
dn.fwd <- dn.fwd[dn.fwd > cutoff]
dn.back <- dn.back[dn.back > cutoff]
dn.df <- data.frame(expansion = c(dn.fwd, dn.back), strand = c(rep("forwad", length(dn.fwd)), rep("reverse", length(dn.back))))
p <- ggplot(dn.df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=100)
p

ggsave(filename = "STCL49_histogram_TG.png")
