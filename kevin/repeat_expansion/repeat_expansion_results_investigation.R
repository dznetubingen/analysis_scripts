###########################
# Analysis of repeat expansion results
# using STRique and tandem-genotypes
###########################
library(stringr)
library(ggplot2)

setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/c9orf72_enrichment/results_fourEnzymes_DN19/")


# Load tandem-genotypes
cutoff = 10
tg <- read.table("TG_results/dn19_tg_disease.txt")
tg <- tg[1,]

# Visualize tandemo-genotypes output
tg.fwd <- as.numeric(unlist(str_split(tg$V7, pattern=",")))
tg.back <- as.numeric(unlist(str_split(tg$V8, pattern=",")))
# remove small expansions
tg.fwd <- tg.fwd[tg.fwd > cutoff]
tg.back <- tg.back[tg.back > cutoff]
tg.df <- data.frame(expansion = c(tg.fwd, tg.back), strand = c(rep("forwad", length(tg.fwd)), rep("reverse", length(tg.back))))
p <- ggplot(tg.df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "DN19_histogram_TG.png")



# Strique results
strique <- read.table("STRique_results_DN19/dn19_hg19.strique.tsv", sep="\t", header=T)
strique <- strique[order(strique$count, decreasing = T),]
# remove small expansions
strique <- strique[strique$count > cutoff,]
st.df <- data.frame(expansion = strique$count, strand = strique$strand)
p <- ggplot(st.df, aes(x=expansion, color=strand, fill=strand)) +
  geom_histogram(alpha=0.5, bins=50)
p
ggsave(filename = "DN19_histogram_STRique.png")






st.df <- st.df[!duplicated(st.df),]
dim(st.df)
