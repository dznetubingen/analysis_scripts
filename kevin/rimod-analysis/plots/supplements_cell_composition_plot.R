library(ggplot2)
library(reshape2)

setwd("~/rimod/paper/figures/supplements//")

# load data
fracs <- read.table("~/rimod/RNAseq/analysis/deconvolution/cdn_predictions.txt", sep="\t", header=T, row.names=1, check.names = F)
fracs <- fracs[,-1]
rownames(fracs) <- gsub("X","", rownames(fracs))

md <- read.table("~/rimod/RNAseq/frontal_samples.txt", sep="\t", header=T)
fracs <- fracs[rownames(fracs) %in% md$ids,]
md <- md[match(rownames(fracs), md$ids),]

fracs$Group <- md$mutated_gene
fracs$Sample <- rownames(fracs)

df <- melt(fracs, id.vars = c("Group", "Sample"))

colnames(df) <- c("Group", "Sample", "Celltype", "Fraction")

df$Group <- factor(as.character(df$Group), levels=c("control", "C9orf72", "MAPT", "GRN"))

ggplot(df, aes(fill=Celltype, y=Fraction, x=Sample)) + 
  geom_bar(stat="identity", position="fill") + 
  theme_minimal() + 
  theme(axis.text.x=element_blank()) +
  facet_grid(~ Group, scales="free", space="free_x")

ggsave("cell_composition_stacked_barplot.png", width=10, height=10)
