##########
# Highest LFC gene plots
#########
library(ggplot2)
library(reshape2)
setwd("~/rimod/paper/figures/supplements/")

mapt <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/MAPT_cell_composition_filterd_DEGs.txt", sep="\t", header=T)
grn <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/GRN_cell_composition_filterd_DEGs.txt", sep="\t", header=T)
c9 <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/C9orf72_cell_composition_filterd_DEGs.txt", sep="\t", header=T)

# add group column
mapt$group <- rep("FTD-MAPT", nrow(mapt))
grn$group <- rep("FTD-GRN", nrow(grn))
c9$group <- rep("FTD-C9orf72", nrow(c9))

mapt <- mapt[!mapt$hgnc_symbol == "",]
grn <- grn[!grn$hgnc_symbol == "",]
c9 <- c9[!c9$hgnc_symbol == "",]

# order according to LFC
mapt <- mapt[order(abs(mapt$log2FoldChange), decreasing = T),]
grn <- grn[order(abs(grn$log2FoldChange), decreasing = T),]
c9 <- c9[order(abs(c9$log2FoldChange), decreasing = T),]

# cut out only top genes
cutoff = 40
mapt <- mapt[1:cutoff,]
grn <- grn[1:cutoff,]
c9 <- c9[1:cutoff,]

# remove some columns
mapt <- mapt[, c(3, 8, 9)]
grn <- grn[, c(3, 8, 9)]
c9 <- c9[, c(3, 8, 9)]

mypal <- c("#7570B3", "#db6e1a","#67e08a")
df <- rbind(melt(mapt), melt(grn), melt(c9))
df <- df[]

ggplot(df, aes(x=hgnc_symbol, y=value, fill=group)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = mypal) +
  facet_wrap(~ group, scale="free") +
  coord_flip() +
  theme_minimal() +
  labs(x="Gene", y="Log Fold Change")


ggsave("top_LFC_genes.png", width=10, height=10)
