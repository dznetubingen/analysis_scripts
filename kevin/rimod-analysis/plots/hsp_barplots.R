###############################
# Figure 4: Barplots for Heat shock proteins
################################
library(stringr)
library(igraph)
library(ggplot2)

setwd("~/rimod/paper/figures/figure4/")


mapt <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/MAPT_cell_composition_filterd_DEGs.txt", sep="\t", header=T, stringsAsFactors = F)
grn <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/GRN_cell_composition_filterd_DEGs.txt", sep="\t", header=T, stringsAsFactors = F)

mapt <- mapt[grepl("HSP", mapt$hgnc_symbol),]
grn <- grn[grepl("HSP", grn$hgnc_symbol),]
mapt <- mapt[order(mapt$log2FoldChange),]

grn <- grn[order(grn$log2FoldChange),]
# PLOTTING
# 

df.mapt <- data.frame(HSP = factor(mapt$hgnc_symbol, levels=as.character(mapt$hgnc_symbol)), LFC = mapt$log2FoldChange, Pvalue = mapt$padj)

p <- ggplot(df, aes(x=HSP, y=LFC, fill = Pvalue)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_minimal()
p
ggsave("MAPT_HSP_barplot.png", width=3, height=3)



#GRN
df.grn <- data.frame(HSP = factor(grn$hgnc_symbol, levels=as.character(grn$hgnc_symbol)), LFC = grn$log2FoldChange, Pvalue = grn$padj)

p <- ggplot(df, aes(x=HSP, y=LFC, fill = Pvalue)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_minimal()
p
ggsave("GRN_HSP_barplot.png", width=3, height=3)

df.mapt$Group <- rep("FTD-MAPT", nrow(df.mapt))
df.grn$Group <- rep("FTD-GRN", nrow(df.grn))

df <- rbind(df.mapt, df.grn)
p <- ggplot(df, aes(x=HSP, y=LFC, fill = Pvalue)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45)) +
  facet_grid(~Group)
p
ggsave("HSP_barplot.png", width=4, height=3)

###
# also make upset plot of promoter shifting genes
###
library(UpSetR)
mapt <- na.omit(read.table("~/rimod/CAGE/cage_analysis/promotor_shifting/frontal/mapt_promotor_shifting_genes_fro.txt"))
grn <- na.omit(read.table("~/rimod/CAGE/cage_analysis/promotor_shifting/frontal/grn_promotor_shifting_genes_fro.txt"))
c9 <- na.omit(read.table("~/rimod/CAGE/cage_analysis/promotor_shifting/frontal/c9_promotor_shifting_genes_fro.txt"))

promshift <- list("FTD-MAPT" = mapt$V1,
                  "FTD-GRN" = grn$V1,
                  "FTD-C9orf72" = c9$V1)

upset(fromList(promshift), order.by = "freq")

