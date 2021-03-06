#################
# Figure 3 Plots
#################
library(VennDiagram)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
setwd("~/rimod/paper/figures/figure3/")


#####
# smRNA-seq DEG barplot
#####
mapt <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_result_mapt.ndc_frontal_smRNAseq.txt",
                   sep="\t", header=T, row.names=1)
grn <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_result_grn.ndc_frontal_smRNAseq.txt",
                   sep="\t", header=T, row.names=1)
c9 <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_result_c9.ndc_frontal_smRNAseq.txt",
                   sep="\t", header=T, row.names=1)
pval <- 0.05
mapt <- mapt[mapt$padj <= pval,]
grn <- grn[grn$padj <= pval,]
c9 <- c9[c9$padj <= pval,]
mir.list <- list("FTD-MAPT"=rownames(mapt),
                 "FTD-GRN"=rownames(grn),
                 "FTD-C9orf72"=rownames(c9))

# color palette
mypal <- brewer.pal(3, "Dark2")

#  make VENN plot
venn.plot <- venn.diagram(mir.list, 
                          filename = "smRNAseq_overlap_venn.png",
                          cex = 1,
                          
                          # image
                          imagetype = "png",
                          resolution = 300,
                          height=800,
                          width=800,
                          
                          # circles
                          lty = 'blank',
                          fill = mypal,
                          
                          # names
                          cat.cex = 1,
                          cat.pos = c(-27, 27, 135),
                          cat.dist = c(0.055, 0.055, 0.055))

#====================================================================#

###
# miRNA Reactome enrichment results
###
mapt <- read.csv("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/reactome_pathway_enrichment/mapt_downTargets_reactome.csv")
grn <- read.csv("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/reactome_pathway_enrichment/grn_downTargets_reactome.csv")
c9 <- read.csv("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/reactome_pathway_enrichment/c9orf72_downTargets_reactome.csv")
mapt <- mapt[1:2,]
grn <- grn[1:2,]
c9 <- c9[1:2,]
# change ordering
mapt <- mapt[order(mapt$negative_log10_of_adjusted_p_value),]
grn <- grn[order(grn$negative_log10_of_adjusted_p_value),]
c9 <- c9[order(c9$negative_log10_of_adjusted_p_value),]

mapt <- mapt[, c(2,5)]
mapt$Group <- rep("FTD-MAPT", nrow(mapt))
grn <- grn[, c(2,5)]
grn$Group <- rep("FTD-GRN", nrow(grn))
c9 <- c9[, c(2,5)]
c9$Group <- rep("FTD-C9orf72", nrow(c9))
df <- rbind(mapt, grn, c9)
colnames(df)[2] <- "negLog"

terms <- as.character(df$term_name)
term_levels <- terms[!duplicated(terms)]
terms <- factor(terms, levels = term_levels)


df$term_name <- terms
df$Group <- factor(df$Group, levels = c("FTD-MAPT", "FTD-GRN", "FTD-C9orf72"))

p <- ggplot(df, aes(x=term_name, y=negLog, fill=Group)) + 
  geom_bar(stat = "identity") + 
  theme_minimal(base_size = 15) +
  facet_grid(. ~ Group) +
  coord_flip() +
  scale_fill_brewer(palette = "Dark2")
p


