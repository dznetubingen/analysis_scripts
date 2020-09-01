##########
# Figure 2 pathway enrichment plots
##########
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
library(reshape2)
library(stringr)
library(gridExtra)
setwd("~/rimod/Methylation/frontal_methylation_0818/pathway_go_analysis/")


# load data
mapt.up <- read.csv("MAPT_up_CpGs_gProfiler.csv")
mapt.down <- read.csv("MAPT_down_CpGs_gProfiler.csv")
grn.up <- read.csv("GRN_up_CpGs_gProfiler.csv")
grn.down <- read.csv("GRN_down_CpGs_gProfiler.csv")
c9.up <- read.csv("C9orf72_up_CpGs_gProfiler.csv")
c9.down <- read.csv("C9orf72_down_CpGs_gProfiler.csv")

make_enrichment_plot <- function(res, cut=20, source="GO:BP", title=""){
  res <- res[, 1:8]
  res <- res[res$source == source,]
  if (nrow(res) >= cut){
    res <- res[1:cut,]  
  }
  res <- res[order(res$adjusted_p_value, decreasing = T),]
  res$term_name <- factor(as.character(res$term_name), levels=as.character(res$term_name))
  
  p <- ggplot(data=res, aes(x=term_name, y=negative_log10_of_adjusted_p_value, color=adjusted_p_value)) + 
    geom_point(aes(size=term_size)) +
    theme_minimal() +
    coord_flip() +
    scale_color_gradient(low="red", high="blue") +
    labs(title = title, x="", y="-log10 adj. P-value")
  return(p)
}



no_terms = 10
# GRN
p1 <- make_enrichment_plot(grn.up,  cut=no_terms, title="FTD-GRN Up")
ggsave("grn_up_enrichment.png", width=6, height = 3)
p2 <- make_enrichment_plot(grn.down, cut=no_terms, title="FTD-GRN Down")
# MAPT
ggsave("grn_down_enrichment.png", width=6, height=3)
p3 <- make_enrichment_plot(mapt.up,  cut=no_terms, title="FTD-MAPT Up")
ggsave("mapt_up_enrichment.png", width=6, height=3)
p4 <- make_enrichment_plot(mapt.down, cut=no_terms, title="FTD-MAPT Down")
ggsave("mapt_down_enrichment.png", width=6, height=3)
# C9orf72
p5 <- make_enrichment_plot(c9.up,  cut=no_terms, title="FTD-C9orf72 Up")
ggsave("c9orf72_up_enrichment.png", width=6, height=3)
p6 <- make_enrichment_plot(c9.down, cut=no_terms, title="FTD-C9orf72 Down")
ggsave("c9orf72_down_enrichment.png", width=6, height=3)


grid.arrange(p1, p2, p3, p4, p5, p6, nrow=4)

# Try to merge all together
filter_terms <- function(res, cut=10, source = "GO:BP", length_cutoff = 50){
  res <- res[res$source == source,]
  if (nrow(res) >= cut){
    res <- res[1:cut,]  
  }
  
  for (i in 1:nrow(res)) {
    if (nchar(res$term_name[i]) > length_cutoff){
      res$term_name[i] <- res$term_id[i]
    }
  }
  
  res <- res[order(res$adjusted_p_value, decreasing = T),]
  res$term_name <- factor(as.character(res$term_name), levels=as.character(res$term_name))
  
  return(res)
}


# load data
mapt.up <- read.csv("MAPT_up_CpGs_gProfiler.csv", stringsAsFactors = F)
mapt.down <- read.csv("MAPT_down_CpGs_gProfiler.csv", stringsAsFactors = F)
grn.up <- read.csv("GRN_up_CpGs_gProfiler.csv", stringsAsFactors = F)
grn.down <- read.csv("GRN_down_CpGs_gProfiler.csv", stringsAsFactors = F)
c9.up <- read.csv("C9orf72_up_CpGs_gProfiler.csv", stringsAsFactors = F)
c9.down <- read.csv("C9orf72_down_CpGs_gProfiler.csv", stringsAsFactors = F)

mapt.up$Group = "MAPT Up"
mapt.down$Group = "MAPT Down"
grn.up$Group = "GRN Up"
grn.down$Group = "GRN Down"
c9.up$Group = "C9orf72 Up"
c9.down$Group = "C9orf72 Down"

mapt.up <- filter_terms(mapt.up)
mapt.down <- filter_terms(mapt.down)
grn.up <- filter_terms(grn.up)
grn.down <- filter_terms(grn.down)
c9.up <- filter_terms(c9.up)
c9.down <- filter_terms(c9.down)

df <- rbind(grn.down, mapt.down, c9.down, grn.up, mapt.up, c9.up)

p <- ggplot(data=df, aes(x=term_name, y=negative_log10_of_adjusted_p_value, color=adjusted_p_value)) + 
  geom_point(aes(size=term_size)) +
  theme_minimal() +
  coord_flip() +
  labs(title = "", x="", y="-log10 adj. P-value") +
  scale_color_gradient(low="red", high="blue")  +
  facet_wrap(~ Group, scales = "free", ncol=2)
p

ggsave("enrichment_plot.png", width=11, height=8)


####
# Full plots for supplements
####
# load data
mapt.up <- read.csv("~/rimod/RNAseq/analysis/pathways_analysis/gprofiler/mapt_ccf_degs_up.csv", stringsAsFactors = F)
mapt.down <- read.csv("~/rimod/RNAseq/analysis/pathways_analysis/gprofiler/mapt_ccf_degs_down.csv", stringsAsFactors = F)
grn.up <- read.csv("~/rimod/RNAseq/analysis/pathways_analysis/gprofiler/grn_ccf_degs_up.csv", stringsAsFactors = F)
grn.down <- read.csv("~/rimod/RNAseq/analysis/pathways_analysis/gprofiler/grn_ccf_degs_down.csv", stringsAsFactors = F)
c9 <- read.csv("~/rimod/RNAseq/analysis/pathways_analysis/gprofiler/c9orf72_all_degs.csv", stringsAsFactors = F)

mapt.up$Group = "MAPT Up"
mapt.down$Group = "MAPT Down"
grn.up$Group = "GRN Up"
grn.down$Group = "GRN Down"

mapt.up <- filter_terms(mapt.up, source="GO:BP", cut=50)
mapt.down <- filter_terms(mapt.down, source="GO:BP", cut=50)
grn.up <- filter_terms(grn.up, source="GO:BP", cut=50)
grn.down <- filter_terms(grn.down, source="GO:BP", cut=50)

df <- rbind(grn.down, mapt.down, grn.up, mapt.up)

p <- ggplot(data=df, aes(x=term_name, y=negative_log10_of_adjusted_p_value, color=adjusted_p_value)) + 
  geom_point(aes(size=term_size)) +
  theme_minimal() +
  coord_flip() +
  labs(title = "", x="", y="-log10 adj. P-value") +
  scale_color_gradient(low="red", high="blue")  +
  facet_wrap(~ Group, scales = "free")
p

ggsave("enrichment_plot_GOBP_cut50.png", width=15, height=15)

