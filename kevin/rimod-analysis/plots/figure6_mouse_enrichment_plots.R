#######
# Plotting of pathway analysis results from mouse analysis for
# Figure 6
####
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
library(reshape2)
library(stringr)
library(gridExtra)

setwd("~/rimod/paper/figures/figure6/")

# load GRN old mouse data
grn <- read.csv("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/mouse/grn_mice/deseq_analysis/gProfiler_GRN_old_DEGs.csv")


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
    theme(text = element_text(size=15)) +
    coord_flip() +
    scale_color_gradient(low="red", high="blue") +
    labs(title = title, x="", y="-log10 adj. P-value")
  return(p)
}



no_terms = 10

p <- make_enrichment_plot(grn, source="GO:BP", cut = 20, title = "Old GRN Knockout Mice")

p
ggsave("GRN_Mouse_old_enrichment_GO:BP.png", width=8, height=5)


####
# Mapt middle aged mice
####
mapt <- read.csv("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/mouse/mapt_mice/deseq_results/gProfiler_MAPT_middle_mice.csv")

no_terms = 10

p <- make_enrichment_plot(mapt, source="GO:BP", cut = 20, title = "Middle-aged MAPT-transgenic Mice")
p
ggsave("MAPT_mouse_middle_enrichment_GO:BP.png", width=8, height=5)
