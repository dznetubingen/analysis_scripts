##########
# Figure 5 pathway enrichment plots
##########
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
library(reshape2)
library(stringr)
library(gridExtra)
setwd("~/rimod/paper/figures/figure5/")


# load data
grn <- read.csv("gProfiler_GRN_M3down.csv")


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

p <- make_enrichment_plot(grn, source="GO:BP", cut = 20, title = "FTD-GRN M3down")

p
ggsave("GRN_M3down_reactome_enrichment.png", width=10, height=5)
