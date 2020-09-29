##########
# Glutamatergic Synapse
# Enrichment plot of TF targets
##########
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
library(reshape2)
library(stringr)
library(gridExtra)
setwd("~/rimod/paper_v2/figures/supplements/glutamatergic_synapse/")


# load data
tgts <- read.csv("EGR3_targets_enrichment.csv")

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
    labs(title = title, x="", y="-log10 adj. P-value") +
    theme(axis.text = element_text(size=16))
  return(p)
}



no_terms = 10

p1 <- make_enrichment_plot(tgts, source="GO:CC", cut=no_terms, title="EGR3 Targets")
ggsave("EGR3_targets_enrichment.png", width=7, height = 5)
