########
# TF target enrichment plots for Figure 4
#######
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
library(reshape2)
library(stringr)
library(gridExtra)


setwd("~/rimod/CAGE/cage_analysis/tf_enrichment_analysis_050420/target_enrichment_results/")


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

# NFKB2
nfkb2 <- read.csv("NFKB2_gProfiler_enrichment.csv")
make_enrichment_plot(nfkb2, cut=10, source = "KEGG")
ggsave("nfkb2_enrichment.png", width=7, height=3)

# RELA
rela <- read.csv("RELA_gProfiler_enrichment.csv")
make_enrichment_plot(rela, cut=10, source = "KEGG")
ggsave("rela_enrichment.png", width=7, height=3)


# SP1
sp1 <- read.csv("SP1_gProfiler_enrichment.csv")
make_enrichment_plot(sp1, cut=10, source="GO:BP")
ggsave("sp1_enrichment.png", width=7, height=3)

# KLF3
klf3 <- read.csv("KLF3_gProfiler_enrichment.csv")
make_enrichment_plot(klf3, cut=10)
ggsave("klf3_enrichment.png", width=7, height=3)
