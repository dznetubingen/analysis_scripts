##########
# enrichment plots for iPSC-neurons miRNA experiments
##########
library(ggplot2)


setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/mirna_mimics/neurons_gprofiler_results/")

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

####
# miR-19b
####
cutoff = 10

# mimic up
df <- read.csv("mir19b_mimic_up.csv")
p <- make_enrichment_plot(df, cut=cutoff, title="miR-19b Mimic up")
ggsave("miR19b_mimic_up.png", width=7, height=4)

# mimic down
df <- read.csv("mir19b_mimic_down.csv")
p <- make_enrichment_plot(df, cut=cutoff, title="miR-19b Mimic down")
ggsave("miR19b_mimic_down.png", width=6, height=4)

# inhib up
df <- read.csv("mir19b_inhib_up.csv")
p <- make_enrichment_plot(df, cut=cutoff, title="miR-19b inhib up")
ggsave("miR19b_inhib_up.png", width=6, height=4)

