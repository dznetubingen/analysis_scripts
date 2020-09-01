############
# M1-down enrichment plots
############
library(ggplot2)
setwd("~/rimod/paper_v2/figures/figure6/")

mapt <- read.csv("mapt_m1down_gProfiler.csv")
grn <- read.csv("grn_m1down_gProfiler.csv")



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



no_terms = 20

p1 <- make_enrichment_plot(grn, cut=no_terms, title="FTD-GRN M1-down")
ggsave("grn_m1_down_GO:BP.png", width=7, height = 5)

p2 <- make_enrichment_plot(mapt, cut=no_terms, title="FTD-MAPT M1-down")
ggsave("mapt_m1_down_GO:BP.png", width=7, height=5)

p3 <- make_enrichment_plot(mapt.up, source="REAC", cut=no_terms, title="FTD-MAPT Up")
p4 <- make_enrichment_plot(mapt.down, source="REAC", cut=no_terms, title="FTD-MAPT Down")

grid.arrange(p1, p2, p3, p4, nrow=4)

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
mapt.up <- read.csv("~/rimod/RNAseq/analysis/pathways_analysis/gprofiler/mapt_ccf_degs_up.csv", stringsAsFactors = F)
mapt.down <- read.csv("~/rimod/RNAseq/analysis/pathways_analysis/gprofiler/mapt_ccf_degs_down.csv", stringsAsFactors = F)
grn.up <- read.csv("~/rimod/RNAseq/analysis/pathways_analysis/gprofiler/grn_ccf_degs_up.csv", stringsAsFactors = F)
grn.down <- read.csv("~/rimod/RNAseq/analysis/pathways_analysis/gprofiler/grn_ccf_degs_down.csv", stringsAsFactors = F)
c9 <- read.csv("~/rimod/RNAseq/analysis/pathways_analysis/gprofiler/c9orf72_all_degs.csv", stringsAsFactors = F)

mapt.up$Group = "MAPT Up"
mapt.down$Group = "MAPT Down"
grn.up$Group = "GRN Up"
grn.down$Group = "GRN Down"

mapt.up <- filter_terms(mapt.up, source="REAC")
mapt.down <- filter_terms(mapt.down, source="REAC")
grn.up <- filter_terms(grn.up, source="REAC")
grn.down <- filter_terms(grn.down, source="REAC")

df <- rbind(grn.down, mapt.down, grn.up, mapt.up)

p <- ggplot(data=df, aes(x=term_name, y=negative_log10_of_adjusted_p_value, color=adjusted_p_value)) + 
  geom_point(aes(size=term_size)) +
  theme_minimal() +
  coord_flip() +
  labs(title = "", x="", y="-log10 adj. P-value") +
  scale_color_gradient(low="red", high="blue")  +
  facet_wrap(~ Group, scales = "free")
p

ggsave("enrichment_plot.png", width=11, height=5)


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

