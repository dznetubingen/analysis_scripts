###########
# Make TF target networks for Figure 4 FTD-GRN immune system seciton
###########
library(biomaRt)
library(igraph)

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

setwd("~/rimod/paper_v2/figures/figure4/")


grn <- read.table("~/rimod/CAGE/cage_analysis/tf_enrichment_analysis_050420/GRN_up_targets.txt")
deg <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/GRN_cell_composition_filterd_DEGs.txt", sep="\t", header=T)
deg <- deg[abs(deg$log2FoldChange) > 1,]
grn <- grn[grn$TargetGene %in% deg$Row.names,]

tfs <- c("KLF3", "NFKB2", "RELA", "SP1", "CEBPD")

for (tf in tfs) {
  # get targets
  tmp <- grn[grn$TF == tf,]
  targets <- as.character(tmp$TargetGene)
  bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters = "ensembl_gene_id", values=targets, mart=ensembl)
  targets <- bm$hgnc_symbol
  targets <- targets[!targets == ""]
  
  
  # make network
  edges <- c()
  for (i in 1:length(targets)) {
    e <- c(tf, targets[i])
    edges <- c(edges, e)
  }
  g <- graph(edges)
  write_graph(g, file = paste0(tf, "_target_network.gml"), format = "gml")
}

