#======================================#
# Celltype Enrichment FTD =============#
#======================================#
library(pheatmap)
library(biomaRt)
library(viridis)
setwd("~/rimod/CAGE/cage_analysis/celltype_enrichment/")

# Load specificity matrix
cells <- read.csv("specificity_values_for_KI_scRNAseq_superset_suppl_schizopreniapaper_level1.csv")
colnames(cells)[1] <- "mgi_symbol"
# Get human orthologs
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
lds <- getLDS(attributes=c("ensembl_gene_id", "mgi_symbol"), filters = "mgi_symbol", values = as.character(cells$mgi_symbol), mart = mouse, attributesL = c("ensembl_gene_id"), martL = human)
colnames(lds) <- c("mmu_gene_id", "mgi_symbol", "hsa_gene_id")
cells <- merge(cells, lds, by.x="mgi_symbol", by.y="mgi_symbol")
# Remove duplicates
cells <- cells[!duplicated(cells$hsa_gene_id),]
rownames(cells) <- cells$hsa_gene_id
cells <- cells[,c(-1, -26, -27)] # remove mouse gene IDs 


# Parameters
deg_file <- "deseq_result_c9.ndc_2018-04-26_14.22.04.txt"
outname <- "c9orf72"
no_of_celltype_genes <- 500 # No. of celltype specific genes to consider
window_size <- 100 # size of the sliding window

# Get celltype specfic list for each celltype
celltypes <- list()
for (i in 1:ncol(cells)) {
  ct <- colnames(cells)[i]
  tmp <- cells[order(cells[,i], decreasing = T),]
  ct_genes <- rownames(tmp)[1:no_of_celltype_genes]
  celltypes$tmp <- ct_genes
  names(celltypes)[i] <- ct
}

# Load DEG data
deg <- read.table(deg_file, sep="\t", header = T, row.names = 1)
deg <- deg[deg$padj <= 0.05,]
deg <- deg[abs(deg$log2FoldChange) >= 0.5,]
# Rank by fold-change
deg <- deg[order(deg$log2FoldChange, decreasing = TRUE),]
deg_genes <- rownames(deg)

# Sliding-window enrichment

window_size_half = window_size/2
no_windows <- nrow(deg) - window_size

df <- data.frame(decoy = rep(1, no_windows))

for (i in 1:length(celltypes)) {
  ct <- names(celltypes)[i]
  print(ct)
  ct_genes <- as.character(unlist(celltypes[i]))
  ct_deg_overlaps <- c()
  
  for (j in 1:no_windows) {
    window_genes <- deg_genes[j:(j+window_size)]
    ovlp <- length(intersect(window_genes, ct_genes))
    ct_deg_overlaps <- c(ct_deg_overlaps, ovlp)
  }
  
  df <- cbind(df, ct_deg_overlaps)
  colnames(df)[i+1] <- ct
  
}
df <- df[,-1]

# Plot as heatmap
lfc_change_col <- nrow(deg[deg$log2FoldChange > 0,])
pheatmap(t(df), cluster_cols = F, color = viridis(200, option="A"), gaps_col = lfc_change_col)

png(paste(outname, "_celltype_enrichment_heatmap.png", sep=""), width=1800, height=1000)
pheatmap(t(df), cluster_cols = F, color = viridis(200, option="A"), gaps_col = lfc_change_col)
dev.off()

