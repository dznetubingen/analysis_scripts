################################################
## Celltype Specific LFC correction          ###
################################################
library(biomaRt)
setwd("~/rimod/CAGE/cage_analysis/celltype_enrichment/ct_foldchange_correction/")

#####
# Preprocessing
######
# Load specificity matrix
cells <- read.csv("~/rimod/CAGE/cage_analysis/celltype_enrichment/specificity_values_for_KI_scRNAseq_superset_suppl_schizopreniapaper_level1.csv")
colnames(cells)[1] <- "mgi_symbol"

# Get human orthologs
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
lds <- getLDS(attributes=c("ensembl_gene_id", "mgi_symbol"), filters = "mgi_symbol", values = as.character(cells$mgi_symbol), mart = mouse, attributesL = c("ensembl_gene_id"), martL = human)
colnames(lds) <- c("mmu_gene_id", "mgi_symbol", "hsa_gene_id")
cells <- merge(cells, lds, by.x="mgi_symbol", by.y="mgi_symbol")

# Remove MGI symbol and put HSA symbols as rownames
cells <- cells[!duplicated(cells$hsa_gene_id),]
rownames(cells) <- cells$hsa_gene_id
cells <- cells[,c(-1, -26, -27)]

# Load count matrix
#exp <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2018-04-26_genewise/deseq_rLog_values_2018-04-26_14.22.04.txt", sep="\t", header = T, row.names=1)
# Load DEG table
deg <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2018-04-26_genewise/deseq_result_mapt.ndc_2018-04-26_14.22.04.txt")
deg <- deg[rownames(deg) %in% rownames(cells),]
outname = "mapt"
#####
# Correlation analysis
#####
# Parameters
no_genes <- 50

ct_fcs <- list()
# Calculate correlations
for (i in 1:ncol(cells)) {
  ct <- colnames(cells)[i]
  print(ct)
  tmp <- cells[order(cells[,i], decreasing = T),]
  tmp <- tmp[1:no_genes,]
  tmp.genes <- rownames(tmp)
  tmp.deg <- na.omit(deg[tmp.genes,])
  ct_fcs <- append(ct_fcs, mean(tmp.deg$log2FoldChange))
  names(ct_fcs)[i] <- ct
}
print(ct_fcs)

ct_fcs <- as.numeric(unlist(ct_fcs))

adj_fcs <- c()

for (i in 1:nrow(deg)) {
  gene <- rownames(deg)[i]
  lfc <- deg[i,]$log2FoldChange
  spec <- cells[gene,]
  sub_lfc <- sum(spec * ct_fcs)
  alfc <- lfc - sub_lfc
  adj_fcs <- c(adj_fcs, alfc)
}

deg$adj_lfc <- adj_fcs

# Filter using p-value and fold change
deg.filt <- deg[deg$padj <= 0.05,]
deg.filt <- deg.filt[abs(deg.filt$adj_lfc) >=1 ,]
deg.filt.up <- deg.filt[deg.filt$adj_lfc > 0, ]
deg.filt.down <- deg.filt[deg.filt$adj_lfc < 0, ]


write.table(rownames(deg.filt.up), paste("sig_Upgenes_ct_corrected_FC_", outname,".txt", sep=""), sep="\t", quote = F, row.names = F, col.names = F)
write.table(rownames(deg.filt.down), paste("sig_Downgenes_ct_corrected_FC_", outname,".txt", sep=""), sep="\t", quote = F, row.names = F, col.names = F)

write.table(rownames(deg.filt), paste("sig_genes_ct_corrected_FC_", outname,".txt", sep=""), sep="\t", quote = F, row.names = F, col.names = F)
write.table(deg.filt, paste("deg_celltype_adjFC_", outname, ".txt", sep=""), sep="\t", quote=F, col.names = NA)

# Some 
# tsne <- Rtsne(t(exp), perplexity=15, check_duplicates = F)
# cols <- colnames(exp)
# cols[grepl("control", cols)] <- "black"
# cols[grepl("MAPT", cols)] <- "blue"
# cols[grepl("GRN", cols)] <- "red"
# cols[grepl("C9orf72", cols)] <- "green"
# plot(tsne$Y, col=cols)

