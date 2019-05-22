##############################################################
# Extract genes specific to Vascular Leptomeningeal Cells    #
# that are up-regulated in FTD                               #
###############################################################
setwd("~/rimod/CAGE/cage_analysis/celltype_enrichment/vlmc_function/")

library(biomaRt)
library(viridis)
library(pheatmap)

# Load specificity matrix
cells <- read.csv("~/rimod/CAGE/cage_analysis/celltype_enrichment/specificity_values_for_KI_scRNAseq_superset_suppl_schizopreniapaper_level1.csv")
colnames(cells)[1] <- "mgi_symbol"
# Get human orthologs
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
lds <- getLDS(attributes=c("ensembl_gene_id", "mgi_symbol"), filters = "mgi_symbol", values = as.character(cells$mgi_symbol), mart = mouse, attributesL = c("ensembl_gene_id"), martL = human)
colnames(lds) <- c("mmu_gene_id", "mgi_symbol", "hsa_gene_id")
cells <- merge(cells, lds, by.x="mgi_symbol", by.y="mgi_symbol")

# Load DEG file
deg_file <- "~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2018-04-26_genewise/deseq_result_mapt.ndc_2018-04-26_14.22.04.txt"

deg <- read.table(deg_file, sep="\t", header = T, row.names = 1)
deg <- deg[deg$padj <= 0.05,]
deg <- deg[deg$log2FoldChange >= 0,]
#deg <- deg[abs(deg$log2FoldChange) >= 1,]

vlmc <- cells[order(cells$Vascular.Leptomeningeal.Cells, decreasing = T),]
spec_cutoff = 0.5
vlmc_cut <- vlmc[vlmc$Vascular.Leptomeningeal.Cells >= spec_cutoff,]
write.table(vlmc_cut$hsa_gene_id, "vlmc_genes_spec0.5.txt", col.names = F, row.names = F, quote =F)
vlmc_genes <- vlmc_cut$hsa_gene_id
table(vlmc_genes %in% rownames(deg))
deg_vlmc <- vlmc_genes[vlmc_genes %in% rownames(deg)]

write.table(deg_vlmc, "c9_vlmc_deg_genes.txt", row.names = F, col.names = F, quote = F)
