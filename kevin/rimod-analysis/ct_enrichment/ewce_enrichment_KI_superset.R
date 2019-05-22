######################################
## Celltype Enrichment using EWCE ####
######################################
library(EWCE)
library(ggplot2)
library(cowplot)
library(limma)
library(readxl)
library(biomaRt)
library(MAGMA.Celltyping)
library(viridis)
source("~/gitRepos/utils/utility_funs.R")



# Parameters
pval_cutoff <- 0.05
lfc_cutoff <- 0.5
reps <- 1000
level <- 2
region <- "frontal"
out_dir <- "~/rimod/CAGE/cage_analysis/celltype_enrichment/ewce_enrichment/frontal_ewce_enrichment/"

setwd(out_dir)

# Load Karolinska Superset
ctd_name_allki <- "~/rimod/CAGE/cage_analysis/celltype_enrichment/snp_enrichment/ctdFiles/ctd_allKI.rda"
load(ctd_name_allki)

# Load geneset
deg_file <- "~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2018-04-26_genewise/deseq_result_mapt.ndc_2018-04-26_14.22.04.txt"
outname <- "mapt"
deg <- read.table(deg_file, sep="\t", header = T, row.names = 1)
deg_ns <- deg
deg <- deg[deg$padj <= pval_cutoff,]
deg <- deg[abs(deg$log2FoldChange) >= lfc_cutoff,]
genes_up <- rownames(deg[deg$log2FoldChange > 0,])
genes_down <- rownames(deg[deg$log2FoldChange < 0,])
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
bm.up <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = genes_up, mart = ensembl)
genes_up_symbol <- as.character(bm.up$hgnc_symbol)
bm.down <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = genes_down, mart = ensembl)
genes_down_symbol <- as.character(bm.down$hgnc_symbol)
# Get all detected genes
all_genes <- rownames(deg_ns)
bm.all <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = all_genes, mart = ensembl)
all_genes_symbol <- as.character(bm.all$hgnc_symbol)

# Generate mouse gene list
data("mouse_to_human_homologs")
m2h = unique(mouse_to_human_homologs[, c("HGNC.symbol", "MGI.symbol")])

# Run EWCE for up-regulated genes
mouse.hits <- unique(m2h[m2h$HGNC.symbol %in% genes_up_symbol, "MGI.symbol"])
mouse.bg  = unique(m2h[m2h$HGNC.symbol %in% all_genes_symbol, "MGI.symbol"])

full_results <- bootstrap.enrichment.test(sct_data = ctd, hits = mouse.hits, bg = mouse.bg, reps = reps, annotLevel = level)
res_up <- full_results$results[order(full_results$results$p),]
print(res_up)
write.table(res_up, paste("EWCE_result_upreg_",outname,"_",region,"_genes_allKIctd_lev2.txt", sep=""), sep="\t", quote=F)
# Plot results
png(paste("EWCE_result_upreg_", outname,"_",region,"_allKictd_lev2.png", sep=""), height=1000, width=2000)
utils.ewce.plot(total_res = full_results$results)
dev.off()


# Run EWCE for down-regulated genes
mouse.hits <- unique(m2h[m2h$HGNC.symbol %in% genes_down_symbol, "MGI.symbol"])
mouse.bg  = unique(m2h[m2h$HGNC.symbol %in% all_genes_symbol, "MGI.symbol"])
full_results <- bootstrap.enrichment.test(sct_data = ctd, hits = mouse.hits, bg = mouse.bg, reps = reps, annotLevel = level)
res_down <- full_results$results[order(full_results$results$p),]
print(res_down)
write.table(res_down, paste("EWCE_result_downreg_",outname,"_",region,"_genes_allKIctd_lev2.txt", sep=""), sep="\t", quote=F)
# Plot results
png(paste("EWCE_result_downreg_", outname,"_",region,"_allKictd_lev2.png", sep=""), height=1000, width=2000)
utils.ewce.plot(total_res = full_results$results)
dev.off()


utils.ewce.plot(total_res = full_results$results)

