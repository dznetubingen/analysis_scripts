######################################################################
# Methylation Cell-type enrichment
# CT enrichment analysis of differentially methylated positions
######################################################################
library(EWCE)
library(ggplot2)
library(cowplot)
library(limma)
library(readxl)
library(MAGMA.Celltyping)
library(viridis)
library(biomaRt)
source("~/gitRepos/utils/utility_funs.R")

# Parameters
pval_cutoff <- 0.05
lfc_cutoff <- 0.5
reps <- 1000
level <- 1
region <- "frontal"
out_dir <- "~/rimod/Methylation/frontal_methylation_0818/celltype_enrichment/"
setwd(out_dir)
outname = "mapt_p301l"

# Load data
dmps <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_maptP301L.ndc_funnNorm.txt")

# Load Karolinska Superset
ctd_name_allki <- "~/rimod/CAGE/cage_analysis/celltype_enrichment/snp_enrichment/ctdFiles/ctd_allKI.rda"
load(ctd_name_allki)

# Extrac gene accessions
getGenes <- function(x){
  genes <- as.character(x$GencodeCompV12_Accession)
  genes <- genes[!genes == ""]
  genes <- as.character(sapply(genes, function(y){strsplit(y, split="[.]")[[1]][[1]]}))
  genes <- genes[!duplicated(genes)]
  return(genes)
}


# Get all CpGs
all_cpgs = getGenes(dmps)

# Get significant up- and down-regulated CpGs
dmps.sig <- dmps[dmps$adj.P.Val <= 0.05,]
up_dmp <- dmps.sig[dmps.sig$logFC > 0, ]
down_dmp <- dmps.sig[dmps.sig$logFC < 0, ]
up_dmp <- getGenes(up_dmp)
down_dmp <- getGenes(down_dmp)

# Use biomaRt to get gene symbols
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
bm_up <- getBM(attributes = c("ensembl_transcript_id", "hgnc_symbol"), filters="ensembl_transcript_id", values=up_dmp, mart=ensembl)
bm_down <- getBM(attributes = c("ensembl_transcript_id", "hgnc_symbol"), filters="ensembl_transcript_id", values=down_dmp, mart=ensembl)
bm_all <- getBM(attributes = c("ensembl_transcript_id", "hgnc_symbol"), filters="ensembl_transcript_id", values=all_cpgs, mart=ensembl)
genes_up_symbol <- as.character(bm_up$hgnc_symbol)
genes_down_symbol <- as.character(bm_down$hgnc_symbol)
genes_all_symbol <- as.character(bm_all$hgnc_symbol)

# Generate mouse gene list
data("mouse_to_human_homologs")
m2h = unique(mouse_to_human_homologs[, c("HGNC.symbol", "MGI.symbol")])

# Run EWCE for up-regulated genes
mouse.hits <- unique(m2h[m2h$HGNC.symbol %in% genes_up_symbol, "MGI.symbol"])
mouse.bg = unique(m2h[m2h$HGNC.symbol %in% genes_all_symbol, "MGI.symbol"])

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
mouse.bg  = unique(m2h$MGI.symbol)
full_results <- bootstrap.enrichment.test(sct_data = ctd, hits = mouse.hits, bg = mouse.bg, reps = reps, annotLevel = level)
res_down <- full_results$results[order(full_results$results$p),]
print(res_down)
write.table(res_down, paste("EWCE_result_downreg_",outname,"_",region,"_genes_allKIctd_lev2.txt", sep=""), sep="\t", quote=F)
# Plot results
png(paste("EWCE_result_downreg_", outname,"_",region,"_allKictd_lev2.png", sep=""), height=1000, width=2000)
utils.ewce.plot(total_res = full_results$results)
dev.off()


utils.ewce.plot(total_res = full_results$results)
