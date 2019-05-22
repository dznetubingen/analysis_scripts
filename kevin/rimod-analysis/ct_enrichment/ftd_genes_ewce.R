##################################
# EWCE analyis with FTD Genes ####
##################################
library(EWCE)
source("~/gitRepos/utils/utility_funs.R")

setwd("~/rimod/CAGE/cage_analysis/celltype_enrichment/ftd_genes_ewce/")


# Parameters
pval_cutoff <- 0.05
lfc_cutoff <- 0.5
reps <- 1000
level <- 1
outname <- "FTD_mendelian"
out_dir <- "~/rimod/CAGE/cage_analysis/celltype_enrichment/ftd_genes_ewce/"

setwd(out_dir)

# Load Karolinska Superset
ctd_name_allki <- "~/rimod/CAGE/cage_analysis/celltype_enrichment/snp_enrichment/ctdFiles/ctd_allKI.rda"
load(ctd_name_allki)

# Load geneset
gene_table <- read.table("FTD_mendelian.sorted.txt")
genes <- as.character(gene_table$V1)

# Generate mouse gene list for hits and background
data("mouse_to_human_homologs")
m2h = unique(mouse_to_human_homologs[, c("HGNC.symbol", "MGI.symbol")])

mouse.hits <- unique(m2h[m2h$HGNC.symbol %in% genes, "MGI.symbol"])
mouse.bg  = unique(m2h$MGI.symbol)

# Run EWCE for up-regulated genes
full_results <- bootstrap.enrichment.test(sct_data = ctd, hits = mouse.hits, bg = mouse.bg, reps = reps, annotLevel = level)
res <- full_results$results[order(full_results$results$p),]
print(res)
write.table(res, paste("EWCE_result_",outname,"_genes_allKIctd.txt", sep=""), sep="\t", quote=F)
# Plot results
png(paste("EWCE_result_upreg_", outname,"_",region,"_allKictd.png", sep=""), height=1000, width=2000)
utils.ewce.plot(total_res = full_results$results)
dev.off()