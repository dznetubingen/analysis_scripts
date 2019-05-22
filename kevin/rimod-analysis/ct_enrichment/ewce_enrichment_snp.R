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

setwd("~/rimod/CAGE/cage_analysis/celltype_enrichment/ewce_enrichment/")


# Load Zeisel, Chen, Campbell dataset
counts <- read.table("mnn_zcc_sce_normcounts.txt", sep="\t", header = T)
rownames(counts) <- counts$X
counts <- counts[,-1]
anno <- read.table("mnn_zcc_sce_celltypes.txt", sep="\t", header = T)
anno$cell_id <- colnames(counts)
anno$level1class <- anno$celltypes
anno$level2class <- anno$celltypes
rownames(anno) <- colnames(counts)

mnn_zcc <- list(exp = counts, annot = anno)

gene="Olig2"
snap_exp <- as.numeric(mnn_zcc$exp["Olig2",])
snap_cells <- as.character(mnn_zcc$annot$level1class)
cellExpDist = data.frame(e=snap_exp, l1=snap_cells)
ggplot(cellExpDist) + geom_boxplot(aes(x=l1,y=e)) + xlab("Cell type") + ylab("Unique Molecule Count") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Generate celltype data for just the cortex/hippocampus data
# Drop genes which do not show significant evidence of varying between celltypes
exp_dropped = drop.uninformative.genes(exp=mnn_zcc$exp, level2annot = mnn_zcc$annot$level1class)
annotLevels = list(level1class=mnn_zcc$annot$level1class,level2class=mnn_zcc$annot$level2class)
# Calculate Celltype specificity matrices
ctd_fname = generate.celltype.data(exp=exp_dropped,annotLevels=annotLevels,groupName="MNN_ZCC")
print(ctd_fname)
ctd_fname = filter.genes.without.1to1.homolog(ctd_fname)
print(ctd_fname)
load(ctd_fname[1])

ctd_name_allki <- "~/rimod/CAGE/cage_analysis/celltype_enrichment/snp_enrichment/ctdFiles/ctd_allKI.rda"
load(ctd_name_allki)

# Some plotting
set.seed(1234)
library(reshape2)
genes = c("Snap25","Mag", "Ascl1")
exp = melt(cbind(ctd[[1]]$mean_exp[genes,],genes),id.vars="genes")
colnames(exp) = c("Gene","Cell","AvgExp")
ggplot(exp)+geom_bar(aes(x=Cell,y=AvgExp),stat="identity")+facet_grid(Gene~.)+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
exp = melt(cbind(data.frame(ctd[[1]]$specificity[genes,]),genes),id.vars="genes")
colnames(exp) = c("Gene","Cell","Expression")
ggplot(exp)+geom_bar(aes(x=Cell,y=Expression),stat="identity")+facet_grid(Gene~.)+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Load geneset
deg_file <- "~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2018-04-26_genewise/deseq_result_mapt.ndc_2018-04-26_14.22.04.txt"
outname <- "mapt"
deg <- read.table(deg_file, sep="\t", header = T, row.names = 1)
deg <- deg[deg$padj <= 0.05,]
deg <- deg[abs(deg$log2FoldChange) >= 0.5,]
genes_up <- rownames(deg[deg$log2FoldChange > 0,])
genes_down <- rownames(deg[deg$log2FoldChange < 0,])
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
bm.up <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = genes_up, mart = ensembl)
genes_up_symbol <- as.character(bm.up$hgnc_symbol)
bm.down <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = genes_down, mart = ensembl)
genes_down_symbol <- as.character(bm.down$hgnc_symbol)
# Generate mouse gene list
data("mouse_to_human_homologs")
m2h = unique(mouse_to_human_homologs[, c("HGNC.symbol", "MGI.symbol")])
# Parameters
reps = 10000
level = 1

# Run EWCE for up-regulated genes
mouse.hits <- unique(m2h[m2h$HGNC.symbol %in% genes_up_symbol, "MGI.symbol"])
mouse.bg  = unique(m2h$MGI.symbol)
full_results <- bootstrap.enrichment.test(sct_data = ctd, hits = mouse.hits, bg = mouse.bg, reps = reps, annotLevel = level)
res_up <- full_results$results[order(full_results$results$p),]
print(res_up)
write.table(res_up, paste("EWCE_result_upreg_",outname,"_genes_allKIctd.txt", sep=""), sep="\t", quote=F)
# Plot results
png(paste("EWCE_result_upreg_", outname, "_allKictd.png", sep=""), height=1000, width=2000)
ewce.plot(total_res = full_results$results)
dev.off()


# Run EWCE for down-regulated genes
mouse.hits <- unique(m2h[m2h$HGNC.symbol %in% genes_down_symbol, "MGI.symbol"])
mouse.bg  = unique(m2h$MGI.symbol)
full_results <- bootstrap.enrichment.test(sct_data = ctd, hits = mouse.hits, bg = mouse.bg, reps = reps, annotLevel = level)
res_down <- full_results$results[order(full_results$results$p),]
print(res_down)
write.table(res_down, paste("EWCE_result_downreg_",outname,"_genes_allKIctd.txt", sep=""), sep="\t", quote=F)
# Plot results
png(paste("EWCE_result_downreg_", outname, "_allKictd.png", sep=""), height=1000, width=2000)
ewce.plot(total_res = full_results$results)
dev.off()



# Do similar thing with SNPs
# Load snps
snps <- read.table("~/rimod/integrative_analysis/snp_integration/neurodegen_dementia.txt", sep="\t", header = T, fill = T)
genes <- as.character(snps$Gene.Locus)
genes <- genes[!genes == ""]
ul_genes <- c()
for (g in genes) {
  sp <- strsplit(g, split=", ")[[1]]
  sp <- as.character(unlist(sp))
  ul_genes <- c(ul_genes, sp)
}
mouse.hits <- unique(m2h[m2h$HGNC.symbol %in% ul_genes, "MGI.symbol"])
mouse.bg  = unique(m2h$MGI.symbol)
reps = 10000
level = 1
full_results <- bootstrap.enrichment.test(sct_data = ctd, hits = mouse.hits, bg = mouse.bg, reps = reps, annotLevel = level)
res_snp <- full_results$results[order(full_results$results$p),]
print(res_snp)
write.table(res_snp, "EWCE_result_snp_genes_allKIctd.txt", sep="\t", quote=F)



