library(biomaRt)

######################################################################
## Comparison of frontal/temporal specific genes to DEGs in FTD
#################################################################

setwd("~/rimod/CAGE/region_analysis/")

deg <- read.table('deseq_result_frotem_rest_controls.txt', sep="\t", header = T, row.names = 1)
all_degs <- deg
deg <- deg[deg$padj <= 0.001,]
deg <- deg[abs(deg$log2FoldChange) >= 1,]

genes <- rownames(deg)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values = genes, mart = ensembl)

mdata <- merge(deg, bm, by.x="row.names", by.y="ensembl_gene_id")
mdata[mdata$hgnc_symbol == "ENSG00000106333",]
deg[rownames(deg) == "ENSG00000106333",]

# Load FTD DEGs
deg.mapt <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2018-04-26_genewise/deseq_result_mapt.ndc_2018-04-26_14.22.04.txt", sep="\t", header = T, row.names=1)
deg.grn <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2018-04-26_genewise/deseq_result_grn.ndc_2018-04-26_14.22.04.txt", sep="\t", header = T, row.names=1)
deg.c9 <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2018-04-26_genewise/deseq_result_c9.ndc_2018-04-26_14.22.04.txt", sep="\t", header = T, row.names=1)

### Mapt intersection
deg.mapt <- deg.mapt[deg.mapt$padj <= 0.05,]
deg.mapt <- deg.mapt[abs(deg.mapt$log2FoldChange) >= 1,]
mapt_inter <- intersect(rownames(deg.mapt), rownames(deg))
#bm.mapt <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values = mapt_inter, mart = ensembl)
fmat <- matrix(c(length(mapt_inter), 
                 nrow(deg.mapt) - length(mapt_inter), 
                 nrow(deg), 
                 nrow(all_degs) - nrow(deg)), nrow=2, ncol=2)
fmat
fs.mapt <- fisher.test(fmat, alternative = "greater")
fs.mapt


### Grn intersection
deg.grn <- deg.grn[deg.grn$padj <= 0.05,]
deg.grn <- deg.grn[abs(deg.grn$log2FoldChange) >= 1,]
grn_inter <- intersect(rownames(deg.grn), rownames(deg))
#bm.grn <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values = grn_inter, mart = ensembl)
fmat <- matrix(c(length(grn_inter), 
                 nrow(deg.grn) - length(grn_inter), 
                 nrow(deg), 
                 nrow(all_degs) - nrow(deg)), nrow=2, ncol=2)
fmat
fs.grn <- fisher.test(fmat, alternative = "greater")
fs.grn


### C9orf72 intersection
deg.c9 <- deg.c9[deg.c9$padj <= 0.05,]
deg.c9 <- deg.c9[abs(deg.c9$log2FoldChange) >= 1,]
c9_inter <- intersect(rownames(deg.c9), rownames(deg))
#bm.c9 <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values = c9_inter, mart = ensembl)
fmat <- matrix(c(length(c9_inter), 
                 nrow(deg.c9) - length(c9_inter), 
                 nrow(deg), 
                 nrow(all_degs) - nrow(deg)), nrow=2, ncol=2)
fs.c9 <- fisher.test(fmat, alternative = "greater")
fs.c9

# Common genes
cmn <- intersect(grn_inter, intersect(mapt_inter, c9_inter))
bm.cmb <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values = cmn, mart = ensembl)

rld <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2018-04-26_genewise/deseq_rLog_values_2018-04-26_14.22.04.txt", sep="\t", header = T, row.names=1)

# some testing
cts <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2018-04-26_genewise/deseq_normalized_counts_2018-04-26_14.22.04.txt", sep="\t", header = T, row.names = 1)
slc <- cts[rownames(cts) == "ENSG00000128245",]
group <- as.character(sapply(colnames(slc), function(x){strsplit(x, split="chrM_")[[1]][[2]]}))
group[!group == "control"] <- "case"
slc <- slc[,order(group)]
group <- group[order(group)]
plot(as.numeric(slc))              



#Further investigate frotem DEGs
deg.up <- deg[deg$log2FoldChange > 0,]
deg.up <- deg.up[order(deg.up$log2FoldChange, decreasing = T),]
degu.genes <- rownames(deg.up)
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = degu.genes, mart = ensembl)
deg.up.annot <- merge(deg.up, bm, by.x="row.names", by.y="ensembl_gene_id")
deg.up.annot <- deg.up.annot[order(deg.up.annot$log2FoldChange,decreasing = T),]


inter_gm <- intersect(grn_inter, mapt_inter)
deg.igm <- deg[rownames(deg) %in% inter_gm,]

# some saing
write.table(rownames(deg.up), "frotem_upregulated_genes.txt", sep="\t", quote = F, row.names=F, col.names=F)
