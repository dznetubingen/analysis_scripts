##############################
# Compare DEGs and expression with cell type specificity
# to filter out DEGs that are likely artifcats of cell 
# composition variability
##############################
library(biomaRt)
library(stringr)

setwd("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/")

# Define cutoffs
cor_cutoff = 0.4
spec_cutoff = 0.5


#=========#
# RNA-SEQ #
#=========#

###
# Data loading and formatting
###

# load fractions
fracs <- read.table("../cdn_predictions.txt", sep="\t", header=T, row.names = 1)
fracs['Neurons'] <- fracs['ExNeurons'] + fracs['InNeurons']
fracs.sd <- apply(fracs, 2, sd)
fracs.mean <- apply(fracs, 2, mean)

# load expression
mat <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-08-12_07.58.35/deseq_vst_values_2019-08-12_07.58.35.txt", sep="\t", header=T ,row.names=1)
rownames(mat) <- str_split(rownames(mat), pattern="[.]", simplify = T)[,1]
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=rownames(mat), mart=ensembl)
mat <- merge(mat, bm, by.x="row.names", by.y="ensembl_gene_id")
mat <- mat[!duplicated(mat$hgnc_symbol),]
mat <- na.omit(mat)
rownames(mat) <- mat$hgnc_symbol
mat <- mat[, c(-1, -ncol(mat))]

# subset fracs
fracs <- fracs[colnames(mat),]

# Load specificity values
spec <- t(read.table("darmanis_cell_type_specificity.txt", sep="\t", header=T, row.names=1, check.names = F))

#============================================================================#

####
# Check for false positives and report them
####

###
# Check for DEGs that are likely caused by cell composition differences
# @param cor_cutoff: cutoff for positive correlation to be called FP
# @param spec_cutoff: specificity cutoff to be considered a candidate
#
checkFalsePositives <- function(deg, cor_cutoff = 0.4, spec_cutoff = 0.4){
  
  # check for cell-type specific genes
  deg.spec <- deg[deg$hgnc_symbol %in% rownames(cts),] # keep only genes with specificity values available
  deg.spec <- spec[deg.spec$hgnc_symbol,]
  deg.cs.max <- apply(deg.spec, 1, max)
  # select candidate genes
  candidate.genes <- deg.spec[deg.cs.max > spec_cutoff,]
  
  # Test candidates for positive correlation
  false_positives <- c()
  for (cg in rownames(candidate.genes)) {
    # Get highest expressing cell type
    tmp <- deg.spec[cg,]
    tmp.ct <- names(which(tmp == max(tmp)))
    # check for correlation
    exp = as.numeric(mat[cg,])
    f = as.numeric(unlist(fracs[tmp.ct]))
    res = cor(exp, f)
    
    # check if gene is positively correlated with cell composition
    if (res > cor_cutoff){
      false_positives <- c(false_positives, cg)
    }
  }
  print(paste0("Found ", length(false_positives), " false positive DEGs"))
  
  return(false_positives)
}


ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

## MAPT
# load DEGs
deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-08-12_07.58.35/deseq_result_mapt.ndc_fro_2019-08-12_07.58.35.txt", sep="\t", header=T, row.names=1)
deg <- deg[deg$padj <= 0.05,]
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=rownames(deg), mart=ensembl)
deg <- merge(deg, bm, by.x="row.names", by.y="ensembl_gene_id")

mapt.fps <- checkFalsePositives(deg, cor_cutoff = cor_cutoff, spec_cutoff = spec_cutoff)
# filter DEGs
deg <- deg[!deg$hgnc_symbol %in% mapt.fps,]
deg.up <- deg[deg$log2FoldChange > 0,]
deg.down <- deg[deg$log2FoldChange < 0,]
write.table(deg, "MAPT_cell_composition_filterd_DEGs.txt", sep="\t", quote=F)
write.table(deg.up$Row.names, "MAPT_filterd_DEGs_up.txt", row.names = F, quote=F, col.names = F)
write.table(deg.down$Row.names, "MAPT_filterd_DEGs_down.txt", row.names = F, quote=F, col.names = F)
# for string
df <- data.frame(gene = deg$Row.names, lfc = deg$log2FoldChange)
write.table(df, "MAPT_ccFiltered_DEGs_STRING.txt", sep="\t", row.names = F, quote=F)

## GRN
deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-08-12_07.58.35/deseq_result_grn.ndc_fro_2019-08-12_07.58.35.txt", sep="\t", header=T, row.names=1)
deg <- deg[deg$padj <= 0.05,]
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=rownames(deg), mart=ensembl)
deg <- merge(deg, bm, by.x="row.names", by.y="ensembl_gene_id")

grn.fps <- checkFalsePositives(deg, cor_cutoff = cor_cutoff, spec_cutoff = spec_cutoff)
deg <- deg[!deg$hgnc_symbol %in% grn.fps,]
deg.up <- deg[deg$log2FoldChange > 0,]
deg.down <- deg[deg$log2FoldChange < 0,]
write.table(deg, "GRN_cell_composition_filterd_DEGs.txt", sep="\t", quote=F)
write.table(deg.up$Row.names, "GRN_filterd_DEGs_up.txt", row.names = F, quote=F, col.names = F)
write.table(deg.down$Row.names, "GRN_filterd_DEGs_down.txt", row.names = F, quote=F, col.names = F)
df <- data.frame(gene = deg$Row.names, lfc = deg$log2FoldChange)
write.table(df, "GRN_ccFiltered_DEGs_STRING.txt", sep="\t", row.names = F, quote=F)


# C9ORF72
deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-08-12_07.58.35/deseq_result_c9.ndc_fro_2019-08-12_07.58.35.txt", sep="\t", header=T, row.names=1)
deg <- deg[deg$padj <= 0.05,]
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=rownames(deg), mart=ensembl)
deg <- merge(deg, bm, by.x="row.names", by.y="ensembl_gene_id")

# Check for FPs
c9.fps <- checkFalsePositives(deg, cor_cutoff = cor_cutoff, spec_cutoff = spec_cutoff)

# Filter DEGs
deg <- deg[!deg$hgnc_symbol %in% c9.fps,]
deg.up <- deg[deg$log2FoldChange > 0,]
deg.down <- deg[deg$log2FoldChange < 0,]
write.table(deg, "C9orf72_cell_composition_filterd_DEGs.txt", sep="\t", quote=F)
write.table(deg.up$Row.names, "C9orf72_filterd_DEGs_up.txt", row.names = F, quote=F, col.names = F)
write.table(deg.down$Row.names, "C9orf72_filterd_DEGs_down.txt", row.names = F, quote=F, col.names = F)
df <- data.frame(gene = deg$Row.names, lfc = deg$log2FoldChange)
write.table(df, "C9orf72_ccFiltered_DEGs_STRING.txt", sep="\t", row.names = F, quote=F)

#=============================================================================================================#

#####
# Perform filtering on CAGE-seq DEGs
#####


setwd("~/rimod/CAGE/cage_analysis/deconvolution_frontal/cell_specificity_filtering/")


###
# Data loading and formatting
###

# load fractions
fracs <- read.table("../cdn_predictions.txt", sep="\t", header=T, row.names = 1)
fracs['Neurons'] <- fracs['ExNeurons'] + fracs['InNeurons']
fracs.sd <- apply(fracs, 2, sd)
fracs.mean <- apply(fracs, 2, mean)

# load expression
mat <- read.table("../rimod_frontal_CAGEseq_VST_aggregated_hgnc.txt", sep="\t", header=T ,row.names=1)

# subset fracs
fracs <- fracs[colnames(mat),]

# Load specificity values
spec <- t(read.table("darmanis_cell_type_specificity.txt", sep="\t", header=T, row.names=1, check.names = F))

#============================================#

# Perform filtering

## MAPT
# load DEGs
deg <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2019-08-15_14.49.08_frontal/deseq_result_mapt.ndc_2019-08-15_14.49.08.txt", sep="\t", header=T, row.names=1)
deg <- deg[deg$padj <= 0.05,]
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=rownames(deg), mart=ensembl)
deg <- merge(deg, bm, by.x="row.names", by.y="ensembl_gene_id")

mapt.fps <- checkFalsePositives(deg, cor_cutoff = cor_cutoff, spec_cutoff = spec_cutoff)

# filter DEGs
deg <- deg[!deg$hgnc_symbol %in% mapt.fps,]
deg.up <- deg[deg$log2FoldChange > 0,]
deg.down <- deg[deg$log2FoldChange < 0,]
write.table(deg, "MAPT_cell_composition_filterd_DEGs.txt", sep="\t", quote=F)
write.table(deg.up$Row.names, "MAPT_filterd_DEGs_up.txt", row.names = F, quote=F, col.names = F)
write.table(deg.down$Row.names, "MAPT_filterd_DEGs_down.txt", row.names = F, quote=F, col.names = F)


## GRN
deg <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2019-08-15_14.49.08_frontal/deseq_result_grn.ndc_2019-08-15_14.49.08.txt", sep="\t", header=T, row.names=1)
deg <- deg[deg$padj <= 0.05,]
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=rownames(deg), mart=ensembl)
deg <- merge(deg, bm, by.x="row.names", by.y="ensembl_gene_id")

grn.fps <- checkFalsePositives(deg, cor_cutoff = cor_cutoff, spec_cutoff = spec_cutoff)
deg <- deg[!deg$hgnc_symbol %in% grn.fps,]
deg.up <- deg[deg$log2FoldChange > 0,]
deg.down <- deg[deg$log2FoldChange < 0,]
write.table(deg, "GRN_cell_composition_filterd_DEGs.txt", sep="\t", quote=F)
write.table(deg.up$Row.names, "GRN_filterd_DEGs_up.txt", row.names = F, quote=F, col.names = F)
write.table(deg.down$Row.names, "GRN_filterd_DEGs_down.txt", row.names = F, quote=F, col.names = F)


# C9ORF72
deg <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2019-08-15_14.49.08_frontal/deseq_result_c9.ndc_2019-08-15_14.49.08.txt", sep="\t", header=T, row.names=1)
deg <- deg[deg$padj <= 0.05,]
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=rownames(deg), mart=ensembl)
deg <- merge(deg, bm, by.x="row.names", by.y="ensembl_gene_id")

# Check for FPs
c9.fps <- checkFalsePositives(deg, cor_cutoff = cor_cutoff, spec_cutoff = spec_cutoff)

# Filter DEGs
deg <- deg[!deg$hgnc_symbol %in% c9.fps,]
deg.up <- deg[deg$log2FoldChange > 0,]
deg.down <- deg[deg$log2FoldChange < 0,]
write.table(deg, "C9orf72_cell_composition_filterd_DEGs.txt", sep="\t", quote=F)
write.table(deg.up$Row.names, "C9orf72_filterd_DEGs_up.txt", row.names = F, quote=F, col.names = F)
write.table(deg.down$Row.names, "C9orf72_filterd_DEGs_down.txt", row.names = F, quote=F, col.names = F)


