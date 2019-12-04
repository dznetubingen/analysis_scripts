##################
# Frontal RNA-seq alternative splicing analysis
##################

setwd("~/rimod/RNAseq/as_analysis/")


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

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

#===============================#
# RNA-SEQ cell composition DEGs #
#===============================#

###
# Data loading and formatting
###

# load fractions
fracs <- read.table("~/rimod/RNAseq/analysis/deconvolution/cdn_predictions.txt", sep="\t", header=T, row.names = 1)
fracs['Neurons'] <- fracs['ExNeurons'] + fracs['InNeurons']
fracs.sd <- apply(fracs, 2, sd)
fracs.mean <- apply(fracs, 2, mean)

# load expression
mat <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_vst_values_2019-10-23_13.33.11.txt", sep="\t", header=T ,row.names=1)
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
spec <- t(read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/darmanis_cell_type_specificity.txt", sep="\t", header=T, row.names=1, check.names = F))



####
# Check for genes strongly associated with cell composition differences
####
keep_vec <- c()

for (i in 1:nrow(mat)) {
  keep <- TRUE
  gene <- rownames(mat)[i]
  # check if specificity is available for gene
  if (gene %in% rownames(spec)){
    gene_spec <- spec[gene,]
    
    # check if specificity is higher than cutoff
    if (gene_spec > spec_cutoff){
      # get cell type with maximum specificity
      tmp <- spec[gene,]
      tmp.ct <- names(which(tmp == max(tmp)))
      # calc correlation
      exp <- as.numeric(mat[i,])
      f = as.numeric(unlist(fracs[tmp.ct]))
      gene_cor <- cor(exp, f)
      
      if (gene_cor > cor_cutoff){
        keep <- FALSE
      }
    }
  }
  keep_vec <- c(keep_vec, keep)
}




###
# Check for DEGs that are likely caused by cell composition differences
# @param cor_cutoff: cutoff for positive correlation to be called FP
# @param spec_cutoff: specificity cutoff to be considered a candidate
#
checkFalsePositives <- function(deg, cor_cutoff = 0.4, spec_cutoff = 0.4){
  
  # check for cell-type specific genes
  deg.spec <- deg[deg$hgnc_symbol %in% rownames(spec),] # keep only genes with specificity values available
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

# Remove cell-composition correlated genes
mat <- mat[keep_vec,]
genes <- rownames(mat)

# Load mapt
mapt <- read.table("majiq/mapt_AS_genes_dPSI_0.2.txt", sep="\t", header=T)
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=mapt$x, mart=ensembl)
mapt <- bm
mapt <- mapt[!mapt$hgnc_symbol == "",]
write.table(mapt$hgnc_symbol, "majiq/mapt_AS_genes_dPSI_0.2_CCF_HGNC.txt", quote=F, col.names = F, row.names = F)

# load GRN
grn <- read.table("majiq/grn_AS_genes_dPSI_0.2.txt", sep="\t", header=T)
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=grn$x, mart=ensembl)
grn <- bm
grn <- grn[!grn$hgnc_symbol == "",]
write.table(grn$hgnc_symbol, "majiq/grn_AS_genes_dPSI_0.2_CCF_HGNC.txt", quote=F, col.names = F, row.names = F)

# load C9orf72
c9 <- read.table("majiq/c9orf72_AS_genes_dPSI_0.2.txt", sep="\t", header=T)
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=c9$x, mart=ensembl)
c9 <- bm
c9 <- c9[!c9$hgnc_symbol == "",]
write.table(c9$hgnc_symbol, "majiq/c9orf72_AS_genes_dPSI_0.2.txt_CCF_HGNC", quote=F, col.names = F, row.names = F)
