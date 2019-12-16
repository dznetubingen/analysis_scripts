#####
# Perform EWCE enrichment analysis with human snRNA-seq data
#####
library(EWCE)
library(stringr)
library(biomaRt)

# Get genes from a module
getModule <- function(modules, mod){
  genes <- str_split(modules[modules$CLUSTER_NAME == mod,]$CLUSTER_GENES, pattern=",")[[1]]
  return(genes)
}
# mart
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")



# Generate Celltype Data from Darmanis dataset
mat <- read.table("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/revision1_september19/rosmap_deconvolution/training_data/processed_data/GSE67835_norm_counts_all.txt",
                  sep="\t", header=T, row.names = 1)
ct <- read.table("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/revision1_september19/rosmap_deconvolution/training_data/processed_data/GSE67835_celltypes.txt",
                 sep="\t", header=T, row.names = 1)
mat <- t(mat)
# keep only genes that are expressed commonly
#rs <- apply(mat, 1, sum)
#keep <- rs > 100
#mat <- mat[keep,]

annotLevel <- list(l1=ct$Celltype)
ct_data <- generate.celltype.data(exp=mat, annotLevels = annotLevel, groupName = "Darmanis")
load(ct_data)

# Define Background set
gene.mat <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_vst_values_2019-10-23_13.33.11.txt", sep="\t", header=T, row.names = 1)
genes <- row.names(gene.mat)
genes <- str_split(genes, pattern="[.]", simplify = T)[,1]
rownames(gene.mat) <- genes

# get HGNC
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=genes, mart=ensembl)
gene.mat <- merge(gene.mat, bm, by.x="row.names", by.y="ensembl_gene_id")
gene.mat <- gene.mat[!duplicated(gene.mat$hgnc_symbol),]
gene.mat <- gene.mat[!gene.mat$hgnc_symbol == "",]
gene.mat <- na.omit(gene.mat)
genes <- gene.mat$hgnc_symbol
rownames(gene.mat) <- gene.mat$hgnc_symbol
gene.mat <- gene.mat[, c(-1, -ncol(gene.mat))]


# Load modules
grn.up.modules <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_grn_filtered_up_modules.txt", header=T, stringsAsFactors = F)
grn.down.modules <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_grn_filtered_down_modules.txt", header=T, stringsAsFactors = F)

mapt.up.modules <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_mapt_filtered_up_modules.txt", header=T, stringsAsFactors = F)
mapt.down.modules <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_mapt_filtered_down_modules.txt", header=T, stringsAsFactors = F)

# Get enrichment results for all modules
perform_module_ewce <- function(modules, enrichment.name){
  mod.enrichment.list <- list()
  for (m in modules$CLUSTER_NAME) {
    mod.genes <- getModule(modules, m)
    mod.genes <- mod.genes[mod.genes %in% genes]
    res <- bootstrap.enrichment.test(ctd, hits=mod.genes, bg = genes, genelistSpecies = "human", sctSpecies = "human", reps=1000)
    res <- res$results
    
    mod.enrichment.list$tmp <- res
    mod.name <- paste0(enrichment.name, "_", m)
    names(mod.enrichment.list)[length(mod.enrichment.list)] <- mod.name
  }
  return(mod.enrichment.list)
}

# GRN
grn.up.enrichment <- perform_module_ewce(grn.up.modules, "GRN_UP")
grn.down.enrichment <- perform_module_ewce(grn.down.modules, "GRN_DOWN")

# MAPT
mapt.up.enrichment <- perform_module_ewce(mapt.up.modules, "MAPT_UP")
mapt.down.enrichment <- perform_module_ewce(mapt.down.modules, "MAPT_DOWN")


# Calculate average correlation of module genes with cell type fractions

# load deconvolution results
fracs <- read.table("~/rimod/RNAseq/analysis/deconvolution/cdn_predictions.txt", sep="\t", header=T, row.names = 1)
fracs <- fracs[colnames(gene.mat),]
fracs <- fracs[, -1]

# Function to calculate the average correlation of 
# a bunch of genes with cell type fractions
meanCor <- function(exp, f){
  cors <- c()
  for (i in 1:nrow(exp)) {
    cors <- c(cors, cor(as.numeric(exp[i,]), f))
  }
  mcor <- mean(na.omit(cors))
  return(mcor)
}


# Calculate correlations of a module to all possible celltypes as available in fractions
calcCellTypeCor <- function(mod.exp, fracs){
  celltypes <- colnames(fracs)
  for (ct in celltypes) {
    f <- as.numeric(unlist(fracs[ct]))
    mcor <- meanCor(mod.exp, f)
    print(paste(ct, mcor, sep=" : "))
  }
}

modules <- mapt.down.modules

for (m in modules$CLUSTER_NAME) {
  cat(paste("\t", m, " \t"))
  print("")
  mod.genes <- getModule(modules, m)
  mod.exp <- gene.mat[mod.genes,]
  
  calcCellTypeCor(mod.exp, fracs)
}




