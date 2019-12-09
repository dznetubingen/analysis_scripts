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
rs <- apply(mat, 1, sum)
keep <- rs > 100
mat <- mat[keep,]

annotLevel <- list(l1=ct$Celltype)
ct_data <- generate.celltype.data(exp=mat, annotLevels = annotLevel, groupName = "Darmanis")
load(ct_data)


# Define Background set
genes <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_vst_values_2019-10-23_13.33.11.txt", sep="\t", header=T)
genes <- as.character(genes$X)
genes <- str_split(genes, pattern="[.]", simplify = T)[,1]

# get HGNC
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=genes, mart=ensembl)
genes <- bm$hgnc_symbol
genes <- genes[!genes == ""]

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
