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
gene.mat <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2020-05-04_15.45.57/deseq_vst_values_2020-05-04_15.45.57.txt", sep="\t", header=T, row.names = 1)
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
    if (length(mod.genes) > 4){
      res <- bootstrap.enrichment.test(ctd, hits=mod.genes, bg = genes, genelistSpecies = "human", sctSpecies = "human", reps=1000)
      res <- res$results
      
      mod.enrichment.list$tmp <- res
      mod.name <- paste0(enrichment.name, "_", m)
      names(mod.enrichment.list)[length(mod.enrichment.list)] <- mod.name
    }
    else {
      print("not enough genes in module")
    }
    
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
  mcor <- median(na.omit(cors))
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

modules <- grn.down.modules

for (m in modules$CLUSTER_NAME) {
  cat(paste("\t", m, " \t"))
  print("")
  mod.genes <- getModule(modules, m)
  mod.exp <- gene.mat[mod.genes,]
  
  calcCellTypeCor(mod.exp, fracs)
}

####
# make heatmaps
####
library(pheatmap)
library(viridis)
setwd("~/rimod/paper/figures/figure3/")
makeDataFrame <- function(enr, dir="up", prefix="GRN_M"){
  df <- data.frame(enr[[1]]$p)
  for (i in 2:length(enr)) {
    tmp <- data.frame(enr[[i]]$p)
    df <- cbind(df, tmp)
  }
  colnames(df) <- paste0(prefix, c(1:length(enr)), "-",dir)
  rownames(df) <- rownames(enr[[1]])
  return(df)
}

pvalue_cutoff = 0.1

# GRN
grn.up.df <- makeDataFrame(grn.up.enrichment, dir="up", prefix="GRN_M")
grn.down.df <- makeDataFrame(grn.down.enrichment, dir="down", prefix="GRN_M")
grn.df <- cbind(grn.up.df, grn.down.df)
# remove Unknown
grn.df <- grn.df[-nrow(grn.df),]

grn.df[grn.df > pvalue_cutoff] <- NA

pheatmap(grn.df, color = viridis(200, option="D"), cluster_rows = F, cluster_cols = F, angle_col = "90",
         height = 3, width = 5, filename = "ewce_heatmap_grn.png")

# MAPT
mapt.up.df <- makeDataFrame(mapt.up.enrichment, dir="up", prefix="MAPT_M")
mapt.down.df <- makeDataFrame(mapt.down.enrichment, dir="down", prefix="MAPT_M")
mapt.df <- cbind(mapt.up.df, mapt.down.df)
mapt.df <- mapt.df[-nrow(mapt.df),]

mapt.df[mapt.df > pvalue_cutoff] <- NA

pheatmap(mapt.df, color = viridis(200, option="D"), cluster_rows = F, cluster_cols = F, angle_col = "90",
         height = 3, width = 5, filename = "ewce_heatmap_mapt.png")
