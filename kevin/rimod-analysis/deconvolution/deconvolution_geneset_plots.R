#######################
# Deconvolution analysis plotting
#######################
library(stringr)
library(biomaRt)
library(pheatmap)
library(viridis)
library(nnls)

setwd("~/rimod/RNAseq/analysis/deconvolution/")

# Load expresion values
#mat <- read.table("frontal_lengthScaledTPM_counts.txt", sep="\t", row.names=1, header=T)
#rownames(mat) <- str_split(rownames(mat), pattern="[.]", simplify = T)[,1]
#colnames(mat) <- gsub("X", "", colnames(mat))
# Get HGNC symbols
#ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
#bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=rownames(mat), mart=ensembl)
#mat <- merge(mat, bm, by.x="row.names", by.y="ensembl_gene_id")
#hgnc <- mat$hgnc_symbol
#mat <- mat[, c(-1, -ncol(mat))]
#mat <- aggregate(mat, by=list(hgnc),  FUN = sum)
#rownames(mat) <- mat$Group.1
#mat <- mat[-1,-1]
#write.table(mat, "frontal_lengthScaledTPM_aggrHGNC.txt", sep="\t", quote=F)

# Load expression values (aggregated and HGNC, see above)
mat <- read.table("frontal_lengthScaledTPM_aggrHGNC.txt", sep="\t", header=T)
colnames(mat) <- gsub("X", "", colnames(mat))


# Load deconvolution results
fracs <- read.table("cdn_predictions.txt", sep="\t", row.names = 1, header=T)
rownames(fracs) <- gsub("X", "", rownames(fracs))
fracs$Neurons <- fracs$InNeurons + fracs$ExNeurons
fracs <- fracs[, c(-1, -2, -7,-8)]
fracs <- as.matrix(fracs)

# Load MD 
md <- read.table("~/rimod/RNAseq/rnaseq_frontal_md.txt", sep="\t", header=T, stringsAsFactors = F)
md <- md[match(rownames(fracs), md$ids),]
md <- md[-19,]

# keep mapt and control
# MAPT
keep <- as.character(unlist(md['mutated_gene'])) %in% c('MAPT')
mapt.ids <- as.character(unlist(md[keep,]['ids']))
mapt.fracs <- fracs[mapt.ids,]
mapt.mat <- mat[, mapt.ids]

# GRN
keep <- as.character(unlist(md['mutated_gene'])) %in% c('GRN')
grn.ids <- as.character(unlist(md[keep,]['ids']))
grn.fracs <- fracs[grn.ids,]
grn.mat <- mat[, grn.ids]

# C9orf72
keep <- as.character(unlist(md['mutated_gene'])) %in% c('C9orf72')
c9.ids <- as.character(unlist(md[keep,]['ids']))
c9.fracs <- fracs[c9.ids,]
c9.mat <- mat[, c9.ids]

# Control
keep <- as.character(unlist(md['mutated_gene'])) %in% c('control')
control.ids <- as.character(unlist(md[keep,]['ids']))
control.fracs <- fracs[control.ids,]
control.mat <- mat[, control.ids]



###
# Function to create a dataframe that can be used
# for heatmap plotting
###
makeHeatmapDataFrame <- function(gene_list, ftd.mat, ftd.fracs, identifier = "ensembl_gene_id"){
  if (identifier != "hgnc_symbol"){
    # get HGNC symbols
    bm <- getBM(attributes = c(identifier, "hgnc_symbol"), filters = identifier, values=gene_list, mart=ensembl)
    gene_list <- bm$hgnc_symbol
  }
  
  # Get MAPT expression
  mapt_list <- list()
  for (g in gene_list){
    exp <- as.numeric(ftd.mat[g,])
    res <- nnls(ftd.fracs, exp)
    
    mapt_list$tmp <- res$x
    names(mapt_list)[length(mapt_list)] <- g
  }
  df.mapt <- data.frame(mapt_list)
  rownames(df.mapt) <- paste(colnames(ftd.fracs), "FTD", sep="_")
  
  
  # Get Control expression
  control_list <- list()
  for (g in gene_list){
    exp <- as.numeric(control.mat[g,])
    res <- nnls(control.fracs, exp)
    
    control_list$tmp <- res$x
    names(control_list)[length(control_list)] <- g
  }
  df.control <- data.frame(control_list)
  rownames(df.control) <- paste(colnames(control.fracs),"NDC", sep="_")
  
  # create common dataframe
  df <- rbind(df.mapt, df.control)
  df <- df[order(rownames(df)),]
  
  return(df)

}
###
# Same as above with 'lm' instead of 'nnls
makeHeatmapDataFrameLM <- function(gene_list, identifier = "ensembl_gene_id"){
  if (identifier != "hgnc_symbol"){
    # get HGNC symbols
    bm <- getBM(attributes = c(identifier, "hgnc_symbol"), filters = identifier, values=gene_list, mart=ensembl)
    gene_list <- bm$hgnc_symbol
  }
  
  # Get MAPT expression
  mapt_list <- list()
  for (g in gene_list){
    exp <- as.numeric(mapt.mat[g,])
    res <- lm(mapt.fracs ~ exp -1)
    
    mapt_list$tmp <- as.numeric(res$coefficients)
    names(mapt_list)[length(mapt_list)] <- g
  }
  df.mapt <- data.frame(mapt_list)
  rownames(df.mapt) <- paste(colnames(mapt.fracs), "FTD", sep="_")
  
  
  # Get Control expression
  control_list <- list()
  for (g in gene_list){
    exp <- as.numeric(control.mat[g,])
    res <- lm(control.fracs ~ exp -1)
    
    control_list$tmp <- as.numeric(res$coefficients)
    names(control_list)[length(control_list)] <- g
  }
  df.control <- data.frame(control_list)
  rownames(df.control) <- paste(colnames(control.fracs),"NDC", sep="_")
  
  # create common dataframe
  df <- rbind(df.mapt, df.control)
  df <- df[order(rownames(df)),]
  
  return(df)
  
}



ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
####
# Generate heatmaps for de-regulated genesets
###

# Oxidative phosphorylation
gene_list <- read.table("investigative_genelists/Oxidative phosphorylation_intersection_genes.txt", header=T)
gene_list <- as.character(gene_list$x)
df <- makeHeatmapDataFrame(gene_list)
pheatmap(df, scale = "column", treeheight_row = 0, treeheight_col = 0,
         color = viridis(200, option = "D"), cluster_rows = F,
         filename = "oxphos_pathway_nnls.png", height = 10, width = 16)


# ECM
gene_list <- read.table("investigative_genelists/extracellular matrix organization_intersection_genes.txt", header=T)
gene_list <- as.character(gene_list$x)
df <- makeHeatmapDataFrameLM(gene_list)
pheatmap(df, scale = "column", treeheight_row = 0, treeheight_col = 0,
         color = viridis(200, option = "D"), cluster_rows = F,
         filename = "ecm_pathway_LM.png", height = 10, width = 16)


# vesicle
gene_list <- read.table("investigative_genelists/vesicle_intersection_genes.txt", header=T)
gene_list <- as.character(gene_list$x)
df <- makeHeatmapDataFrameLM(gene_list)
pheatmap(df, scale = "column", treeheight_row = 0, treeheight_col = 0,
         color = viridis(200, option = "D"), cluster_rows = F,
         filename = "vesicle_pathway_LM.png", height = 10, width = 16)

# membrane trafficking
gene_list <- read.table("investigative_genelists/Membrane Trafficking_intersection_genes.txt", header=T)
gene_list <- as.character(gene_list$x)
df <- makeHeatmapDataFrameLM(gene_list, identifier = "refseq_mrna")
pheatmap(df, scale = "column", treeheight_row = 0, treeheight_col = 0,
         color = viridis(200, option = "D"), cluster_rows = F,
         filename = "membrane_trafficking_pathway_LM.png", height = 10, width = 16)

# mir-19b
gene_list <- read.table("investigative_genelists/hsa-miR-19b-3p_intersection_genes.txt", header=T)
gene_list <- as.character(gene_list$x)
df <- makeHeatmapDataFrameLM(gene_list, identifier = "refseq_mrna")
pheatmap(df, scale = "column", treeheight_row = 0, treeheight_col = 0,
         color = viridis(200, option = "D"), cluster_rows = F,
         filename = "mir19b_pathway_LM.png", height = 10, width = 16)


# Gene list testing

gene_list = c("MAPT", "GRN", "C9orf72")
df.mapt <- makeHeatmapDataFrame(gene_list, mapt.mat, mapt.fracs, identifier = "hgnc_symbol")
df.grn <- makeHeatmapDataFrame(gene_list, grn.mat, grn.fracs, identifier = "hgnc_symbol")
df.c9 <- makeHeatmapDataFrame(gene_list, c9.mat, c9.fracs, identifier = "hgnc_symbol")


pheatmap(df, scale = "column", treeheight_row = 0, treeheight_col = 0,
         color = viridis(200, option = "D"), cluster_rows = F)

