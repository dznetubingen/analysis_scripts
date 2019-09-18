#######################
# Deconvolution analysis plotting
#######################
library(stringr)
library(biomaRt)
library(pheatmap)
library(viridis)

setwd("~/rimod/RNAseq/analysis/deconvolution/")

# Load expresion values
mat <- read.table("frontal_lengthScaledTPM_counts.txt", sep="\t", row.names=1, header=T)
rownames(mat) <- str_split(rownames(mat), pattern="[.]", simplify = T)[,1]
colnames(mat) <- gsub("X", "", colnames(mat))

# Get HGNC symbols
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=rownames(mat), mart=ensembl)


 # Load deconvolution results
fracs <- read.table("cdn_predictions.txt", sep="\t", row.names = 1, header=T)
rownames(fracs) <- gsub("X", "", rownames(fracs))
fracs$Neurons <- fracs$InNeurons + fracs$ExNeurons
fracs <- fracs[, c(-2, -8)]

# Load MD 
md <- read.table("~/rimod/RNAseq/rnaseq_frontal_md.txt", sep="\t", header=T, stringsAsFactors = F)
md <- md[match(rownames(fracs), md$ids),]
md <- md[-19,]

# keep mapt and control
keep <- as.character(unlist(md['mutated_gene'])) %in% c('MAPT', 'control')
# testing for differential expression
md <- md[keep,]
fracs = fracs[as.character(unlist(md['ids'])),]
mat <- mat[, as.character(unlist(md['ids']))]

# remove unknown from fracs
fracs <- fracs[,-1]
# make as matrix
fracs <- as.matrix(fracs)
###
# testing
###
test_genes <- rownames(mat)[1:5]

res_list <- list()
for (g in test_genes){
  exp <- as.numeric(mat[g,])
  res <- nnls(fracs, exp)
  
  res_list$tmp <- res$x
  names(res_list)[length(res_list)] <- g
}
# turn into dataframe
df <- data.frame(res_list)
rownames(df) <- colnames(fracs)
pheatmap(df, scale = "column", treeheight_row = 0, treeheight_col = 0, color = viridis(200, option = "A"))
