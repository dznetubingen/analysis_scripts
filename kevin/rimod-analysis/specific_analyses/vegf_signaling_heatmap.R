############
# Make Heatmap of genes important to vasculature growth
############
library(stringr)
library(biomaRt)
library(pheatmap)
library(viridis)

setwd("~/rimod/integrative_analysis/vegf_signaling/")

# Load Reactome pathway
vegf <- read.table("reactome_vegf_signaling.tsv", sep="\t", header=T)
genes <- as.character(vegf$MoleculeName)
genes <- str_split(genes, pattern=" ", simplify=T)[,2]

# Load gene expression
rna <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_vst_values_2019-10-23_13.33.11.txt", sep="\t", header=T, row.names=1)
rownames(rna) <- str_split(rownames(rna), pattern="[.]", simplify=T)[,1]

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=rownames(rna), mart=ensembl)
bm <- bm[!bm$hgnc_symbol == "",]
rna <- merge(rna, bm, by.x="row.names", by.y="ensembl_gene_id")
rna <- rna[!duplicated(rna$hgnc_symbol),]
rownames(rna) <- rna$hgnc_symbol
rna <- rna[, c(-1, -ncol(rna))]

# format colnames
colnames(rna) <- gsub("X", "", colnames(rna))
colnames(rna) <- str_split(colnames(rna), pattern="_", simplify = T)[,1]
colnames(rna) <- str_pad(colnames(rna), width=5, side="left", pad=0)

# load metadata
md <- read.table("~/rimod/RNAseq/rimod_frontal_rnaseq_metadata.txt", sep="\t", header=T)
md <- md[match(colnames(rna), md$SampleID),]

###
# Make Heatmap
###
rna <- na.omit(rna[genes,])

# make annotation df
mutation <- as.character(md$Mutation)
mutation[!mutation == "P301L"] <- "MAPT-other"
mutation[md$Disease.Code == "FTD-GRN"] <- "GRN"
mutation[md$Disease.Code == "FTD-C9"] <- "C9orf72"
mutation[md$Disease.Code == "Control"] <- "Control"
colanno <- data.frame(mutation)
rownames(colanno) <- colnames(rna)

pheatmap(rna, scale="row", color = viridis(200, option = "A"), annotation_col = colanno)


# Make a PCA
pca <- prcomp(t(rna))
pca <- as.data.frame(pca$x)
pca$mutation <- mutation

library(ggplot2)
ggplot(pca, aes(x=PC1, y=PC2, color=mutation)) + 
  geom_text(label = mutation)





##########
# WGCNA VEGF gene module membership
##########

mms <- read.table("~/rimod/RNAseq/analysis/wgcna_modules/WGCNA_gene_module_membership.txt", sep="\t", header=T)
mms_all <- mms
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=mms$geneID, mart=ensembl)
mms <- merge(mms, bm, by.x="geneID", by.y="ensembl_gene_id")

mms <- mms[mms$hgnc_symbol %in% genes,]
res <- data.frame(table(mms$moduleColors))
res <- res[order(res$Freq, decreasing = T),]

res_all <- data.frame(table(mms_all$moduleColors))
res_all <- res_all[match(res$Var1, res_all$Var1),]

# calculate ratio
res$ratio <- res$Freq / res_all$Freq
