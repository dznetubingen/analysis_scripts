############################################
# Analysis of RiMod RNA-seq data (frontal) #
############################################


##########################
## PARAMETER SECTION

# load libs
library(DESeq2)
library(stringr)
library(viridis)
library(pheatmap)
library(biomaRt)
library(fgsea)


### Hard-coded section
script_name = "rnaseq_exploratory_basic_analysis.R"
date = Sys.Date()
current_time = gsub(":", ".", gsub(" ", "_", Sys.time()))
####

# parameters parsing
row_sum_cutoff = 1
metadata = "/home/kevin/rimod/files/FTD_Brain.csv"
count_file = "/home/kevin/rimod/RNAseq/RiMOD_RNAseq_frontal_count.matrix.txt"
analysis_dir = "/home/kevin/rimod/RNAseq/analysis/"
region <- "fro"

###########################

# set working directory
setwd(analysis_dir)

# Create sub-folder for current analysis
dir.create(paste("RNAseq_analysis","_",region, "_", current_time, sep=""))
setwd(paste("RNAseq_analysis", "_", region, "_",current_time, sep=""))

# Save parameters in config file
params <- c(current_time, as.character(row_sum_cutoff), metadata, count_file, analysis_dir, script_name, region)
param.names <- c("Time", "Row_sum_cutoff", "Metadata", "Count_file", "Analysis_directory", "Script_name", "Region")
params <- data.frame(param.name = param.names, param = params)
write.table(params, paste("config_file", current_time, sep="_"), quote = F, row.names = F, sep="\t")


# Read cage genewise count table (created by Tenzin)
cts <- read.table(count_file, sep="\t", header=T, row.names = 1, stringsAsFactors = F)
# TODO CHECK if this is valid
cts <- round(cts)

# Load metadata
md <- read.csv(metadata, stringsAsFactors = FALSE)
md$SAMPLEID <- as.character(sapply(md$SAMPLEID, function(x){strsplit(x, split="_")[[1]][[1]]}))
md$SAMPLEID <- str_pad(md$SAMPLEID, width = 5, side = "left", pad = "0") # fill sample ids to 5 digits


# bring counts and md in similar format
rna.samples <- as.character(sapply(colnames(cts), function(x){strsplit(x, split="_")[[1]][[1]]}))
rna.samples <- str_pad(gsub("X", "", rna.samples), width=5, side='left', pad='0')
md <- md[md$SAMPLEID %in% rna.samples,]
md <- md[match(rna.samples, md$SAMPLEID),]


#### REMOVE ALL SPORADIC CASES ####
disease.codes <- c("FTD-C9", "FTD-MAPT", "FTD-GRN", "control")
keep <- md$DISEASE.CODE %in% disease.codes
md <- md[keep,]
cts <- cts[,keep]
md$DISEASE.CODE <- gsub("-", "_", md$DISEASE.CODE) # make disease code names safe
# Split Age covariate into bins
age_bins = 5
md$AGE.BIN <- make.names(cut(md$AGE, breaks=age_bins))

#===========================================#
# DESeq2 analysis
# Generate DDS object
dds <- DESeqDataSetFromMatrix(cts,
                              colData = md,
                              design = ~ AGE.BIN + GENDER + DISEASE.CODE)

# Specify control group
dds$DISEASE.CODE <- relevel(dds$DISEASE.CODE, ref = "control")

# apply prefiltering
keep <- rowSums(counts(dds)) >= row_sum_cutoff
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds)
resnames <- resultsNames(dds)

#== Extract results ==#
### MAPT - control
res.mapt <- results(dds, c("DISEASE.CODE", "FTD_MAPT", "control"))
res.mapt <- na.omit(res.mapt)
deg.mapt <- res.mapt[res.mapt$padj <= 0.05,]

### GRN - control
res.grn <- results(dds, c("DISEASE.CODE", "FTD_GRN", "control"))
res.grn <- na.omit(res.grn)
deg.grn <- res.grn[res.grn$padj <= 0.05,]
### C9orf72 - control
res.c9 <- results(dds, c("DISEASE.CODE", "FTD_C9", "control"))
res.c9 <- na.omit(res.c9)
deg.c9 <- res.c9[res.c9$padj <= 0.05,]

###########
## Save results

write.table(res.mapt, paste("deseq_result_mapt.ndc", "_", region, "_",current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(res.grn, paste("deseq_result_grn.ndc",  "_", region, "_", current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(res.c9, paste("deseq_result_c9.ndc", "_", region, "_",current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)



########################################
## Generate count table and rLog table
########################################

# normalized count values
norm.counts <- counts(dds, normalized=TRUE)
write.table(norm.counts, paste("deseq_normalized_counts", "_", current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)

# reg log transformed values
rld <- rlog(dds, blind=FALSE)
rld.mat <- assay(rld)
write.table(rld.mat, paste("deseq_rLog_values","_", current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)

################################
## Plotting section ############

## PCA
pca <- plotPCA(rld, intgroup = "DISEASE.CODE")
png(paste("pca_group_deseq_rLogvals", "_", current_time, ".png", sep=""), width = 1200, height = 900)
pca
dev.off()

# Make more PCAs
# separate PCAs for the different mutation
library(factoextra)
dc <- md$DISEASE.CODE
mapt.rld <- rld.mat[,dc %in% c('control', 'FTD_MAPT')]
grn.rld <- rld.mat[,dc %in% c('control', 'FTD_GRN')]
c9.rld <- rld.mat[,dc %in% c('control', 'FTD_C9')]


# MAPT - control
mapt.pca <- prcomp(t(mapt.rld), retx=T)
mapt.dc <- dc[dc %in% c('control', 'FTD_MAPT')]
mapt.gene <- md$GENE[dc %in% c('control', 'FTD_MAPT')]
mapt.gene[mapt.dc == 'control'] <- 'control'
fviz_eig(mapt.pca)
mapt.x <- as.data.frame(mapt.pca$x)
mapt.x$Disease_code <- mapt.gene
mpca <- ggplot(mapt.x, aes(x=PC1, y=PC2, color=Disease_code)) +
  geom_point(size=3) +
  stat_ellipse()
png(paste("pca_mapt_rlog", "_", current_time, ".png", sep=""), width = 1200, height = 900)
mpca
dev.off()


# GRN - control
grn.pca <- prcomp(t(grn.rld), retx=T)
grn.dc <- dc[dc %in% c('control', 'FTD_GRN')]
fviz_eig(grn.pca)
grn.x <- as.data.frame(grn.pca$x)
grn.x$Disease_code <- grn.dc
gpca <- ggplot(grn.x, aes(x=PC1, y=PC2, color=Disease_code)) +
  geom_point(size=3) +
  stat_ellipse()
png(paste("pca_grn_rlog", "_", current_time, ".png", sep=""), width = 1200, height = 900)
gpca
dev.off()

# C9 - control
c9.pca <- prcomp(t(c9.rld), retx=T)
c9.dc <- dc[dc %in% c('control', 'FTD_C9')]
fviz_eig(c9.pca)
c9.x <- as.data.frame(c9.pca$x)
c9.x$Disease_code <- c9.dc
cpca <- ggplot(c9.x, aes(x=PC1, y=PC2, color=Disease_code)) +
  geom_point(size=3) +
  stat_ellipse()
png(paste("pca_c9orf72_rlog", "_", current_time, ".png", sep=""), width = 1200, height = 900)
cpca
dev.off()


# PCA for all samples
all.pca <- prcomp(t(rld.mat), retx = T)
fviz_eig(all.pca)
all.x <- as.data.frame(all.pca$x)
all.x$Disease_code <- dc
pca <- ggplot(all.x, aes(x=PC1, y=PC2, color=Disease_code)) +
  geom_point(size=3) 
pca


#####################################
## fGSEA analysis
#####################################

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
pathways <- gmtPathways("~/resources/genesets/h.all.v6.1.entrez.gmt")
pval_filter <- 0.05

## MAPT FGSEA
mapt <- as.data.frame(res.mapt)
bm <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = "hgnc_symbol", values = rownames(mapt), mart = ensembl)
mapt <- merge(mapt, bm, by.x='row.names', by.y='hgnc_symbol')
mapt <- mapt[order(mapt$log2FoldChange),]
ranks <- mapt[,3]
names(ranks) <- mapt$entrezgene
mapt.gsea <- fgsea(pathways, ranks, minSize=15, maxSize=500, nperm=1000)
mapt.gsea <- mapt.gsea[order(mapt.gsea$pval)]
mapt.gsea <- as.data.frame(mapt.gsea)
mapt.gsea <- mapt.gsea[, -ncol(mapt.gsea)] # get rid of last column
write.table(mapt.gsea, "fGSEA_results_hallmark_MAPT.txt", sep="\t", quote=F)


## GRN FGSEA
grn <- as.data.frame(res.grn)
bm <- getBM(attributes = c("hgnc_symbol", "entrezgene"), filters = "hgnc_symbol", values = rownames(grn), mart = ensembl)
grn <- merge(grn, bm, by.x='row.names', by.y='hgnc_symbol')
grn <- grn[order(grn$log2FoldChange),]
ranks <- grn[,3]
names(ranks) <- grn$entrezgene
grn.gsea <- fgsea(pathways, ranks, minSize=15, maxSize=500, nperm=1000)
grn.gsea <- grn.gsea[order(grn.gsea$pval)]
grn.gsea <- as.data.frame(grn.gsea)
grn.gsea <- grn.gsea[, -ncol(grn.gsea)] # get rid of last column
write.table(grn.gsea, "fGSEA_results_hallmark_GRN.txt", sep="\t", quote=F)

## C9 FGSEA
c9 <- as.data.frame(res.c9)
bm <- getBM(attributes = c("hgnc_symbol", "entrezgene"), filters = "hgnc_symbol", values = rownames(c9), mart = ensembl)
c9 <- merge(c9, bm, by.x='row.names', by.y='hgnc_symbol')
c9 <- c9[order(c9$log2FoldChange),]
ranks <- c9[,3]
names(ranks) <- c9$entrezgene
c9.gsea <- fgsea(pathways, ranks, minSize=15, maxSize=500, nperm=1000)
c9.gsea <- c9.gsea[order(c9.gsea$pval),]
c9.gsea <- as.data.frame(c9.gsea)
c9.gsea <- c9.gsea[, -ncol(c9.gsea)] # get rid of last column
write.table(c9.gsea, "fGSEA_results_hallmark_c9orf72.txt", sep="\t", quote=F)





