##################################################
### Analysis of CAGE-seq data                 ####
##################################################


##########################
## PARAMETER SECTION

# load libs
library(DESeq2)
library(stringr)
library(viridis)
library(pheatmap)
library(fgsea)
### Hard-coded section
script_name = "CAGE_deseq_differential_expression_analysis_v1.0.R"
date = Sys.Date()
current_time = gsub(":", ".", gsub(" ", "_", Sys.time()))
####

# parameters parsing
row_sum_cutoff = 1
metadata = "/home/kevin/rimod/files/FTD_Brain.csv"
count_file = "/home/kevin/rimod/CAGE/cage_data/cage_all7regions_3kbgr_aggr.txt"
analysis_dir = "/home/kevin/rimod/CAGE/cage_analysis/"
region <- "fro"

###########################

# set working directory
setwd(analysis_dir)
# Create sub-folder for current analysis
dir.create(paste("CAGE_deseq_analysis","_",region, "_", current_time, sep=""))
setwd(paste("CAGE_deseq_analysis", "_", region, "_",current_time, sep=""))

# Save parameters in config file
params <- c(current_time, as.character(row_sum_cutoff), metadata, count_file, analysis_dir, script_name, region)
param.names <- c("Time", "Row_sum_cutoff", "Metadata", "Count_file", "Analysis_directory", "Script_name", "Region")
params <- data.frame(param.name = param.names, param = params)
write.table(params, paste("config_file", current_time, sep="_"), quote = F, row.names = F, sep="\t")


# Read cage genewise count table (created by Tenzin)
cage <- read.table(count_file, sep="\t", header=T, row.names = 1, stringsAsFactors = F)
# TODO: think about how to remove this hard-coded part ...
# Remove sample sample_09218_froR
# Keep only desired region
cage <- cage[,grepl(region, colnames(cage))]
# Remove 'froR' samples
cage <- cage[,!grepl("froR", colnames(cage))]


# Load metadata
md <- read.csv(metadata, stringsAsFactors = FALSE)
md$SAMPLEID <- str_pad(md$SAMPLEID, width = 5, side = "left", pad = "0") # fill sample ids to 5 digits

# bring counts and md in similar format
cage.samples <- as.character(gsub("sample_","",colnames(cage)))
cage.samples <- as.character(sapply(cage.samples, function(x){strsplit(x, split=paste("_", region, sep=""))[[1]][[1]]}))
md <- md[md$SAMPLEID %in% cage.samples,]
md <- md[match(cage.samples, md$SAMPLEID),]


#### REMOVE ALL SPORADIC CASES ####
disease.codes <- c("FTD-C9", "FTD-MAPT", "FTD-GRN", "control")
keep <- md$DISEASE.CODE %in% disease.codes
md <- md[keep,]
cage <- cage[,keep]
md$DISEASE.CODE <- gsub("-", "_", md$DISEASE.CODE) # make disease code names safe

# PH
ph <- as.numeric(md$PH)
ph.mean <- mean(na.omit(ph))
ph[is.na(ph)] <- ph.mean
md$PH <- ph

rownames(md) <- colnames(cage)
#===========================================#
# DESeq2 analysis
# Generate DDS object
dds <- DESeqDataSetFromMatrix(cage,
                              colData = md,
                              design = ~ PH +  GENDER + DISEASE.CODE)

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
pval_cut <- 0.05
res.mapt <- results(dds, c("DISEASE.CODE", "FTD_MAPT", "control"))
res.mapt <- na.omit(res.mapt)
deg.mapt <- res.mapt[res.mapt$padj <=pval_cut,]
print(dim(deg.mapt))
### GRN - control
res.grn <- results(dds, c("DISEASE.CODE", "FTD_GRN", "control"))
res.grn <- na.omit(res.grn)
deg.grn <- res.grn[res.grn$padj <= pval_cut,]
print(dim(deg.grn))

### C9orf72 - control
res.c9 <- results(dds, c("DISEASE.CODE", "FTD_C9", "control"))
res.c9 <- na.omit(res.c9)
deg.c9 <- res.c9[res.c9$padj <= pval_cut,]
print(dim(deg.c9))

###########
## Save results

write.table(res.mapt, paste("deseq_result_mapt.ndc", "_",current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(res.grn, paste("deseq_result_grn.ndc",  "_", current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(res.c9, paste("deseq_result_c9.ndc", "_", current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)

# Save only significant genes for online tools
write.table(rownames(deg.mapt), paste("DEGs_mapt.ndc", "_", region, "_",current_time, ".txt", sep=""), sep="\t", quote=F, row.names=F)
write.table(rownames(deg.grn), paste("DEGs_grn.ndc", "_", region, "_",current_time, ".txt", sep=""), sep="\t", quote=F, row.names=F)
write.table(rownames(deg.c9), paste("DEGs_c9.ndc", "_", region, "_",current_time, ".txt", sep=""), sep="\t", quote=F, row.names=F)

# Divde in up and down regulated genes
# MAPT
mapt.up <- deg.mapt[deg.mapt$log2FoldChange > 0,]
mapt.down <- deg.mapt[deg.mapt$log2FoldChange < 0,]
write.table(rownames(mapt.up), "DEGs_UP_mapt.ndc.txt", quote=F, row.names=F)
write.table(rownames(mapt.down), "DEGs_Down_mapt.ndc.txt", quote=F, row.names=F)
# GRN
grn.up <- deg.grn[deg.grn$log2FoldChange > 0,]
grn.down <- deg.grn[deg.grn$log2FoldChange < 0,]
write.table(rownames(grn.up), "DEGs_UP_grn.ndc.txt", quote=F, row.names=F)
write.table(rownames(grn.down), "DEGs_Down_grn.ndc.txt", quote=F, row.names=F)
# C9orf72
c9.up <- deg.c9[deg.c9$log2FoldChange > 0,]
c9.down <- deg.c9[deg.c9$log2FoldChange < 0,]
write.table(rownames(c9.up), "DEGs_UP_c9.ndc.txt", quote=F, row.names=F)
write.table(rownames(c9.down), "DEGs_Down_c9.ndc.txt", quote=F, row.names=F)

########################################
## Generate count table and rLog table
########################################

# normalized count values
norm.counts <- counts(dds, normalized=TRUE)
write.table(norm.counts, paste("deseq_normalized_counts", "_", current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)

# reg log transformed values
rld <- vst(dds, blind=FALSE)
rld.mat <- assay(rld)
write.table(rld.mat, paste("deseq_vst_values","_", current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)

################################
## Plotting section ############

## PCA
pca <- plotPCA(rld, intgroup = "DISEASE.CODE")
png(paste("pca_group_deseq_rLogvals", "_", current_time, ".png", sep=""), width = 1200, height = 900)
pca
dev.off()

## Make more PCAs
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
  geom_point(size=3) +
  stat_ellipse()
pca

#####################################
## fGSEA analysis
#####################################
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
pathways <- gmtPathways("~/resources/genesets/h.all.v6.1.entrez.gmt")
pval_filter <- 0.05

## MAPT FGSEA
mapt <- as.data.frame(res.mapt)
bm <- getBM(attributes = c("ensembl_gene_id", "entrezgene"), filters = "ensembl_gene_id", values = rownames(mapt), mart = ensembl)
mapt <- merge(mapt, bm, by.x='row.names', by.y='ensembl_gene_id')
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
bm <- getBM(attributes = c("ensembl_gene_id", "entrezgene"), filters = "ensembl_gene_id", values = rownames(grn), mart = ensembl)
grn <- merge(grn, bm, by.x='row.names', by.y='ensembl_gene_id')
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
bm <- getBM(attributes = c("ensembl_gene_id", "entrezgene"), filters = "ensembl_gene_id", values = rownames(c9), mart = ensembl)
c9 <- merge(c9, bm, by.x='row.names', by.y='ensembl_gene_id')
c9 <- c9[order(c9$log2FoldChange),]
ranks <- c9[,3]
names(ranks) <- c9$entrezgene
c9.gsea <- fgsea(pathways, ranks, minSize=15, maxSize=500, nperm=1000)
c9.gsea <- c9.gsea[order(c9.gsea$pval),]
c9.gsea <- as.data.frame(c9.gsea)
c9.gsea <- c9.gsea[, -ncol(c9.gsea)] # get rid of last column
write.table(c9.gsea, "fGSEA_results_hallmark_c9orf72.txt", sep="\t", quote=F)



#== Prioritize Genes for Humanbase Use ==#
# Save only significant genes for online tools
hb.mapt <- deg.mapt[abs(deg.mapt$log2FoldChange) > 0.8,]
hb.grn <- deg.grn[abs(deg.grn$log2FoldChange) > 0.8,]
hb.c9 <- deg.c9[abs(deg.c9$log2FoldChange) > 0.8,]

write.table(rownames(hb.mapt), paste("DEGs_HB_lfc0.8_mapt.ndc", "_", region, "_",current_time, ".txt", sep=""), sep="\t", quote=F, row.names=F)
write.table(rownames(hb.grn), paste("DEGs_HB_lfc0.8_grn.ndc", "_", region, "_",current_time, ".txt", sep=""), sep="\t", quote=F, row.names=F)
write.table(rownames(hb.c9), paste("DEGs_HB_lfc0.8_c9.ndc", "_", region, "_",current_time, ".txt", sep=""), sep="\t", quote=F, row.names=F)



