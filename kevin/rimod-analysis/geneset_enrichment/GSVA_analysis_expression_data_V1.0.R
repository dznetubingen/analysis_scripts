################################
### GSVA for Expression data ###
################################
# Author: Kevin Menden
# Performs GSVA on normalized expression data
# Used to analyze data from RiMod (CAGE and RNA-seq)

##########################
## PARAMETER SECTION

## Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 3) {
  args <- c("--help")
}
## Help section
if("-h" %in% args ||"--help" %in% args) {
  cat("
      Please pass 6 arguments to the script in the following order:
      <expression_values> the normalized expression values
      <design file> the design file
      <analysis dir> the directory to store the results in (a subdirectory will be created)")
  stop()
}


# Load Libraries
library(edgeR)
library(GSVA)
library(GSEABase)
library(pheatmap)
library(limma)
library(viridis)
library(biomaRt)
library(stringr)

### Hard-coded section
script_name = "GSVA_analysis_expression_data_V1.0.R"
date = Sys.Date()
current_time = gsub(":", ".", gsub(" ", "_", Sys.time()))
gmt_file <- "~/resources/genesets/c2.all.v6.1.symbols.gmt"
####

# Testing

# parameters parsing
expression_values = args[1]
design_file = args[2]
analysis_dir = args[3]

###########################



# set working directory
setwd(analysis_dir)
# Create sub-folder for current analysis
dir.create(paste("GSVA_analysis","_", current_time, sep=""))
setwd(paste("GSVA_analysis", "_",current_time, sep=""))

# Save parameters in config file
params <- c(current_time, expression_values, gmt_file, analysis_dir, script_name)
param.names <- c("Time", "Expression values", "GMT file", "Analysis_directory", "Script_name")
params <- data.frame(param.name = param.names, param = params)
write.table(params, paste("config_file", current_time, sep="_"), quote = F, row.names = F, sep="\t")

# Read rLog expression values generated using DESeq2 (or other normalized values)
# The rLog values are normalized and prefiltered
rna <- read.table(expression_values, sep="\t", header=T, check.names = F, row.names=1)
rownames(rna) <- as.character(sapply(rownames(rna), function(x){strsplit(x, split="[.]")[[1]][[1]]})) # remove ensembl version number

# Use biomaRt to get symbols for gene identifiers
# TODO: work on ENTREZ identifier insteads
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- rownames(rna)
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = genes, mart = ensembl)
keep.genes.hgnc <- !bm$hgnc_symbol == ""
print("Conversion of ENSEMBL IDs to symbols ")
table(keep.genes.hgnc)

mdata <- merge(rna, bm, by.x="row.names", by.y="ensembl_gene_id")
mdata <- mdata[!mdata$hgnc_symbol == "",]
mdata <- mdata[!duplicated(mdata$Row.names),]
rownames(mdata) <- mdata$hgnc_symbol
mdata <- mdata[,c(-1,-ncol(mdata))]
rna <- mdata

# Load design
design <- read.csv(design_file)

# Load hallmark colleation ("H") 
gmt <- getGmt(gmt_file, geneIdType = SymbolIdentifier(), collectionType = BroadCollection(category = "c2"))

# Run GSVA
gsva.res <- gsva(as.matrix(rna), gmt, kcdf = "Gaussian", parallel.sz = 3)


# Create annotation data frame for heatmap plotting 
design$SAMPLEID <- str_pad(design$SAMPLEID, width = 5, side = "left", pad = "0") # fill sample ids to 5 digits
samples <- as.character(sapply(colnames(rna), function(x){strsplit(x, split="_")[[1]][[2]]}))
samples[samples == "A144"] <- "A144_12" # manual adjustment of sample A144 (hard-coded ..)
# bring counts and md in similar format

design <- design[design$SAMPLEID %in% samples,]
design <- design[match(samples, design$SAMPLEID),]
col.annot <- data.frame(age = design$AGE, pmd = design$PMD, gender = design$GENDER, group = design$DISEASE.CODE)
rownames(col.annot) <- colnames(gsva.res)

# Plot the results from GSVA in a heatmap
png(paste("gsva_heatmap_", current_time, "_all_genesets.png", sep=""), width=1500, height = 1000)
pheatmap(gsva.res, show_rownames = F, annotation_col = col.annot, color = viridis(200, option="B"))
dev.off()

write.table(gsva.res, paste("GSVA_result_all_", current_time, ".txt", sep=""), sep="\t", quote =F)

# The advantage of GSVA is that you can perform differential analysis on the GSVA enrichment scores
# using complex comparisons. Using limma one can easily find differentially expressed genesets based
# on der enrichment scores
# Limma analysis on GSVA scores
G <- factor(design$DISEASE.CODE)
dm <- model.matrix(~ -1 + G)
colnames(dm) <- c("control", "C9", "GRN", "MAPT")
cont = c("C9-control", "GRN-control", "MAPT-control", "(MAPT+GRN+C9)-control")
fit = lmFit(gsva.res, design = dm)
cm <- makeContrasts(contrasts = cont, levels = dm)
cont.fit <- contrasts.fit(fit, contrasts = cm)
fit <- eBayes(cont.fit)
res <- decideTests(fit)
summary(res)

# Based on the results from the limma analysis, plot only the differentially enrichred
# genesets for the different comparisons

# Case control
res.all <- topTable(fit, coef=4, number = Inf, p.value = 0.01)
res.mat.all <- gsva.res[rownames(gsva.res) %in% rownames(res.all),]
write.table(res.mat.all,paste("GSVA_FTD.NDC_sig_scores_", current_time, ".txt", sep=""), sep="\t", quote=F)
png(paste("gsva_heatmap_", current_time, "_FTD.NDC_pval001.png", sep=""), width=1500, height = 1000)
pheatmap(res.mat.all, show_rownames = T, annotation_col = col.annot, color = viridis(200, option="B"))
dev.off()


# GRN
res.grn <- topTable(fit, coef="GRN-control", number = Inf, p.value = 0.01)
res.mat <- gsva.res[rownames(gsva.res) %in% rownames(res.grn),]
write.table(res.mat,paste("GSVA_GRN.NDC_sig_scores_", current_time, ".txt", sep=""), sep="\t", quote=F)
png(paste("gsva_heatmap_", current_time, "_FTDGRN.NDC_pval001.png", sep=""), width=1500, height = 1000)
pheatmap(res.mat, show_rownames = T, annotation_col = col.annot, color = viridis(200, option="B"))
dev.off()

# MAPT
res.mapt <- topTable(fit, coef="MAPT-control", number = Inf, p.value = 0.01)
res.mat <- gsva.res[rownames(gsva.res) %in% rownames(res.mapt),]
write.table(res.mat,paste("GSVA_MAPT.NDC_sig_scores_", current_time, ".txt", sep=""), sep="\t", quote=F)
png(paste("gsva_heatmap_",current_time,"_FTDMAPT.NDC_pval001.png", sep=""), width=1500, height = 1000)
pheatmap(res.mat, show_rownames = T, annotation_col = col.annot, color = viridis(200,option = "B"))
dev.off()

# C9orf72
res.c9 <- topTable(fit, coef="C9-control", number = Inf, p.value = 0.01)
res.mat.c9 <- gsva.res[rownames(gsva.res) %in% rownames(res.c9),]
write.table(res.mat,paste("GSVA_C9.NDC_sig_scores_", current_time, ".txt", sep=""), sep="\t", quote=F)
png(paste("gsva_heatmap_", current_time,"_FTDC9.NDC_pval001.png", sep=""), width=1500, height = 1000)
pheatmap(res.mat.c9, show_rownames = T, annotation_col = col.annot, color = viridis(200, option = "B"))
dev.off()

print("Script finished :-)")