##################################################################################
### Differential expression analysis of smallRNA-seq data using DESeq2
### Done on 23.02.2018
### Comparison: All
##################################################################################
# Author: Kevin Menden

# load libs
library(DESeq2)
library(stringr)
library(viridis)
library(pheatmap)
library(limma)
source("~/scripts/utility_funs.R")

## SCRIPT NAME
# change for new version!
script_name = "deseq_smallRNA_analysis_script_v1.0_230218.R"

##########################
## PARAMETER SECTION

## Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 4) {
  args <- c("--help")
}
## Help section
if("-h" %in% args ||"--help" %in% args) {
  cat("
        Please pass 4 arguments to the script in the following order:
        <row_sum_cutoff> the row sum cutoff for the analysis
        <metadata file> the metadata file
        <count file> the file containing the count data
        <analysis dir> the directory to store the results in (a subdirectory will be created)")
  stop()
}

date = Sys.Date()
current_time = gsub(":", ".", gsub(" ", "_", Sys.time()))
row_sum_cutoff = args[1]
metadata = args[2]
count_file = args[3]
analysis_dir = args[4]

###########################

# go to analysis directory
setwd(analysis_dir)
# Create sub-folder for current analysis
dir.create(paste("deseq_smallRNA_analysis", current_time, sep="_"))
setwd(paste("deseq_smallRNA_analysis", current_time, sep="_"))

# Save parameters in config file
params <- c(current_time, as.character(row_sum_cutoff), metadata, count_file, analysis_dir, script_name)
param.names <- c("Time", "Row_sum_cutoff", "Metadata", "Count_file", "Analysis_directory", "Script_name")
params <- data.frame(param.name = param.names, param = params)
write.table(params, paste("config_file", current_time, sep="_"), quote = F, row.names = F, sep="\t")


# parse metadata
md <- read.table(metadata, sep="\t", header=T)
md$sample <- str_pad(md$sample, 5, side="left", pad="0")
write.table(md, paste("design_file_", current_time, ".txt", sep=""), quote=F, sep="\t")
# Exclude bad samples
md <- md[md$include == 1,]

# parse counts
cts <- read.table(count_file, sep="\t", row.names=1, header=T)
cols <- colnames(cts)
cols <- as.character(gsub("RNAomeTb", "", cols))
cols <- substr(cols, 1,5) # get only samples name
colnames(cts) <- cols

# Bring  metadata and counts in correct order and format
cts <- cts[,colnames(cts) %in% md$sample]
md <- md[md$sample %in% colnames(cts),]
cts <- cts[,match(md$sample, colnames(cts))]
colnames(cts) <- paste("sample", colnames(cts), sep="_")
md$group <- make.names(md$group)
rownames(md) <- colnames(cts)





#### DESeq2 Analysis
# Create DESeqData object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = md,
                              design = ~ age + gender + group)


# apply prefiltering
keep <- rowSums(counts(dds)) >= row_sum_cutoff
dds <- dds[keep,]

# Specify control group
dds$group <- relevel(dds$group, ref = "NDC")

# Run DESeq
dds <- DESeq(dds)
resnames <- resultsNames(dds)

#== Extract results ==#
### MAPT - control
res.mapt <- results(dds, c("group", "FTD.MAPT", "NDC"))
res.mapt <- na.omit(res.mapt)
deg.mapt <- res.mapt[res.mapt$padj <= 0.05,]

### GRN - control
res.grn <- results(dds, c("group", "FTD.GRN", "NDC"))
res.grn <- na.omit(res.grn)
deg.grn <- res.grn[res.grn$padj <= 0.05,]
### C9orf72 - control
res.c9 <- results(dds, c("group", "FTD.C9", "NDC"))
res.c9 <- na.omit(res.c9)
deg.c9 <- res.c9[res.c9$padj <= 0.05,]

###########
## Save results

write.table(res.mapt, paste("deseq_result_mapt.ndc_", current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(res.grn, paste("deseq_result_grn.ndc_", current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(res.c9, paste("deseq_result_c9.ndc_", current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)

########################################
## Generate count table and rLog table
########################################

# normalized count values
norm.counts <- counts(dds, normalized=TRUE)
write.table(norm.counts, paste("deseq_normalized_counts_", current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)

# reg log transformed values
rld <- rlog(dds, blind=FALSE)
rld.mat <- assay(rld)
write.table(rld.mat, paste("deseq_rLog_values", current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)

################################
## Plotting section ############

## PCA
pca <- plotPCA(rld, intgroup = "group")
png(paste("pca_group_deseq_rLogvals", current_time, ".png", sep=""), width = 1200, height = 900)
pca
dev.off()

# PCA with sample names for outlier detection
png(paste("pca_sample_name_rLog_", current_time, ".png", sep=""), width= 1200, height = 900)
plotMDS(rld.mat, labels = md$sample)
dev.off()


