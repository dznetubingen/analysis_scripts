##################################################
### Analysis of CAGE-seq data                 ####
##################################################


##########################
## PARAMETER SECTION

## Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 5) {
  args <- c("--help")
}
## Help section
if("-h" %in% args ||"--help" %in% args) {
  cat("
      Please pass 6 arguments to the script in the following order:
      <row_sum_cutoff> the row sum cutoff for the analysis
      <metadata file> the metadata file (FTD_Brain.csv)
      <count file> the file containing the count data
      <analysis dir> the directory to store the results in (a subdirectory will be created)
      <region> The brain region to use, e.g. 'fro', 'temp'")
  stop()
}

# Edit this section only for adjustment


# load libs
library(DESeq2)
library(stringr)
library(viridis)
library(pheatmap)
#source("~/scripts/utility_funs.R")

### Hard-coded section
script_name = "CAGE_deseq_differential_expression_analysis_v1.0.R"
date = Sys.Date()
current_time = gsub(":", ".", gsub(" ", "_", Sys.time()))
####

# parameters parsing
row_sum_cutoff = args[1]
metadata = args[2]
count_file = args[3]
analysis_dir = args[4]
region <- args[5]

###########################

# set working directory
setwd(analysis_dir)
# Create sub-folder for current analysis
dir.create(paste("CAGE_deseq_analysis","_", current_time, sep=""))
setwd(paste("CAGE_deseq_analysis", "_",current_time, sep=""))

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
# Split Age covariate into bins
age_bins = 5
md$AGE.BIN <- make.names(cut(md$AGE, breaks=age_bins))

#===========================================#
# DESeq2 analysis
# Generate DDS object
dds <- DESeqDataSetFromMatrix(cage,
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

write.table(res.mapt, paste("deseq_result_mapt.ndc", "_",current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(res.grn, paste("deseq_result_grn.ndc",  "_", current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(res.c9, paste("deseq_result_c9.ndc", "_", current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)



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
