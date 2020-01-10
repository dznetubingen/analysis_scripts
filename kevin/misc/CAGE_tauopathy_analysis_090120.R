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
region <- "tem"

###########################

# set working directory
setwd(analysis_dir)
# Create sub-folder for current analysis
dir.create(paste("CAGE_deseq_analysis","_tauopathy",region, sep=""))
setwd(paste("CAGE_deseq_analysis", "_tauopathy", region, sep=""))

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


#####################################
library(stringr)
pca <- plotPCA(rld, intgroup = "DISEASE.CODE", returnData = TRUE)
pca$sample <- str_split(pca$name, pattern="_", simplify = T)[,2] 
ggplot(pca, aes(x=PC1, y=PC2, color=DISEASE.CODE)) + 
  geom_text(label=pca$sample)



