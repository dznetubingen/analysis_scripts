###########################################################
# Analysis of Salmon quantified RiMod frontal RNA-seq data
# Usin SVA instead of using AGE and Gender as covariates
##########################################################
library(tximport)
library(DESeq2)
library(GenomicFeatures)
library(stringr)
library(pheatmap)
library(IHW)
library(biomaRt)
library(fgsea)
library(sva)
setwd("~/rimod/RNAseq/results_salmon/")


#### Hard-coded section
script_name = "rnaseq_salmon_analysis_rimod_frontal.R"
date = Sys.Date()
current_time = gsub(":", ".", gsub(" ", "_", Sys.time()))
####

# parameters parsing
row_sum_cutoff = 10
row_sum_samples_nr = 5
metadata = "/home/kevin/rimod/files/FTD_Brain.csv"
analysis_dir = "/home/kevin/rimod/RNAseq/analysis/"
region <- "fro"
salmon_files = "/home/kevin/rimod/RNAseq/analysis/txi_salmon/frontal_lengthScaledTPM_counts.txt"
#====================================================================#



# set working directory
setwd(analysis_dir)

# Create sub-folder for current analysis
dir.create(paste("RNAseq_analysis","_",region, "_", current_time, sep=""))
setwd(paste("RNAseq_analysis", "_", region, "_",current_time, sep=""))

# Save parameters in config file
params <- c(current_time, as.character(row_sum_cutoff), metadata, salmon_files, analysis_dir, script_name, region)
param.names <- c("Time", "Row_sum_cutoff", "Metadata", "salmon_files", "Analysis_directory", "Script_name", "Region")
params <- data.frame(param.name = param.names, param = params)
write.table(params, paste("config_file", current_time, sep="_"), quote = F, row.names = F, sep="\t")

#=================================================================#

# Load metadata
md <- read.csv(metadata, stringsAsFactors = FALSE)
md$SAMPLEID <- as.character(sapply(md$SAMPLEID, function(x){strsplit(x, split="_")[[1]][[1]]}))
md$SAMPLEID <- str_pad(md$SAMPLEID, width = 5, side = "left", pad = "0") # fill sample ids to 5 digits

# load counts
cts <- read.table(salmon_files, sep="\t", header=T, row.names=1)

# bring counts and md in similar format
rna.samples <- as.character(sapply(colnames(cts), function(x){strsplit(x, split="_")[[1]][[1]]}))
rna.samples <- str_pad(gsub("X", "", rna.samples), width=5, side='left', pad='0')
md <- md[md$SAMPLEID %in% rna.samples,]
md <- md[match(rna.samples, md$SAMPLEID),]

#### REMOVE ALL SPORADIC CASES ####
disease.codes <- c("FTD-C9", "FTD-MAPT", "FTD-GRN", "control")
keep <- md$DISEASE.CODE %in% disease.codes
md <- md[keep,]
# subset TXI
cts <- cts[,keep]

md$DISEASE.CODE <- gsub("-", "_", md$DISEASE.CODE) # make disease code names safe

###
# SUBSET TO ONLY CONTAIN MAPT STUFF FOR A MINUTE
###
#keep <- md$DISEASE.CODE %in% c("FTD_MAPT", "control")
#cts <- cts[,keep]
#md <- md[keep,]

##
# Calculate surrogate variables
##
library(edgeR)
edata <- cpm(cts)
keep <- rowSums(edata >= 1) >= row_sum_samples_nr
edata <- as.matrix(edata[keep,])
#edata <- log2(edata + 1)
mod1 <- model.matrix(~md$CASE.CONTROL)
mod0 <- cbind(mod1[,1])
svs <- svaseq(edata, mod1, mod0)$sv

colnames(svs) <- paste0("SV", c(1,2,3,4,5,6,7,8))
md <- cbind(md, svs)


#===========================================#
# DESeq2 analysis
# Generate DDS object
cts <- round(cts) # round to integer counts
dds <- DESeqDataSetFromMatrix(cts,
                              colData = md,
                              design = ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + DISEASE.CODE)


# Save DDS object
#saveRDS(dds, file = "frontal_dds_object.rds")

# Specify control group
dds$DISEASE.CODE <- relevel(dds$DISEASE.CODE, ref = "control")

# apply prefiltering
dds <- estimateSizeFactors(dds)
keep <- rowSums((counts(dds, normalized=TRUE) >= row_sum_cutoff)) >= row_sum_samples_nr
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds)
resnames <- resultsNames(dds)

#== Extract results ==#
pval_cut <- 0.05
### MAPT - control
res.mapt <- results(dds, c("DISEASE.CODE", "FTD_MAPT", "control"), filterFun = ihw)
res.mapt <- na.omit(res.mapt)
rownames(res.mapt) <- str_split(rownames(res.mapt), pattern="[.]", simplify = T)[,1]
deg.mapt <- res.mapt[res.mapt$padj <= pval_cut,]
print(dim(deg.mapt))

### GRN - control
res.grn <- results(dds, c("DISEASE.CODE", "FTD_GRN", "control"), filterFun = ihw)
res.grn <- na.omit(res.grn)
rownames(res.grn) <- str_split(rownames(res.grn), pattern="[.]", simplify = T)[,1]
deg.grn <- res.grn[res.grn$padj <= pval_cut,]
print(dim(deg.grn))

### C9orf72 - control
res.c9 <- results(dds, c("DISEASE.CODE", "FTD_C9", "control"), filterFun = ihw)
res.c9 <- na.omit(res.c9)
rownames(res.c9) <- str_split(rownames(res.c9), pattern="[.]", simplify = T)[,1]
deg.c9 <- res.c9[res.c9$padj <= pval_cut,]
print(dim(deg.c9))

