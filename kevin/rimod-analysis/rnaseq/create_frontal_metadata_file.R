###########################################################
# Generat metadata file for draft upload of RiMod data
##########################################################
library(tximport)
library(DESeq2)
library(GenomicFeatures)
library(stringr)
library(pheatmap)
setwd("~/rimod/RNAseq/")

# parameters parsing
row_sum_cutoff = 10
row_sum_samples_nr = 5
metadata = "/home/kevin/rimod/files/FTD_Brain.csv"
analysis_dir = "/home/kevin/rimod/RNAseq/analysis/"
region <- "fro"
salmon_files = "/home/kevin/rimod/RNAseq/results_salmon/salmon/"
#====================================================================#


###
# Load data as TXI object
samples <- list.files(salmon_files)
files <- file.path(salmon_files, samples, "quant.sf")
names(files) <- samples
# create txdb
txdb <- makeTxDbFromGFF("~/resources/gencode.v31.annotation.gff3")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
# load counts
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar=TRUE)
#============================#



#== Save TXI stuff ==#
abundance <- txi$abundance
counts <- txi$counts
length <- txi$length


#======================#


# Load metadata
md <- read.csv(metadata, stringsAsFactors = FALSE)
md$SAMPLEID <- as.character(sapply(md$SAMPLEID, function(x){strsplit(x, split="_")[[1]][[1]]}))
md$SAMPLEID <- str_pad(md$SAMPLEID, width = 5, side = "left", pad = "0") # fill sample ids to 5 digits


# bring counts and md in similar format
cts <- txi$counts
rna.samples <- as.character(sapply(colnames(cts), function(x){strsplit(x, split="_")[[1]][[1]]}))
rna.samples <- str_pad(gsub("X", "", rna.samples), width=5, side='left', pad='0')
md <- md[md$SAMPLEID %in% rna.samples,]
md <- md[match(rna.samples, md$SAMPLEID),]

# Generate sample template
temp <- data.frame('Disease Code'=md$DISEASE.CODE, Age=md$AGE, Gender=md$GENDER, Mutation=md$GENE, SampleID=md$SAMPLEID)
temp$Disease.Code <- gsub("control", "Control", temp$Disease.Code)

write.table(temp, "rimod_frontal_rnaseq_metadata.txt", sep="\t", quote=F, row.names = F)
