###################################################
# Generation of a lengthScaled TPM matrix 
# from Salmon output for use with other packages
###################################################

library(tximport)
library(GenomicFeatures)


#### Hard-coded section
script_name = "rnaseq_salmon_analysis_rimod_frontal.R"
date = Sys.Date()
current_time = gsub(":", ".", gsub(" ", "_", Sys.time()))
####

analysis_dir = "/home/kevin/rimod/RNAseq/analysis/"
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
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar=TRUE, countsFromAbundance = "lengthScaledTPM")

counts <- txi$counts
setwd(analysis_dir)
write.table(counts, "frontal_lengthScaeldTPM_counts.txt", sep="\t", quote=F, col.names = NA)
#============================#