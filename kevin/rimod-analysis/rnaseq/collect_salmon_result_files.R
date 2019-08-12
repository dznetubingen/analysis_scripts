############################################
# Collect results from Salmon experiment
#########################################
library(tximport)
library(DESeq2)
library(GenomicFeatures)
setwd("~/rimod/RNAseq/results_salmon/")

samples <- list.files("salmon/")
files <- file.path("salmon", samples, "quant.sf")
names(files) <- samples

txdb <- makeTxDbFromGFF("~/resources/gencode.v31.annotation.gff3")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar=TRUE)
counts <- txi.salmon$counts


# TODO: adjust to also collect the TPMs
tpm.txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = T, countsFromAbundance = "lengthScaledTPM")
