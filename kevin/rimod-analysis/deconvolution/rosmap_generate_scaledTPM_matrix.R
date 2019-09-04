##########################
# Generatin of a lengthScadel TPM matrix 
# from Salmon output for ROSMAP data
##############################

library(tximport)
library(GenomicFeatures)


analysis_dir = "/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/deconvolution/"
salmon_files = "/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/deconvolution/rosmap_data/"
#====================================================================#
setwd(analysis_dir)

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
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")

counts <- txi$counts
setwd(analysis_dir)

#============================#

##
# ID mapping
##
library(biomaRt)
library(stringr)
genes <- rownames(counts)
genes <- str_split(genes, patter="[.]", simplify = T)[,1]
rownames(counts) <- genes

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=genes, mart=ensembl)

counts <- merge(counts, bm, by.x="row.names", by.y="ensembl_gene_id")

symbols <- as.character(counts$hgnc_symbol)
counts <- counts[, c(-1, -ncol(counts))]
counts <- aggregate(counts, by=list(symbols), FUN=sum) # aggregate duplicates
rownames(counts) <- counts$Group.1
counts <- counts[-1,]
counts <- counts[,-1]


write.table(counts, "ROSMAP_lengthScaeldTPM_counts.txt", sep="\t", quote=F, col.names = NA)

