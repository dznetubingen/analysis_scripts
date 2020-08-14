#######################################################################
# Analysis of smRNA-seq in combination with RNA-seq data for better 
# validation of actual miRNA-target regulation
#######################################################################
library(biomaRt)
library(stringr)
setwd("~/rimod/smallRNA/frontal/analysis/")

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="uswest.ensembl.org", ensemblRedirect = FALSE)

####
# Formatting of expression matrices
####

# sRNA expression
srna <- read.table("analysis_0719/deseq_rLog_values_temporal_smRNA.txt", sep="\t", header=T, row.names=1)
short_samples <- str_sub(gsub("sample_", "", gsub("X", "", colnames(srna))), 1, 5)

# mrna expression
mrna <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2020-05-04_15.45.57/deseq_vst_values_2020-05-04_15.45.57.txt", sep="\t", header=T, row.names=1)
long_samples <- str_sub(gsub("X", "", colnames(mrna)), 1, 5)

# Only use intersection of samples
samples <- intersect(short_samples, long_samples)
short_keep <- short_samples %in% samples
long_keep <- long_samples %in% samples
short_samples <- short_samples[short_keep]
long_samples <- long_samples[long_keep]
srna <- srna[,short_keep]
mrna <- mrna[,long_keep]
# bring in correct order
order_long <- match(short_samples, long_samples)
mrna <- mrna[, order_long]
long_samples <- long_samples[order_long]

# get rid of dot
rownames(mrna) <- str_split(rownames(mrna), pattern = "[.]", simplify = T)[,1]

# get refseq IDs for matching with targets
bm <- getBM(attributes = c("ensembl_gene_id", "refseq_mrna"), filters="ensembl_gene_id", mart = ensembl, values = rownames(mrna))

mrna <- merge(mrna, bm, by.x="row.names", by.y="ensembl_gene_id")
mrna <- mrna[!duplicated(mrna$refseq_mrna),]
mrna <- mrna[!is.na(mrna$refseq_mrna),]
rownames(mrna) <- mrna$refseq_mrna
mrna <- mrna[, c(-1, -(ncol(mrna)))]

##===========================================================#

# Correlate expression of hsa-miR-19b-3p with all its' predicted targets
cor_cutoff = -0.4

tgt <- read.table("mirna_target_analysis_0719/miR_19b_3p_targets.txt", sep="\t", header=T, stringsAsFactors = F)
mir <- as.character(tgt$V1)[1]

target_list <- list()
mir.exp <- as.numeric(srna[mir,])
target_vec <- c()


# calculate correlation for each mir-target pair
for (i in 1:nrow(tgt)) {
  gene <- tgt[i,2]
  gene.exp <- as.numeric(mrna[gene,])
  mir.gene.cor <- cor(mir.exp, gene.exp, method='pearson') # calculate correlation 
  
  # check if correlation is above threshold
  if (!is.na(mir.gene.cor)){
    if (mir.gene.cor < cor_cutoff){
      target_vec <- c(target_vec, gene)
    }
  }
}


bm <- getBM(attributes = c("hgnc_symbol", "refseq_mrna"), filters="refseq_mrna", mart = ensembl, values = target_vec)

targets_hgnc <- bm$hgnc_symbol
targets_hgnc <- targets_hgnc[!duplicated(targets_hgnc)]

write.table(target_vec, "miR_19b_correlated_targets_refseq.txt", row.names = F, quote=F, col.names = F)
write.table(targets_hgnc, "miR_19b_correlated_targets_HGNC.txt", row.names = F, quote=F, col.names = F)


