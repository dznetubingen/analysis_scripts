#######################################################################
# Analysis of smRNA-seq in combination with RNA-seq data for better 
# validation of actual miRNA-target regulation
#######################################################################
library(biomaRt)
library(stringr)
setwd("~/rimod/smallRNA/frontal/analysis/")

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

####
# Formatting of expression matrices
####

# sRNA expression
srna <- read.table("analysis_0719/deseq_rLog_values_frontal_smRNA.txt", sep="\t", header=T, row.names=1)
short_samples <- str_sub(gsub("sample_", "", gsub("X", "", colnames(srna))), 1, 5)

# mrna expression
mrna <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-08-12_07.58.35/deseq_vst_values_2019-08-12_07.58.35.txt", sep="\t", header=T, row.names=1)
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



###
# Function to generate a list of targets that are negatively correlated with
# their targeting miRNAs.
###
generateTargetList <- function(deg.srna.path, targets.path, mrna, srna, 
                               pval = 0.05, lfc = 0.8, cor_cutoff = -0.4){
  
  # Cutoffs
  pval <- 0.05
  lfc <- 0.8
  cor_cutoff = -0.4
  
  
  # Load miRNA DEGs
  deg.srna <- read.table(deg.srna.path, sep="\t", header=T)
  deg.srna <- deg.srna[deg.srna$padj <= pval,]
  deg.srna <- deg.srna[abs(deg.srna$log2FoldChange) >= lfc,]
  
  # Load targets
  targets <- read.table(targets.path, sep="\t", header=T, stringsAsFactors = F)
  
  ## Up-regulated miRNAs
  target_list <- list()

  count = 0
  for (mir in as.character(deg.srna$X)) {
    # Some printing while processing
    print(mir)
    count = count + 1
    print(count)
    # Get the predicted targets for this miRNA
    tgt <- targets[targets$V1 == mir,]
    tgt <- tgt[tgt$V2 %in% rownames(mrna),]
    mir.exp <- as.numeric(srna[mir,]) # expression of this miRNA
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
    # Add targetes to list
    target_list$tmp <- target_vec
    names(target_list)[length(target_list)] <- mir
  }
  
  return(target_list)
}
#======00===================================================================#

###############
# Calculate the negatively correlated targets
###############
setwd("/home/kevin/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/")
# Cutoffs

pval <- 0.05
lfc <- 0.8
cor_cutoff = -0.4

## MAPT
mapt_targets <- generateTargetList(deg.srna.path = "../analysis_0719/deseq_result_mapt.ndc_frontal_smRNAseq.txt",
                                   targets.path = "../mirna_target_analysis_0719/MAPT_DEG_targets.txt",
                                   mrna = mrna,
                                   srna = srna,
                                   pval = pval,
                                   lfc = lfc,
                                   cor_cutoff = cor_cutoff)

# GRN
grn_targets <- generateTargetList(deg.srna.path = "../analysis_0719/deseq_result_grn.ndc_frontal_smRNAseq.txt",
                                   targets.path = "../mirna_target_analysis_0719/GRN_DEG_targets.txt",
                                   mrna = mrna,
                                   srna = srna,
                                   pval = pval,
                                   lfc = lfc,
                                   cor_cutoff = cor_cutoff)


# C9ORF72
c9_targets <- generateTargetList(deg.srna.path = "../analysis_0719/deseq_result_c9.ndc_frontal_smRNAseq.txt",
                                   targets.path = "../mirna_target_analysis_0719/C9_DEG_targets.txt",
                                   mrna = mrna,
                                   srna = srna,
                                   pval = pval,
                                   lfc = lfc,
                                   cor_cutoff = cor_cutoff)

#=====================================================================#

#####
# Save output
######


### MAPT
deg <- read.table("../analysis_0719/deseq_result_mapt.ndc_frontal_smRNAseq.txt", sep="\t", header=T)
deg <- deg[deg$padj <= pval,]
deg.up <- deg[deg$log2FoldChange > lfc,]
deg.down <- deg[deg$log2FoldChange < -lfc,]
deg.up <- mapt_targets[names(mapt_targets) %in% as.character(deg.up$X)]
deg.down <- mapt_targets[names(mapt_targets) %in% as.character(deg.down$X)]
mapt.up <- as.character(unlist(deg.up))
mapt.down <- as.character(unlist(deg.down))
mapt.up <- mapt.up[!duplicated(mapt.up)]
mapt.down <- mapt.down[!duplicated(mapt.down)]

write.table(mapt.up, "MAPT_upMir_correlated_targets_Refseq.txt", quote=F, row.names=F)
write.table(mapt.down, "MAPT_downMir_correlated_targets_Refseq.txt", quote=F, row.names=F)

### GRN
deg <- read.table("../analysis_0719/deseq_result_grn.ndc_frontal_smRNAseq.txt", sep="\t", header=T)
deg <- deg[deg$padj <= pval,]
deg.up <- deg[deg$log2FoldChange > lfc,]
deg.down <- deg[deg$log2FoldChange < -lfc,]
deg.up <- grn_targets[names(grn_targets) %in% as.character(deg.up$X)]
deg.down <- grn_targets[names(grn_targets) %in% as.character(deg.down$X)]
grn.up <- as.character(unlist(deg.up))
grn.down <- as.character(unlist(deg.down))
grn.up <- grn.up[!duplicated(grn.up)]
grn.down <- grn.down[!duplicated(grn.down)]

write.table(grn.up, "GRN_upMir_correlated_targets_Refseq.txt", quote=F, row.names=F)
write.table(grn.down, "GRN_downMir_correlated_targets_Refseq.txt", quote=F, row.names=F)


### C9orf72
deg <- read.table("../analysis_0719/deseq_result_c9.ndc_frontal_smRNAseq.txt", sep="\t", header=T)
deg <- deg[deg$padj <= pval,]
deg.up <- deg[deg$log2FoldChange > lfc,]
deg.down <- deg[deg$log2FoldChange < -lfc,]
deg.up <- c9_targets[names(c9_targets) %in% as.character(deg.up$X)]
deg.down <- c9_targets[names(c9_targets) %in% as.character(deg.down$X)]
c9.up <- as.character(unlist(deg.up))
c9.down <- as.character(unlist(deg.down))
c9.up <- c9.up[!duplicated(c9.up)]
c9.down <- c9.down[!duplicated(c9.down)]


write.table(c9.up, "C9_upMir_correlated_targets_Refseq.txt", quote=F, row.names=F)
write.table(c9.down, "C9_downMir_correlated_targets_Refseq.txt", quote=F, row.names=F)


###
# Generate miRNA-target (edge) tables
###

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# mapt
mirs <- c()
targets <- c()
count <- 1
for (n in names(mapt_targets)){
  print(count)
  count <- count + 1
  # get targets
  tgts <- as.character(unlist(mapt_targets[n]))
  # get gene symbols for targets
  bm <- getBM(attributes = c("hgnc_symbol", "refseq_mrna"), filters="refseq_mrna", values=tgts, mart=ensembl)
  tgts <- bm$hgnc_symbol
  tgts <- tgts[!duplicated(tgts)]
  # add to vectors
  mirs <- c(mirs, rep(n, length(tgts)))
  targets <- c(targets, tgts)
}
mapt.df <- data.frame(mirna = mirs, targets = targets)

# grn
mirs <- c()
targets <- c()
count <- 1
for (n in names(grn_targets)){
  print(count)
  count <- count + 1
  # get targets
  tgts <- as.character(unlist(grn_targets[n]))
  # get gene symbols for targets
  bm <- getBM(attributes = c("hgnc_symbol", "refseq_mrna"), filters="refseq_mrna", values=tgts, mart=ensembl)
  tgts <- bm$hgnc_symbol
  tgts <- tgts[!duplicated(tgts)]
  # add to vectors
  mirs <- c(mirs, rep(n, length(tgts)))
  targets <- c(targets, tgts)
}
grn.df <- data.frame(mirna = mirs, targets = targets)

# c9orf72
mirs <- c()
targets <- c()
count <- 1
for (n in names(c9_targets)){
  print(count)
  count <- count + 1
  # get targets
  tgts <- as.character(unlist(c9_targets[n]))
  # get gene symbols for targets
  bm <- getBM(attributes = c("hgnc_symbol", "refseq_mrna"), filters="refseq_mrna", values=tgts, mart=ensembl)
  tgts <- bm$hgnc_symbol
  tgts <- tgts[!duplicated(tgts)]
  # add to vectors
  mirs <- c(mirs, rep(n, length(tgts)))
  targets <- c(targets, tgts)
}
c9.df <- data.frame(mirna = mirs, targets = targets)

# save tables
write.table(mapt.df, "MAPT_miRNA_target_edge_table.txt", sep="\t", row.names = F, quote=F)
write.table(grn.df, "GRN_miRNA_target_edge_table.txt", sep="\t", row.names = F, quote=F)
write.table(c9.df, "C9_miRNA_target_edge_table.txt", sep="\t", row.names = F, quote=F)
