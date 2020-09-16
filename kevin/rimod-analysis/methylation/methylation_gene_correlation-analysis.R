#############
# Perform correlation of genes with CpGs for
# all different groups
# Correlatione expressin of OPTN with
# methylation of cg01532982
library(stringr)
library(biomaRt)

#####
# FTD-MAPT
#####
# load rna and methylation data
rna <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2020-05-04_15.45.57/deseq_vst_values_2020-05-04_15.45.57.txt", sep="\t",  header=T, row.names = 1)
rownames(rna) <- str_split(rownames(rna), pattern="[.]", simplify = T)[,1]
beta <- read.table("~/rimod/Methylation/frontal_methylation_0818/betaVals_matrix_frontal_metyhlation.txt", sep="\t", header=T)
beta <- read.table("~/rimod/Methylation/frontal_methylation_0818/mVals_matrix_frontal_methylation.txt", sep="\t", header=T)
met <- beta
# get CpGs
grn <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_grn.ndc_quant.txt", sep="\t", header=T)
mapt <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_mapt.ndc_quant.txt", sep="\t", header=T)
c9 <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_c9orf72.ndc_quant.txt", sep="\t", header=T)
grn <- grn[grn$adj.P.Val <= 0.05,]
grn <- grn[grn$logFC > 0.6,]
mapt <- mapt[mapt$adj.P.Val <= 0.05,]
mapt <- mapt[mapt$logFC > 0.6,]
c9 <- c9[c9$adj.P.Val <= 0.05,]
c9 <- c9[c9$logFC > 0.6,]



grn <- grn[, c(-41, -40, -39, -38, -37, -36)]
mapt <- mapt[, c(-41, -40, -39, -38, -37, -36)]
c9 <- c9[, c(-41, -40, -39, -38, -37, -36)]
cpg <- rbind(mapt, grn, c9)
cpg <- cpg[!duplicated(rownames(cpg)),]

# Transform row names to HGNC symbol (of RNA-seq data)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=rownames(rna), mart=ensembl)
rna <- merge(rna, bm, by.x="row.names", by.y="ensembl_gene_id")
rna <- rna[!duplicated(rna$hgnc_symbol),]
rownames(rna) <- rna$hgnc_symbol
rna <- rna[, c(-1, -ncol(rna))]


# format colnames
colnames(rna) <- gsub("X", "", colnames(rna))
colnames(rna) <- str_split(colnames(rna), pattern="_", simplify = T)[,1]
colnames(met) <- gsub("X", "", colnames(met))
colnames(met) <- str_split(colnames(met), pattern="_", simplify = T)[,1]
cmn <- intersect(colnames(met), colnames(rna))
rna <- rna[, colnames(rna) %in% cmn]
met <- met[, colnames(met) %in% cmn]
rna <- rna[, match(colnames(met), colnames(rna))]


# Calc correlation for every gene
gene_cor <- c()
cpg_list <- c()

for (i in 1:nrow(rna)) {
  if (i %% 100 == 0){
    print(i)
  }
  gene <- rownames(rna)[i]
  exp <- as.numeric(rna[i,])
  
  tmp <- cpg[grepl(gene, cpg$GencodeBasicV12_NAME),]
  if (nrow(tmp) > 0){
    for (j in 1:nrow(tmp)) {
      tmp.cpg <- rownames(tmp)[j]
      exp.cpg <- as.numeric(met[tmp.cpg,])
      
      res <- cor(exp, exp.cpg)
      gene_cor <- c(gene_cor, res)
      cpg_list <- c(cpg_list, tmp.cpg)
    }
  }
}

# remove duplicats
keep <- !duplicated(cpg_list)
cpg_list <- cpg_list[keep]
gene_cor <- gene_cor[keep]

keep <- !is.na(gene_cor)
cpg_list <- cpg_list[keep]
gene_cor <- gene_cor[keep]
