#############
# Correlatione expressin of OPTN with
# methylation of cg01532982
library(stringr)
library(biomaRt)

# get genes
genes <- read.table("~/rimod/paper/figures/figure5/GRN_M3down_hyperMethylatin.txt", sep="\t")
genes <- genes$V1
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="hgnc_symbol", values=genes, mart=ensembl)

# get CpGs
cpg <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_grn.ndc_quant.txt", sep="\t", header=T)


# load RNA-seq
exp <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2020-05-04_15.45.57/deseq_vst_values_2020-05-04_15.45.57.txt", sep="\t", header=T)
exp$X <- str_split(exp$X, pattern="[.]", simplify = T)[,1]
rownames(exp) <- exp$X


# load methylation
#met <- read.table("~/rimod/Methylation/frontal_methylation_0818/mVals_matrix_frontal_methylation.txt", sep="\t", header=T)
beta <- read.table("~/rimod/Methylation/frontal_methylation_0818/betaVals_matrix_frontal_metyhlation.txt", sep="\t", header=T)
met <- beta


# format colnames
colnames(exp) <- gsub("X", "", colnames(exp))
colnames(exp) <- str_split(colnames(exp), pattern="_", simplify = T)[,1]
colnames(met) <- gsub("X", "", colnames(met))
colnames(met) <- str_split(colnames(met), pattern="_", simplify = T)[,1]
cmn <- intersect(colnames(met), colnames(exp))
exp <- exp[, colnames(exp) %in% cmn]
met <- met[, colnames(met) %in% cmn]
exp <- exp[, match(colnames(met), colnames(exp))]


sig_genes <- c()
for (i in 1:nrow(bm)) {
  hg <- bm$hgnc_symbol[i]
  ens <- bm$ensembl_gene_id[i]
  
  if (ens %in% rownames(exp)){
    e <- exp[ens,]
    rownames(e) <- e$X
    e <- as.numeric(e)
    
    tmp <- cpg[grepl(hg, cpg$GencodeBasicV12_NAME),]
    for (cp in rownames(tmp)) {
      #print(cp)
      m <- met[cp,]
      m <- as.numeric(m)
      
      res <- cor.test(m, e)
      pval <- res$p.value
      res.cor <- res$estimate
      
      if (pval < 0.05){
        print(hg)
        print(cp)
        print(res.cor)
        print(pval)
        sig_genes <- c(sig_genes, hg)
      }
    }
  }
  

}







