library(pathview)
library(biomaRt)

setwd("~/rimod/integrative_analysis/immune_system_pathway_analysis/")
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

data(paths.hsa)

for (i in 264:length(paths.hsa)) {
  print(i)
  pw <- gsub("hsa", "", names(paths.hsa)[i])
  pname <- paths.hsa[i]
  pname <- gsub(" ", "", gsub("/", "", pname))
  
  
  
  # C9orf72
  deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2020-05-04_15.45.57/deseq_result_c9.ndc_fro_2020-05-04_15.45.57.txt", sep="\t", header=T)
  #deg <- deg[deg$padj <= 0.05,]
  bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=deg$X, mart=ensembl)
  deg <- merge(deg, bm, by.x="X", by.y="ensembl_gene_id")
  deg <- deg[!deg$hgnc_symbol == "",]
  
  gene.data <- deg$log2FoldChange
  names(gene.data) <- deg$hgnc_symbol
  
  v <- pathview(gene.data = gene.data,
                pathway.id = pw,
                gene.idtype = "SYMBOL",
                out.suffix = paste("C9orf72", pname, sep=""),
                node.sum = 'mean')
  
  # GRN
  deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2020-05-04_15.45.57/deseq_result_grn.ndc_fro_2020-05-04_15.45.57.txt", sep="\t", header=T)
  #deg <- deg[deg$padj <= 0.05,]
  bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=deg$X, mart=ensembl)
  deg <- merge(deg, bm, by.x="X", by.y="ensembl_gene_id")
  deg <- deg[!deg$hgnc_symbol == "",]
  
  gene.data <- deg$log2FoldChange
  names(gene.data) <- deg$hgnc_symbol
  
  v <- pathview(gene.data = gene.data,
                pathway.id = pw,
                gene.idtype = "SYMBOL",
                out.suffix = paste("GRN", pname, sep=""),
                node.sum = 'mean')
  
  # MAPT
  deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2020-05-04_15.45.57/deseq_result_mapt.ndc_fro_2020-05-04_15.45.57.txt", sep="\t", header=T)
  #deg <- deg[deg$padj <= 0.05,]
  bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=deg$X, mart=ensembl)
  deg <- merge(deg, bm, by.x="X", by.y="ensembl_gene_id")
  deg <- deg[!deg$hgnc_symbol == "",]
  
  gene.data <- deg$log2FoldChange
  names(gene.data) <- deg$hgnc_symbol
  
  v <- pathview(gene.data = gene.data,
                pathway.id = pw,
                gene.idtype = "SYMBOL",
                out.suffix = paste("MAPT", pname, sep=""),
                node.sum = 'mean')
  
}



####
# Individual pathways for inspection
####
print(i)
pw <- '04724'
pname <- "GlutamatergicSynapse"


# C9orf72
deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2020-05-04_15.45.57/deseq_result_c9.ndc_fro_2020-05-04_15.45.57.txt", sep="\t", header=T)
#deg <- deg[deg$padj <= 0.05,]
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=deg$X, mart=ensembl)
deg <- merge(deg, bm, by.x="X", by.y="ensembl_gene_id")
deg <- deg[!deg$hgnc_symbol == "",]

gene.data <- deg$log2FoldChange
names(gene.data) <- deg$hgnc_symbol


v <- pathview(gene.data = gene.data,
              pathway.id = pw,
              out.suffix = paste("C9orf72", pname, sep=""),
              gene.idtype = "SYMBOL",
              node.sum = 'mean')

# GRN
deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2020-05-04_15.45.57/deseq_result_grn.ndc_fro_2020-05-04_15.45.57.txt", sep="\t", header=T)
#deg <- deg[deg$padj <= 0.05,]
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=deg$X, mart=ensembl)
deg <- merge(deg, bm, by.x="X", by.y="ensembl_gene_id")
deg <- deg[!deg$hgnc_symbol == "",]

gene.data <- deg$log2FoldChange
names(gene.data) <- deg$hgnc_symbol

v <- pathview(gene.data = gene.data,
              pathway.id = pw,
              gene.idtype = "SYMBOL",
              out.suffix = paste("GRN", pname, sep=""),
              node.sum = 'mean')

# MAPT
deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2020-05-04_15.45.57/deseq_result_mapt.ndc_fro_2020-05-04_15.45.57.txt", sep="\t", header=T)
#deg <- deg[deg$padj <= 0.05,]
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=deg$X, mart=ensembl)
deg <- merge(deg, bm, by.x="X", by.y="ensembl_gene_id")
deg <- deg[!deg$hgnc_symbol == "",]

gene.data <- deg$log2FoldChange
names(gene.data) <- deg$hgnc_symbol

v <- pathview(gene.data = gene.data,
              pathway.id = pw,
              gene.idtype = "SYMBOL",
              out.suffix = paste("MAPT", pname, sep=""),
              node.sum = 'mean')

