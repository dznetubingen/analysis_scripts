library(stringr)

setwd("/Users/kevin/dzne/")

gtf <- read.table("gencode.mirna.gtf", sep="\t", stringsAsFactors = F)

mir <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/analysis_0719/deseq_result_mapt.ndc_frontal_smRNAseq.txt", sep="\t", header=T, row.names=1, stringsAsFactors = F)
mir <- mir[mir$padj <= 0.05,]

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
ats <- listAttributes(ensembl)[,1]
ats[grepl("mir", ats)]

bm <- getBM(attributes = c("ensembl_gene_id", "mirbase_id", "hgnc_symbol"), filters = "mirbase_id", values=rownames(mir), mart=ensembl)

# manually modify miRNA names
mirs <- rownames(mir)
mirs <- gsub("-5p", "", gsub("-3p", "", mirs))
mirs <- gsub("hsa", "", gsub("-", "", mirs))
mirs <- gsub("p", "", mirs)
mirs <- toupper(mirs)

mir$names <- mirs

ids <- str_split(gtf$V9, pattern=";", simplify = T)
genes <- c()
for (i in 1:nrow(ids)){
  for (j in 1:ncol(ids)){
    if (grepl("gene_name", ids[i,j])){
      genes <- c(genes, str_split(ids[i,j], pattern=" ", simplify=T)[,3])
    }
  }
}
gtf$gene <- genes

mir <- mir[mir$names %in% gtf$gene,]
gtf <- gtf[gtf$gene %in% mir$names,]
gtf <- gtf[gtf$V3 == "gene",]

# load methylation data
met <- read.table("/Users/kevin/dzne/rimod_package/frontal_methylation_0818/DMPs_mapt.ndc_quant.txt", sep="\t", header=T, stringsAsFactors = F)
met <- met[met$adj.P.Val <= 0.05,]

# for every miRNA, look for matching CpGs
cutoff <- 2000
hits <- c()
hit_mirs <- c()
for (i in 1:nrow(gtf)) {
  m <- gtf$gene[i]
  chr <- gtf$V1[i]
  start = gtf$V4[i]
  end = gtf$V5[i]
  strand = gtf$V7[i]
  print(m)

  for (j in 1:nrow(met)){
    if (chr == met$chr[j]){
      if (abs(met$pos[j] - start) < cutoff){
        print("hit")
        hits <- c(hits, j)
        hit_mirs <- c(hit_mirs, m)
      }
      else if (abs(met$pos[j] - end) < cutoff){
        print("hit")
        hits <- c(hits, j)
        hit_mirs <- c(hit_mirs, m)
      }
    }
  }
  
}
hits <- met[hits,]
hits$mir <- hit_mirs


met <- read.table("/Users/kevin/dzne/rimod_package/frontal_methylation_0818/DMPs_c9orf72.ndc_quant.txt", sep="\t", header=T, stringsAsFactors = F)
met <- met[met$adj.P.Val <= 0.05,]
dim(met)

