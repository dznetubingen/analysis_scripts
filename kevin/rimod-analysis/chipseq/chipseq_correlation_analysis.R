################
# correlation ChIP-seq - CAGE-seq
################
library(stringr)
library(edgeR)
setwd("~/rimod/integrative_analysis/")

chip <- read.table("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/chipseq/brain_cemil/RiMod_ChIPseq_annotated_count_table.txt", sep="\t", header=T)
cage <- read.table("~/rimod/CAGE/results_annotation/RiMod_aggrGeneCounts_CAGEseq_fro.txt", sep="\t", header=T)

# formatting
colnames(cage) <- gsub("_fro", "", colnames(cage))
rownames(chip) <- str_split(rownames(chip), pattern="[.]", simplify = T)[,1]
# subset to common samples
cmn <- intersect(colnames(chip), colnames(cage))
chip <- chip[cmn]
cage <- cage[cmn]

# subset to common genes
genes <- intersect(rownames(cage), rownames(chip))

chip <- chip[genes,]
cage <- cage[genes,]


# library size normalization
chip <- cpm(chip)
cage <- cpm(cage)

# Sample-wise correlation
sample_cor <- c()
for (i in 1:ncol(cage)) {
  res <- cor(cage[,i], chip[,i])
  sample_cor <- c(sample_cor, res)
}
print(mean(sample_cor))


# Gene-wise correlation
gene_cor <- c()
for (i in 1:nrow(cage)) {
  res <- cor(as.numeric(cage[i,]), as.numeric(chip[i,]))
  gene_cor <- c(gene_cor, res)
}
print(mean(na.omit(gene_cor)))

