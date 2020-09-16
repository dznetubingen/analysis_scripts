setwd("~/rimod/")
library(stringr)

rna <- read.table("RNAseq/analysis/RNAseq_analysis_fro_2020-05-04_15.45.57/deseq_vst_values_2020-05-04_15.45.57.txt", sep="\t", row.names = 1, header = T)
rownames(rna) <- str_split(rownames(rna), patter="[.]", simplify = T)[,1]

cage <- read.table("CAGE/cage_analysis/CAGE_deseq_analysis_2020-05-05_17.35.50_frontal/deseq_vst_values_2020-05-05_17.35.50.txt",
                   sep="\t", row.names=1, header=T)

genes <- intersect(rownames(rna), rownames(cage))

rna <- rna[genes,]
cage <- cage[genes,]

# colnames
rna.samples <- colnames(rna)
rna.samples <- gsub("X", "", rna.samples)
rna.samples <- str_sub(rna.samples, 1, 5)
colnames(rna) <- rna.samples

cage.samples <- colnames(cage)
cage.samples <- gsub("sample_", "", cage.samples)
cage.samples <- str_sub(cage.samples, 1, 5)
colnames(cage) <- cage.samples

samples <- intersect(rna.samples, cage.samples)

rna <- rna[, samples]
cage <- cage[, samples]

# # calculate correlation for every gene
cor_list <- c()
for (i in 1:nrow(rna)){
   print(i)
   cage.tmp <- as.numeric(cage[i,])
   rna.tmp <- as.numeric(rna[i,])
   tmp <- cor(cage.tmp, rna.tmp)
   cor_list <- c(cor_list, tmp)
}
print(mean(cor_list))


sample_cor <- c()
for (j in 1:ncol(rna)){
  cage.tmp <- as.numeric(cage[,j])
  rna.tmp <- as.numeric(rna[,j])
  tmp <- cor(cage.tmp, rna.tmp)
  sample_cor <- c(sample_cor, tmp)
}
print(mean(sample_cor))

library(ggplot2)
library(reshape2)

df = melt(cage)
df2 = melt(rna)
df$RNAseq <- df2$value
colnames(df) <- c("Sample", "CAGEseq", "RNAseq")


# density plot
p <- ggplot(df, aes(x=RNAseq, y=CAGEseq)) + 
  geom_point() + 
  geom_density_2d() + 
  stat_density_2d(aes(fill = ..level..), geom="polygon") +
  scale_fill_gradient(low="blue", high="red") +
  theme_minimal(base_size = 15)
p

outname = "~/rimod/paper_v2/figures/supplements/rna_cage_correlation.png"
png(outname, height=300, width=300)
p
dev.off()

df <- data.frame(Correlation = sample_cor, Sample = samples)
p <- ggplot(df, aes(x=Sample, y=Correlation)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
p

outname = "~/rimod/paper_v2/figures/supplements/barplot_sample_correlation_RNAseq_CAGEseq.png"
png(outname, height=300, width=300)
p
dev.off()

########################
# Fold Change testing of RNA-seq and CAGE-seq DE analysis results
########################

# MAPT
rna <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_result_mapt.ndc_fro_2019-10-23_13.33.11.txt", sep="\t", header=T, row.names = 1)
cage <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_fro_2019-10-23_13.30.53/deseq_result_mapt.ndc_2019-10-23_13.30.53.txt", sep="\t", header=T, row.names = 1)
common.genes <- intersect(rownames(rna), rownames(cage))
rna <- rna[common.genes,]
cage <- cage[common.genes,]

# check percentage of genes with common fold change direction
hits <- c()
nohits <- c()
for (i in 1:nrow(rna)) {
  rna.fc <- rna[i,2]
  cage.fc <- cage[i,2]
  if (rna.fc > 0 && cage.fc > 0) {
    hits <- c(hits, i)
  }
  else if (rna.fc < 0 && cage.fc < 0){
    hits <- c(hits, i)
  }
  else  {
    nohits <- c(nohits, i)
  }
}
hit.pct <- length(hits) / nrow(rna)
print(hit.pct)

# average fold changes
# RNA
tmp <- rna[hits,]
print(mean(abs(tmp$log2FoldChange)))
tmp <- rna[nohits,]
print(mean(abs(tmp$log2FoldChange)))

# CAGE
tmp <- cage[hits,]
print(mean(abs(tmp$log2FoldChange)))
tmp <- cage[nohits,]
print(mean(abs(tmp$log2FoldChange)))

# Save common (cell composition filtered) DEGs
rna.up <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/MAPT_filterd_DEGs_up.txt", sep="\t")
cage.up <- read.table("~/rimod/CAGE/cage_analysis/deconvolution_frontal/cell_specificity_filtering/MAPT_filterd_DEGs_up.txt", sep="\t")
up <- intersect(rna.up$V1, cage.up$V1)
rna.down <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/MAPT_filterd_DEGs_down.txt", sep="\t")
cage.down <- read.table("~/rimod/CAGE/cage_analysis/deconvolution_frontal/cell_specificity_filtering/MAPT_filterd_DEGs_down.txt", sep="\t")
down <- intersect(rna.down$V1, cage.down$V1)

write.table(up, "~/rimod/integrative_analysis/rnaseq_cageseq_degs/MAPT_up_ccf_intersect.txt", row.names = F, col.names = F, quote=F)
write.table(down, "~/rimod/integrative_analysis/rnaseq_cageseq_degs/MAPT_down_ccf_intersect.txt", row.names = F, col.names = F, quote=F)


# GRN
rna <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_result_grn.ndc_fro_2019-10-23_13.33.11.txt", sep="\t", header=T, row.names = 1)
cage <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_fro_2019-10-23_13.30.53/deseq_result_grn.ndc_2019-10-23_13.30.53.txt", sep="\t", header=T, row.names = 1)
common.genes <- intersect(rownames(rna), rownames(cage))
rna <- rna[common.genes,]
cage <- cage[common.genes,]

# check percentage of genes with common fold change direction
hits <- c()
nohits <- c()
for (i in 1:nrow(rna)) {
  rna.fc <- rna[i,2]
  cage.fc <- cage[i,2]
  if (rna.fc > 0 && cage.fc > 0) {
    hits <- c(hits, i)
  }
  else if (rna.fc < 0 && cage.fc < 0){
    hits <- c(hits, i)
  }
  else  {
    nohits <- c(nohits, i)
  }
}
hit.pct <- length(hits) / nrow(rna)
print(hit.pct)

# average fold changes
# RNA
tmp <- rna[hits,]
print(mean(abs(tmp$log2FoldChange)))
tmp <- rna[nohits,]
print(mean(abs(tmp$log2FoldChange)))

# CAGE
tmp <- cage[hits,]
print(mean(abs(tmp$log2FoldChange)))
tmp <- cage[nohits,]
print(mean(abs(tmp$log2FoldChange)))

# Save common (cell composition filtered) DEGs
rna.up <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/GRN_filterd_DEGs_up.txt", sep="\t")
cage.up <- read.table("~/rimod/CAGE/cage_analysis/deconvolution_frontal/cell_specificity_filtering/GRN_filterd_DEGs_up.txt", sep="\t")
up <- intersect(rna.up$V1, cage.up$V1)
rna.down <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/GRN_filterd_DEGs_down.txt", sep="\t")
cage.down <- read.table("~/rimod/CAGE/cage_analysis/deconvolution_frontal/cell_specificity_filtering/GRN_filterd_DEGs_down.txt", sep="\t")
down <- intersect(rna.down$V1, cage.down$V1)

write.table(up, "~/rimod/integrative_analysis/rnaseq_cageseq_degs/GRN_up_ccf_intersect.txt", row.names = F, col.names = F, quote=F)
write.table(down, "~/rimod/integrative_analysis/rnaseq_cageseq_degs/GRN_down_ccf_intersect.txt", row.names = F, col.names = F, quote=F)



# C9orf72
rna <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_result_c9.ndc_fro_2019-10-23_13.33.11.txt", sep="\t", header=T, row.names = 1)
cage <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_fro_2019-10-23_13.30.53/deseq_result_c9.ndc_2019-10-23_13.30.53.txt", sep="\t", header=T, row.names = 1)
common.genes <- intersect(rownames(rna), rownames(cage))
rna <- rna[common.genes,]
cage <- cage[common.genes,]

# check percentage of genes with common fold change direction
hits <- c()
nohits <- c()
for (i in 1:nrow(rna)) {
  rna.fc <- rna[i,2]
  cage.fc <- cage[i,2]
  if (rna.fc > 0 && cage.fc > 0) {
    hits <- c(hits, i)
  }
  else if (rna.fc < 0 && cage.fc < 0){
    hits <- c(hits, i)
  }
  else  {
    nohits <- c(nohits, i)
  }
}
hit.pct <- length(hits) / nrow(rna)
print(hit.pct)

# average fold changes
# RNA
tmp <- rna[hits,]
print(mean(abs(tmp$log2FoldChange)))
tmp <- rna[nohits,]
print(mean(abs(tmp$log2FoldChange)))

# CAGE
tmp <- cage[hits,]
print(mean(abs(tmp$log2FoldChange)))
tmp <- cage[nohits,]
print(mean(abs(tmp$log2FoldChange)))

# Save common (cell composition filtered) DEGs
rna.up <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/C9orf72_filterd_DEGs_up.txt", sep="\t")
cage.up <- read.table("~/rimod/CAGE/cage_analysis/deconvolution_frontal/cell_specificity_filtering/C9orf72_filterd_DEGs_up.txt", sep="\t")
up <- intersect(rna.up$V1, cage.up$V1)
rna.down <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/C9orf72_filterd_DEGs_down.txt", sep="\t")
cage.down <- read.table("~/rimod/CAGE/cage_analysis/deconvolution_frontal/cell_specificity_filtering/C9orf72_filterd_DEGs_down.txt", sep="\t")
down <- intersect(rna.down$V1, cage.down$V1)

write.table(up, "~/rimod/integrative_analysis/rnaseq_cageseq_degs/C9orf72_up_ccf_intersect.txt", row.names = F, col.names = F, quote=F)
write.table(down, "~/rimod/integrative_analysis/rnaseq_cageseq_degs/C9orf72_down_ccf_intersect.txt", row.names = F, col.names = F, quote=F)
