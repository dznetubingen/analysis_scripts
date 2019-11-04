setwd("~/rimod/")
library(stringr)

rna <- read.table("RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_vst_values_2019-10-23_13.33.11.txt", sep="\t", row.names = 1, header = T)
rownames(rna) <- str_split(rownames(rna), patter="[.]", simplify = T)[,1]

cage <- read.table("CAGE/cage_analysis/CAGE_deseq_analysis_2019-08-15_14.49.08_frontal/deseq_rLog_values_2019-08-15_14.49.08.txt",
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
# cor_list <- c()
# for (i in 1:nrow(rna)){
#   print(i)
#   cage.tmp <- as.numeric(cage[i,])
#   rna.tmp <- as.numeric(rna[i,])
#   tmp <- cor(cage.tmp, rna.tmp)
#   cor_list <- c(cor_list, tmp)
# }
# print(mean(cor_list))


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

outname = "~/rimod/paper/figures/figure2/rna_cage_correlation.png"
png(outname, height=300, width=300)
p
dev.off()




