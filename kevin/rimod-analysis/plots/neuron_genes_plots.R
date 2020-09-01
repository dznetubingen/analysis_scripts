################
# Make publication ready figures of iPSC-neuronal gene expression
#################
library(stringr)
library(reshape2)
library(ggplot2)
library(biomaRt)

setwd("~/rimod/Neurons/rimod_neurons/")

# Load data
mat <- read.table("RiMod_iPSCNeurons_vst_values.txt", sep="\t", header=T, row.names=1, check.names = F)
mat <- read.table("frontal_lengthScaeldTPM_counts.txt", sep="\t", header=T, row.names=1, check.names = F)
rownames(mat) <- str_split(rownames(mat), pattern = "[.]", simplify = T)[,1]
md <- read.table("iPSCNeurons_mRNAseq_metadata.txt", sep="\t", header=T, stringsAsFactors = F)
md <- md[md$sample %in% colnames(mat),]
md <- md[match(colnames(mat), md$sample),]


# Make HGNC row.names
ensembl <- useEnsembl("ensembl", mirror = 'asia', dataset = "hsapiens_gene_ensembl")
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=rownames(mat), mart=ensembl)
mat <- merge(mat, bm, by.x="row.names", by.y="ensembl_gene_id")
mat <- mat[!duplicated(mat$hgnc_symbol),]
rownames(mat) <- mat$hgnc_symbol
mat <- mat[, -1]

# genes to plot
genes <- c("GRN", "C9orf72", "MAPT", "TIMP2", "SOX2", "MMP2", "PLGLB2", "CHCHD10", "MMP15")

df <- melt(mat[genes,])
df <- merge(df, md, by.x="variable", by.y="sample")
df$group <- factor(df$group, levels = c("control", "GRN", "MAPT", "C9"))
mypal <- c("#616665", "#db6e1a","#67e08a", "#7570B3")

ggplot(df, aes(x=group, y=value, fill=group)) + 
  geom_boxplot() +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  scale_fill_manual(values = mypal) +
  facet_wrap(~ hgnc_symbol, scales="free_y", ncol=3)

ggsave("~/rimod/paper/figures/figure6/boxplot_neurons.png", width=6, height=4)


# next lot
# genes to plot
genes <- c("TBL1X", "GPSM1", "WNT9A", "SOX2", "RFX2", "ADGRB1")
df <- melt(mat[genes,])
df <- merge(df, md, by.x="variable", by.y="sample")
df$group <- factor(df$group, levels = c("control", "GRN", "MAPT", "C9"))
ggplot(df, aes(x=group, y=value, fill=group)) + 
  geom_boxplot() +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  scale_fill_manual(values = mypal) +
  facet_wrap(~ hgnc_symbol, scales="free_y")