#####################
# Make boxplot for CHCHD10 (any maybe some other genes)
#####################
library(ggplot2)
library(viridis)
library(stringr)

setwd("~/rimod/paper_v2/figures/figure5//")


# color palette only for disease groups
mypal <- c("#616665", "#7570B3", "#db6e1a","#67e08a")

# load expression matrix
mat <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2020-05-04_15.45.57/deseq_vst_values_2020-05-04_15.45.57.txt", sep="\t", header=T, row.names = 1)
rownames(mat) <- str_split(rownames(mat), pattern="[.]", simplify = T)[,1]
colnames(mat) <- gsub("X", "", colnames(mat))
colnames(mat) <- str_pad(str_split(colnames(mat), pattern="_", simplify = T)[,1], width=5, side="left", pad="0")

# load md
md <- read.table("~/rimod/RNAseq/rimod_frontal_rnaseq_metadata.txt", sep="\t", header=T)
md <- md[match(colnames(mat), md$SampleID),]

# filter out sporadic
keep <- !md$Disease.Code == "Sporadic-TDP"
md <- md[keep,]
mat <- mat[,keep]

# load DEG results
mapt <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/MAPT_cell_composition_filterd_DEGs.txt", sep="\t", header=T)
grn <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/GRN_cell_composition_filterd_DEGs.txt", sep="\t", header=T)
c9 <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/C9orf72_cell_composition_filterd_DEGs.txt", sep="\t", header=T)

# order by fold changes
mapt <- mapt[order(abs(mapt$log2FoldChange), decreasing = T),]
grn <- grn[order(abs(grn$log2FoldChange), decreasing = T),]
c9 <- c9[order(abs(c9$log2FoldChange), decreasing = T),]

# CHCHD10
CHCHD10 = "ENSG00000250479"
test <- as.numeric(mat[CHCHD10,])
df <- data.frame(Expression = test, Group = md$Disease.Code)

p <- ggplot(df, aes(x=Group, y=Expression, color=Group)) + 
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_manual(values = mypal) + 
  ggtitle("CHCHD10")
p
ggsave("CHCHD10_plot.png", width=3, height=2)


