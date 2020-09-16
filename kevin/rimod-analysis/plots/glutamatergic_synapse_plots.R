#####################
# Make boxplot for CHCHD10 (any maybe some other genes)
#####################
library(ggplot2)
library(viridis)
library(stringr)

setwd("~/rimod/paper_v2/figures/figure3/")


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

# EAAT5
EAAT5 = "ENSG00000162383"
test <- as.numeric(mat[EAAT5,])
df <- data.frame(Expression = test, Group = md$Disease.Code)

p <- ggplot(df, aes(x=Group, y=Expression, color=Group)) + 
  geom_boxplot() + geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_manual(values = mypal) + 
  ggtitle("EAAT5")
p
ggsave("EAAT5_plot.png", width=3, height=2)


# GRIA1
GRIA1 = "ENSG00000155511"
test <- as.numeric(mat[GRIA1,])
df <- data.frame(Expression = test, Group = md$Disease.Code)

p <- ggplot(df, aes(x=Group, y=Expression, color=Group)) + 
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + theme(legend.position = "none") +
  scale_fill_manual(values = mypal) + 
  ggtitle("GRIA1")
p
ggsave("GRIA1_plot.png", width=1.8, height=2)

# GRIA2
GRIA2 = "ENSG00000120251"
test <- as.numeric(mat[GRIA2,])
df <- data.frame(Expression = test, Group = md$Disease.Code)

p <- ggplot(df, aes(x=Group, y=Expression, color=Group)) + 
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + theme(legend.position = "none") +
  scale_fill_manual(values = mypal) + 
  ggtitle("GRIA2")
p
ggsave("GRIA2_plot.png", width=1.8, height=2)


# GRIA3
GRIA3 = "ENSG00000125675"
test <- as.numeric(mat[GRIA3,])
df <- data.frame(Expression = test, Group = md$Disease.Code)

p <- ggplot(df, aes(x=Group, y=Expression, color=Group)) + 
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + theme(legend.position = "none") +
  scale_fill_manual(values = mypal) + 
  ggtitle("GRIA3")
p
ggsave("GRIA3_plot.png", width=1.8, height=2)


# GRIA4
GRIA4 = "ENSG00000152578"
test <- as.numeric(mat[GRIA4,])
df <- data.frame(Expression = test, Group = md$Disease.Code)

p <- ggplot(df, aes(x=Group, y=Expression, color=Group)) + 
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + theme(legend.position = "none") +
  scale_fill_manual(values = mypal) + 
  ggtitle("GRIA4")
p
ggsave("GRIA4_plot.png", width=1.8, height=2)



test <- "ENSG00000060982"
test <- as.numeric(mat[test,])
df <- data.frame(Expression = test, Group = md$Disease.Code)

p <- ggplot(df, aes(x=Group, y=Expression, color=Group)) + 
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + theme(legend.position = "none") +
  scale_fill_manual(values = mypal) + 
  ggtitle("test")
p