#####################
# Make boxplots of highest DE genes
#####################
library(ggplot2)
library(viridis)
library(stringr)

setwd("~/rimod/paper_v2/figures/supplements/glutamatergic_synapse/")


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


# color palette only for disease groups
mypal <- c("#616665", "#7570B3", "#db6e1a","#67e08a")

# Prepare colors

# color palette only for disease groups
mypal <- c("#616665", "#7570B3", "#db6e1a","#67e08a")

# GLUA1
glua1 = "ENSG00000155511"
glua1 <- as.numeric(mat[glua1,])
df <- data.frame(Expression = glua1, Group = md$Disease.Code)

p <- ggplot(df, aes(x=Group, y=Expression, fill=Group)) + 
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_manual(values = mypal) + theme(legend.position = "none") +
  ggtitle("GLUA1") 
p
ggsave("GLUA1_plot.png", width=2, height=2)

# GLUA2
glua2 = "ENSG00000120251"
glua2 <- as.numeric(mat[glua2,])
df <- data.frame(Expression = glua2, Group = md$Disease.Code)

p <- ggplot(df, aes(x=Group, y=Expression, fill=Group)) + 
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_manual(values = mypal) + theme(legend.position = "none") +
  ggtitle("GLUA2") 
p
ggsave("GLUA2_plot.png", width=2, height=2)


# GLUA3
glua3 = "ENSG00000125675"
glua3 <- as.numeric(mat[glua3,])
df <- data.frame(Expression = glua3, Group = md$Disease.Code)

p <- ggplot(df, aes(x=Group, y=Expression, fill=Group)) + 
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_manual(values = mypal) + theme(legend.position = "none") +
  ggtitle("GLUA3") 
p
ggsave("GLUA3_plot.png", width=2, height=2)

# GLUA4
glua4 = "ENSG00000152578"
glua4 <- as.numeric(mat[glua4,])
df <- data.frame(Expression = glua4, Group = md$Disease.Code)

p <- ggplot(df, aes(x=Group, y=Expression, fill=Group)) + 
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_manual(values = mypal) + theme(legend.position = "none") +
  ggtitle("GLUA4") 
p
ggsave("GLUA4_plot.png", width=2, height=2)


# GRIP1
grip = "ENSG00000155974"
grip <- as.numeric(mat[grip,])
df <- data.frame(Expression = grip, Group = md$Disease.Code)

p <- ggplot(df, aes(x=Group, y=Expression, fill=Group)) + 
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_manual(values = mypal) + theme(legend.position = "none") +
  ggtitle("GRIP1") 
p
ggsave("GRIP1_plot.png", width=2, height=2)
