##########
# miRNA boxplots
##########

setwd("~/rimod/paper_v2/figures/figure7/")

# load expression
mat <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_vst_values_frontal_smRNA.txt", sep="\t", header=T, row.names = 1, check.names = F)

# load metadata
md <- read.table("~/rimod/smallRNA/frontal/rimod_human_frontal_smRNAseq_metadata.txt", sep="\t", header=T)
md <- md[match(colnames(mat), md$sample),]

# color palette only for disease groups
mypal <- c("#616665", "#7570B3", "#db6e1a","#67e08a")

# hsa-mir-150-5p
tmp <- as.numeric(mat["hsa-miR-150-5p",])
df <- data.frame(Expression = tmp, Group = md$dc)

p <- ggplot(df, aes(x=Group, y=Expression, color=Group)) + 
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_manual(values = mypal) + 
  ggtitle("hsa-miR-150-5p")
p
ggsave("hsa-miR-150-5p_plot.png", width=3.2, height=2)


# hsa-miR-193a-3p
tmp <- as.numeric(mat['hsa-miR-193a-3p',])
df <- data.frame(Expression = tmp, Group = md$dc)

p <- ggplot(df, aes(x=Group, y=Expression, color=Group)) + 
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + theme(legend.position = "none") +
  scale_fill_manual(values = mypal) + 
  ggtitle("hsa-miR-193a-3p")
p
ggsave("hsa-miR-193a-3p_plot.png", width=2, height=2)


# hsa-miR-19b-3p
tmp <- as.numeric(mat['hsa-miR-19b-3p',])
df <- data.frame(Expression = tmp, Group = md$dc)

p <- ggplot(df, aes(x=Group, y=Expression, color=Group)) + 
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + theme(legend.position = "none") +
  scale_fill_manual(values = mypal) + 
  ggtitle("hsa-miR-19b-3p")
p
ggsave("hsa-miR-19b-3p_plot.png", width=2, height=2)

