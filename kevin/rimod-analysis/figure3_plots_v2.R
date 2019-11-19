#####
# Visualizatin of HumanBase modules
####
setwd("~/rimod/RNAseq/analysis/human_base/")


df <- read.table("rnaseq_grn_filtered_down_enrichment.txt", sep="\t", header=T)
mods <- as.character(levels(df$CLUSTER_NAME))

# module 1
m1 <- mods[mods$CLUSTER_NAME == "M4",]
m1 <- m1[order(m1$TERM_Q_VALUE),]
