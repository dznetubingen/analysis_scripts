###############
# Methylation module analysis
###############
library(stringr)
library(ggplot2)
library(reshape)
library(extrafont)
loadfonts()
setwd("~/rimod/paper/figures/figure4/")


####
# MAPT
####


# Load modules
mod.down <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_mapt_filtered_down_modules.txt", header=T)
mod.up <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_mapt_filtered_up_modules.txt", header=T)
#mod.down <- as.character(str_split(mod.down$CLUSTER_GENES, pattern=",", simplify = T))


# load methylation
met.up <- read.table("~/rimod/Methylation/frontal_methylation_0818/MAPT_UP_CpGs_genes.txt", stringsAsFactors = F)$V1
met.down <- read.table("~/rimod/Methylation/frontal_methylation_0818/MAPT_DOWN_CpGs_genes.txt", stringsAsFactors = F)$V1
met <- c(met.up, met.down)

# data frame
df <- data.frame("Module" = "dummy", "DMPs" = 0, "Size" = 0)

# Up-regulated modules
for (m in mod.up$CLUSTER_NAME) {
  tmp <- mod.up[mod.up$CLUSTER_NAME == m,]
  tmp <- as.character(str_split(tmp$CLUSTER_GENES, pattern=",", simplify = T))
  size = length(tmp)
  print(size)
  dmps = length(intersect(tmp, met))
  print(dmps)
  tmp <- data.frame("Module" = paste(m, "up", sep="-"), "DMPs" = dmps, "Size" = size)
  df = rbind(df, tmp)
}


# Down-regulated modules
for (m in mod.down$CLUSTER_NAME) {
  tmp <- mod.down[mod.down$CLUSTER_NAME == m,]
  tmp <- as.character(str_split(tmp$CLUSTER_GENES, pattern=",", simplify = T))
  size = length(tmp)
  dmps = length(intersect(tmp, met))
  tmp <- data.frame("Module" = paste(m, "down", sep="-"), "DMPs" = dmps, "Size" = size)
  df = rbind(df, tmp)
}

df <- df[-1,]
df <- melt(df)
colnames(df) <- c("Module", "DMPs", "Count")
dmps <- as.character(df$DMPs)
dmps[dmps == "DMPs"] <- "Yes"
dmps[dmps == "Size"] <- "No"
df$DMPs = dmps


# Plotting section
p <- ggplot(df, aes(x=Module, y=Count, fill=DMPs)) +
  geom_bar(stat="identity", position="fill") +
  theme_minimal() + 
  coord_flip()
p

ggsave("MAPT_methylation_module_plot.png", width=3, height=3)
#======= END MAPT ==============#



####
# GRN
####


# Load modules
mod.down <- read.table("~/rimod/RNAseq//analysis/human_base/rnaseq_grn_filtered_down_modules.txt", header=T)
mod.down <- mod.down[!mod.down$CLUSTER_NAME == "",]
mod.up <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_grn_filtered_up_modules.txt", header=T)


# load methylation
met.up <- read.table("~/rimod/Methylation/frontal_methylation_0818/GRN_DMPs_UP.txt", stringsAsFactors = F)$V1
met.down <- read.table("~/rimod/Methylation/frontal_methylation_0818/GRN_DMPs_DOWN.txt", stringsAsFactors = F)$V1
met <- c(met.up, met.down)

# data frame
df <- data.frame("Module" = "dummy", "DMPs" = 0, "Size" = 0)

# Up-regulated modules
for (m in mod.up$CLUSTER_NAME) {
  tmp <- mod.up[mod.up$CLUSTER_NAME == m,]
  tmp <- as.character(str_split(tmp$CLUSTER_GENES, pattern=",", simplify = T))
  size = length(tmp)
  dmps = length(intersect(tmp, met))
  tmp <- data.frame("Module" = paste(m, "up", sep="-"), "DMPs" = dmps, "Size" = size)
  df = rbind(df, tmp)
}


# Down-regulated modules
for (m in mod.down$CLUSTER_NAME) {
  tmp <- mod.down[mod.down$CLUSTER_NAME == m,]
  tmp <- as.character(str_split(tmp$CLUSTER_GENES, pattern=",", simplify = T))
  size = length(tmp)
  dmps = length(intersect(tmp, met))
  tmp <- data.frame("Module" = paste(m, "down", sep="-"), "DMPs" = dmps, "Size" = size)
  df = rbind(df, tmp)
}

df <- df[-1,]
df <- melt(df)
colnames(df) <- c("Module", "DMPs", "Count")
dmps <- as.character(df$DMPs)
dmps[dmps == "DMPs"] <- "Yes"
dmps[dmps == "Size"] <- "No"
df$DMPs = dmps


# Plotting section
p <- ggplot(df, aes(x=Module, y=Count, fill=DMPs)) +
  geom_bar(stat="identity", position="fill") +
  theme_minimal() + 
  coord_flip()
p

ggsave("GRN_methylation_module_plot.png", width=3, height=3)

#======= END GRN ==============#