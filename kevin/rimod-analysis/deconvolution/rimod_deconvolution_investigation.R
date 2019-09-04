###############################
# Deconvolution results investigation
# Script to check if the results of deconvolution for RiMod make sens
################################
library(stringr)
library(ggplot2)
library(reshape2)

setwd("~/rimod/RNAseq/analysis/deconvolution/")

# load fractions
fracs <- read.table("cdn_predictions.txt", sep="\t", header=T, row.names=1)
rownames(fracs) <- gsub("X", "", rownames(fracs))
samples <- str_split(rownames(fracs), pattern = "_", simplify = T)[,1]
samples <- str_pad(samples, width=5, side='left', pad="0")

# load MD file
md <- read.csv("~/rimod/files/FTD_Brain.csv")
md <- md[md$REGION == "frontal",]
md$id <- str_pad(md$SAMPLEID, width=5, side='left', pad="0")
md <- md[md$id %in% samples,]
df <- data.frame(sample = md$id, dc = md$DISEASE.CODE)

# remove sample we don't have in MD
keep <- samples %in% df$sample
samples <- samples[keep]
fracs <- fracs[keep,]

# match
df <- df[match(samples, df$sample),]

fracs$dc <- df$dc

# merge neurons
fracs$Neurons <- fracs$InNeurons + fracs$ExNeurons
fracs <- fracs[fracs$dc %in% c("control", "FTD-C9", "FTD-MAPT", "FTD-GRN"),]

# Plotting
df <- melt(fracs, id.vars = "dc")

# neurons
neuro <- df[df$variable == "Neurons",]

p <- ggplot(neuro, aes(x=dc, y=value, color=dc)) +
  geom_boxplot()
p
