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
df <- data.frame(sample = md$id, dc = md$DISEASE.CODE, age=md$AGE, gender=md$GENDER, gene=md$GENE)

# remove sample we don't have in MD
keep <- samples %in% df$sample
samples <- samples[keep]
fracs <- fracs[keep,]

# match
df <- df[match(samples, df$sample),]

fracs$dc <- df$dc
fracs$age <- df$age
fracs$gender <- df$gender
fracs$gene <- df$gene

# merge neurons
fracs$Neurons <- fracs$InNeurons + fracs$ExNeurons
fracs <- fracs[fracs$dc %in% c("control", "FTD-C9", "FTD-MAPT", "FTD-GRN"),]

# Plotting
df <- melt(fracs, id.vars = c("dc", "age", "gender", "gene"))

# neurons
neuro <- df[df$variable == "Neurons",]

p <- ggplot(neuro, aes(x=dc, y=value, color=dc)) +
  geom_violin() +
  geom_point()
ggsave("~/tmp_stuff/pics/rimod_dc.png")
p

p <- ggplot(neuro, aes(x=dc, y=value, color=age, alpha=0.7)) +
  geom_point(size=10)
ggsave("~/tmp_stuff/pics/rimod_age.png")
p

p <- ggplot(neuro, aes(x=dc, y=value, color=gender, alpha=0.7)) +
  geom_point(size=5)
ggsave("~/tmp_stuff/pics/rimod_gender.png")
p

p <- ggplot(neuro, aes(x=dc, y=value, color=gene, alpha=0.7)) +
  geom_point(size=5)
ggsave("~/tmp_stuff/pics/rimod_gene.png")
p




