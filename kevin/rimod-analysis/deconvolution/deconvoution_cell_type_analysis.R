##################
# Test for enrichment of certain cell types on RNA-seq data
##################
library(stringr)

setwd("~/rimod/RNAseq/analysis/deconvolution/")


####
# Cell composition plot
####

fracs <- read.table("~/rimod/RNAseq/analysis/deconvolution/cdn_predictions.txt", sep="\t", header=T)
colnames(fracs)[1] <- "sample"
fracs$sample <- gsub("X", "", fracs$sample)
fracs$sample <- str_split(fracs$sample, pattern="_", simplify = T)[,1]
fracs$Neurons <- fracs$InNeurons + fracs$ExNeurons
fracs <- fracs[, c(-3, -9)]



# get design matrix
md <- read.csv("~/rimod/files/FTD_Brain.csv")
md <- md[md$REGION == "frontal",]
md$sample <- str_split(md$GIVENSAMPLENAME, pattern="_", simplify = T)[,1]
md <- md[md$sample %in% fracs$sample,]
md <- md[match(fracs$sample, md$sample),]
fracs$group <- as.character(md$DISEASE.CODE)

fracs$group[as.character(md$GENE) == "P301L"] <- "MAPT-P301L"

# remove sporadic
fracs <- fracs[!fracs$group == "Sporadic-TDP",]

# remove Unknown
fracs <- fracs[,-1]

####
# Calculate average increase/decrase
####
mean_fun <- median
fracs <- fracs[, -1]
mapt <- fracs[fracs$group == "FTD-MAPT",]
grn <- fracs[fracs$group == "FTD-GRN",]
c9 <- fracs[fracs$group == "FTD-C9",]
mp3 <- fracs[fracs$group == "MAPT-P301L",]

control <- fracs[fracs$group == "control",]
mapt <- mapt[, -ncol(mapt)]
mapt <- apply(mapt, 2, mean_fun)
mp3 <- mp3[, -ncol(mp3)]
mp3 <- apply(mp3, 2, mean_fun)
grn <- grn[, -ncol(grn)]
grn <- apply(grn, 2, mean_fun)
control <- control[, -ncol(control)]
control <- apply(control, 2, mean_fun)
c9 <- c9[, -ncol(c9)]
c9 <- apply(c9, 2, mean_fun)




# MAPT
mapt_list <- c()
for (i in 1:length(control)) {
  ct.mean <- control[i]
  diff <- mapt[i] - ct.mean
  pct <- diff / ct.mean
  print(pct)
  mapt_list <- c(mapt_list, pct)
}

# MAPT-P301L
mp3_list <- c()
for (i in 1:length(control)) {
  ct.mean <- control[i]
  diff <- mp3[i] - ct.mean
  pct <- diff / ct.mean
  print(pct)
  mp3_list <- c(mp3_list, pct)
}


# GRN
grn_list <- c()
for (i in 1:length(control)) {
  ct.mean <- control[i]
  diff <- grn[i] - ct.mean
  pct <- diff / ct.mean
  print(pct)
  grn_list <- c(grn_list, pct)
}

# C9ORF72
c9_list <- c()
for (i in 1:length(control)) {
  ct.mean <- control[i]
  diff <- c9[i] - ct.mean
  pct <- diff / ct.mean
  print(pct)
  c9_list <- c(c9_list, pct)
}

# make data frame
grn <- t(data.frame(grn_list))
mapt <- t(data.frame(mapt_list))
mp3 <- t(data.frame(mp3_list))
c9 <- t(data.frame(c9_list))
df <- data.frame(rbind(mapt, mp3, grn, c9))
df$group <- c("MAPT", "MAPT-P301L", "GRN", "C9ORF72")


# plotting
library(reshape2)
library(ggplot2)

df <- melt(df)

ggplot(data=df, aes(x=variable, y=value, fill=group)) +
  geom_bar(stat='identity', position = position_dodge()) +
  theme_minimal()

