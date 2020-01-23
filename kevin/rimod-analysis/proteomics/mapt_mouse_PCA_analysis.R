setwd("~/rimod/Proteomics/mouse/")

mat <- read.csv("P301L_mice_LFC.csv", sep=";", stringsAsFactors = F)
mat <- as.data.frame(mat)
mat <- mat[, -1]

for (i in 1:nrow(mat)){
  for (j in 1:ncol(mat)){
    mat[i,j] <- gsub(",", ".", mat[i,j])
    if (mat[i,j] == "NaN"){
      mat[i,j] <- 0
    }
    mat[i,j] <- as.numeric(mat[i,j])
  }
}

m <- mapply(mat, FUN=as.numeric)
m <- matrix(data=m, ncol=ncol(mat), nrow=nrow(mat))


pca_df <- prcomp(t(m))

library(ggplot2)
library(stringr)

df <- data.frame(PC1=pca_df$x[,1], PC2=pca_df$x[,2], sample=colnames(mat))
group <- rep("control", nrow(df))
group[grepl("Tg", df$sample)] <- "TG"
df$group <- group

age <- str_split(df$sample, pattern="_", simplify = T)[,2]
df$age <- age

ggplot(df, aes(x=PC1, y=PC2, color=group, shape=age)) +
  geom_point(size=10)

