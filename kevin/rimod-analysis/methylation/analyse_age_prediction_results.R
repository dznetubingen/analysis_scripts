#####
# Analysis of age prediction results
#####
library(ggplot2)
library(stringr)
setwd("~/rimod/Methylation/age_prediction/")

beta <- read.table("rimod_frontal_pheno_data.txt")

pred <- read.table("predictions.csv")
beta$prediction <- pred$X0

md <- read.csv("~/rimod/files/FTD_Brain.csv")
md <- md[md$REGION == "temporal",]
md$sample <- str_pad(md$SAMPLEID, width=5, side="left", pad="0")

beta$sample <- str_sub(beta$Sample_Name, 1, 5)

# merging
md <- md[md$sample %in% beta$sample,]
beta <- beta[beta$sample %in% md$sample,]

md <- md[match(beta$sample, md$sample),]

beta$age <- md$AGE
beta$dc <- md$DISEASE.CODE


mad <- function(x, y){
  return(mean(abs(x - y)))
}
rmse <- function(x, y){
  return(sqrt(mean((x - y)^2)))
}

# control
cont <- beta[beta$dc == "control",]
print(mad(cont$prediction, cont$age))
print(rmse(cont$prediction, cont$age))
print(mean(cont$age - cont$prediction))

# mapt
ftd <- beta[beta$dc == "FTD-MAPT",]
print(mad(ftd$prediction, ftd$age))
print(rmse(ftd$prediction, ftd$age))
print(mean(ftd$age - ftd$prediction))

# mapt
ftd <- beta[beta$dc == "FTD-GRN",]
print(mad(ftd$prediction, ftd$age))
print(rmse(ftd$prediction, ftd$age))
print(mean(ftd$age - ftd$prediction))

# C9orf72
ftd <- beta[beta$dc == "FTD-C9",]
print(mad(ftd$prediction, ftd$age))
print(rmse(ftd$prediction, ftd$age))
print(mean(ftd$age - ftd$prediction))

ftd <- beta[beta$dc %in% c("FTD-C9", "FTD-MAPT", "FTD-GRN"),]
print(mad(ftd$prediction, ftd$age))
print(rmse(ftd$prediction, ftd$age))

print(mad(beta$prediction, beta$age))
print(rmse(beta$prediction, beta$age))

# Do some plotting
df1 <- data.frame(sample = beta$sample, group = beta$dc, age = beta$prediction, category = rep("Prediction", nrow(beta)))
df2 <- data.frame(sample = beta$sample, group = beta$dc, age = beta$age, category = rep("Truth", nrow(beta)))
df <- rbind(df1, df2)


p <- ggplot(df, aes(x=sample, y=age, group=category)) + 
  geom_line(aes(color=category)) +
  facet_wrap(~group)
p

p <- ggplot(df, aes(x=group, y=age, fill=category)) + 
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45))
p
ggsave("~/rimod/paper/figures/figure4/age_prediction_boxplot.png", width=3, height = 3)
