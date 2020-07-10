#########
# RiMod table
#########
library(stringr)
setwd("~/rimod/")

md <- read.csv("files/FTD_Brain_corrected.csv", stringsAsFactors = F)
md <- md[md$REGION == "frontal",]
md <- md[!grepl("spor", md$DISEASE.CODE),]

# match with used cage data
cage <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2020-05-05_17.35.50_frontal/deseq_vst_values_2020-05-05_17.35.50.txt", sep="\t", header=T, row.names = 1)
samples <- colnames(cage)
samples <- gsub("sample_", "", gsub("_fro", "", samples))
md$sample <- str_pad(as.character(md$SAMPLEID), width=5, pad="0", side="left")
md <- md[md$sample %in% samples,]

# Generate metrics for all datasets
g <- "FTD-MAPT"
df <- md[md$DISEASE.CODE == g,]
age.mean <- mean(df$AGE)
age.sd <- sd(df$AGE)

pmi <- na.omit(as.numeric(df$PMD.MIN.))
pmi.mean <- mean(pmi)
pmi.sd <- sd(pmi)

ph <- na.omit(as.numeric(df$PH))
ph.mean <- mean(ph)
ph.sd <- sd(ph)

rin <- na.omit(as.numeric(df$RIN))
rin.mean <- mean(rin)
rin.sd <- sd(rin)

gender <- data.frame(table(df$GENDER))
print(gender)
pct.male <- gender[gender$Var1 == "M",]$Freq / nrow(df)
metrics.df <- data.frame(MeanAge = age.mean, SDAge = age.sd, MeanPMI = pmi.mean, SDPMI = pmi.sd, MeanPH = ph.mean, SDPH = ph.sd,
                  MeanRIN = rin.mean, SDRIN = rin.sd, PctMale = pct.male, Group = g)

groups <- c("FTD-GRN", "FTD-C9", "control")
for (g in groups) {
  df <- md[md$DISEASE.CODE == g,]
  
  df <- md[md$DISEASE.CODE == g,]
  age.mean <- mean(df$AGE)
  age.sd <- sd(df$AGE)
  pmi <- na.omit(as.numeric(df$PMD.MIN.))
  pmi.mean <- mean(pmi)
  pmi.sd <- sd(pmi)
  
  ph <- na.omit(as.numeric(df$PH))
  ph.mean <- mean(ph)
  ph.sd <- sd(ph)
  
  rin <- na.omit(as.numeric(df$RIN))
  rin.mean <- mean(rin)
  rin.sd <- sd(rin)
  
  gender <- data.frame(table(df$GENDER))
  print(gender)
  pct.male <- gender[gender$Var1 == "M",]$Freq / nrow(df)
  tmp <- data.frame(MeanAge = age.mean, SDAge = age.sd, MeanPMI = pmi.mean, SDPMI = pmi.sd, MeanPH = ph.mean, SDPH = ph.sd,
                    MeanRIN = rin.mean, SDRIN = rin.sd, PctMale = pct.male, Group = g)
  
  metrics.df <- rbind(metrics.df, tmp)
}
