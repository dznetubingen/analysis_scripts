#######
# Methylation Data pre-processing for age prediction
#######
setwd("~/rimod/Methylation/age_prediction/")
library(RnBeads)
library(RnBeads.hg38)
library(impute)
# Data paths
data.dir <- "/home/kevin/rimod/Methylation/age_prediction/rnbeads_analysis/idata_data/"
sample.sheet.path <- "/home/kevin/rimod/Methylation/age_prediction/rnbeads_analysis/8632_Samplesheet_2018-130-ILL_MET_rnbeads.csv"
report.dir = "/home/kevin/rimod/Methylation/age_prediction/rnbeads_analysis/rnbeads_results/"


###
# Specificy RNB Options
rnb.options(
  assembly = "hg38",
  filtering.sex.chromosomes.removal = TRUE,
  normalization.method = "bmiq",
  normalization.background.method = "methylumi.noob"
)

## Import the data
dataSource <- c(data.dir, sample.sheet.path)
result <- rnb.run.import(dataSource, data.type="infinium.idat.dir", 
                         dir.reports=report.dir)
rnbs <- result$rnb.set

##
# Perform preprocessing and normalization
rnb.options(
  filtering.sex.chromosomes.removal = TRUE,
  normalization.method = "bmiq",
  normalization.background.method = "methylumi.noob",
  filtering.snp = "3"
)
result <- rnb.run.preprocessing(rnbs, report.dir)
rnbs <- result$rnb.set

## Save data
mm <- meth(rnbs, row.names = T)
ph <- pheno(rnbs)

# Subset with Wenda-CpGs
wenda <- read.table("CpG_sites_reduced.txt", sep="\t")
wenda_cpgs <- as.character(wenda$V1)
common_cps <- intersect(wenda_cpgs, rownames(mm))
mm <- mm[common_cps,]

## Perform imputation
imp.res <- impute.knn(mm)
mm.imputed <- imp.res$data

# make dummy age for Wenda
ph$age <- rep(50, nrow(ph))

write.table(mm.imputed, "rimod_frontal_methylation_data.txt")
write.table(ph, "rimod_frontal_pheno_data.txt")



####
# reduce the wenda datasets to common CpGs
# test data
test <- read.csv("wenda_data/3-TestReducedMeth.csv", sep=" ", header=T)
test <- test[common_cps,]
write.table(test, "wenda_data/test_reduced_epic.csv")
# train data
train <- read.csv("wenda_data/3-TrainingReducedMeth.csv", sep=" ", header=T)
train <- train[common_cps,]
write.table(train, "wenda_data/train_reduced_epic.csv")

# testing
test <- read.csv("wenda_data/3-TestReducedPheno.csv", sep=" ")
