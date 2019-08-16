library(CAGEr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
setwd("/data")


# Check if CAGE TSS already saved  - otherwise load from scratch
if (!file.exists("cage_object_frontal.rds")){
  files <- list.files("all_ctss/", full.names = T)
  files <- files[grepl("_tem_", files)]
  cage <- new("CAGEset",
              genomeName = "BSgenome.Hsapiens.UCSC.hg38",
              inputFiles = files,
              inputFilesType = "ctss",
              sampleLabels = gsub(".ctss", "", basename(files)))
  
  # load data
  getCTSS(cage)
  
  saveRDS(cage, "cage_object_tem.rds")
  gc()
}
cage <- readRDS("cage_object_tem.rds")

# Plot Reverse Cumultativ curves for all samples
#plotReverseCumulatives(cage, fitInRange = c(5, 1000), onePlot = TRUE)

# Perform normlization (TPM)
normalizeTagCount(cage, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.2, T = 10^6)

# clustering
clusterCTSS(cage, threshold = 1, thresholdIsTpm = T,
            nrPassThreshold = 1, method = "distclu", maxDist = 20,
            removeSingletons = TRUE, keepSingletonsAbove = 5)

# Further downstream analyses
cumulativeCTSSdistribution(cage, clusters = "tagClusters")

quantilePositions(cage, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)

plotInterquantileWidth(cage, clusters = "tagClusters", tpmThreshold = 3, qLow = 0.1, qUp = 0.9)

aggregateTagClusters(cage, tpmThreshold = 1, qLow = 0.1, qUp = 0.9, maxDist = 100)

# Score Promotor shift
cumulativeCTSSdistribution(cage, clusters = "consensusClusters")
saveRDS(cage, "aggregated_clusters_TPM_tem.Rds")

# Define groups of samples
library(stringr)
md <- read.csv("FTD_Brain.csv")
md$samples <- str_pad(as.character(md$SAMPLEID), width=5, pad="0", side='left')
md$samples <- paste0("sample_", md$samples)
md$samples <- paste0(md$samples, "_tem_no_chrM")
md <- md[md$samples %in% as.character(sampleLabels(cage)),]
md <- md[!duplicated(md$samples),]
mapt <- md[md$DISEASE.CODE == "FTD-MAPT",]$samples
grn <- md[md$DISEASE.CODE == 'FTD-GRN',]$samples
c9 <- md[md$DISEASE.CODE == 'FTD-C9',]$samples
ndc <- md[md$DISEASE.CODE == 'control',]$samples

# MAPT shifting
mapt.shift <- cage
scoreShift(mapt.shift, groupX = mapt, groupY = ndc)
mapt.shifting.promotors <- getShiftingPromoters(mapt.shift)
write.table(mapt.shifting.promotors, "mapt_shifting_promotors_tem.txt", sep="\t", quote=F, row.names=F)

# GRN shifting
grn.shift <- cage
scoreShift(grn.shift, groupX = grn, groupY = ndc)
grn.shifting.promotors <- getShiftingPromoters(grn.shift)
write.table(grn.shifting.promotors, "grn_shifting_promotors_tem.txt", sep="\t", quote=F, row.names=F)


# C9orf72 shifting
c9.shift <- cage
scoreShift(c9.shift, groupX = c9, groupY = ndc)
c9.shifting.promotors <- getShiftingPromoters(c9.shift)
write.table(c9.shifting.promotors, "c9_shifting_promotors_tem.txt", sep="\t", quote=F, row.names=F)






