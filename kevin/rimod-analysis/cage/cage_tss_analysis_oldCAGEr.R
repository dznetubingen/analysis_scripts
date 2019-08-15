library(CAGEr)
library(BSgenome.Hsapiens.UCSC.hg38)
setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/CAGEseq/")


# Check if CAGE TSS already saved  - otherwise load from scratch
if (!file.exists("cage_object_frontal.rds")){
  files <- list.files("all_ctss/", full.names = T)
  files <- files[grepl("_fro_", files)]
  cage <- new("CAGEset",
              genomeName = "BSgenome.Hsapiens.UCSC.hg38",
              inputFiles = files,
              inputFilesType = "ctss",
              sampleLabels = gsub(".ctss", "", basename(files)))
  
  # load data
  getCTSS(cage)
  
  saveRDS(cage, "cage_object_frontal.rds")
  gc()
}
cage <- readRDS("cage_object_frontal.rds")

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

aggregateTagClusters(cage, tpmThreshold = 1, qlow = 0.1, qUp = 0.9, maxDist = 100)