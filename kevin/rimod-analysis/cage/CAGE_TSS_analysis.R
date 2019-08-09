############################################################
# CTSS focused analysis of CAGE-seq data
############################################################
library(CAGEr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/CAGEseq/")

# Load data
inputFiles = list.files(path = "all_ctss/", pattern = "ctss$", full.names=TRUE)

# Only perform analysis of frontal samples
inputFiles <- inputFiles[grepl("_fro_", inputFiles)]


cage <- new("CAGEset",
            genomeName = "BSgenome.Hsapiens.UCSC.hg38",
            inputFiles = inputFiles,
            inputFilesType = "ctss",
            sampleLabels = gsub(".ctss", "", basename(inputFiles)))

# load data
getCTSS(cage)

plotReverseCumulatives(cage, fitInRange = c(5, 1000), onePlot = TRUE)

# normalization
normalizeTagCount(cage, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.2, T = 10^6)

# clustering
clusterCTSS(cage, threshold = 1, thresholdIsTpm = T,
            nrPassThreshold = 1, method = "distclu", maxDist = 20,
            removeSingletons = TRUE, keepSingletonsAbove = 5)

# stuff
cumulativeCTSSdistribution(cage, clusters = "tagClusters")

quantilePositions(cage, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)

plotInterquantileWidth(cage, clusters = "tagClusters", tpmThreshold = 3, qLow = 0.1, qUp = 0.9)

aggregateTagClusters(cage, tpmThreshold = 1, qlow = 0.1, qUp = 0.9, maxDist = 100)






