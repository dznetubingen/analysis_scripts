####################################################
# Analysis of RiMod ChIP-seq data using ChIPseeker
####################################################
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/chipseq/")

# GRN
files = list.files("results_grn/bwa/mergedLibrary/macs/broadPeak/", pattern="*.broadPeak", full.names = T)

peak = readPeakFile(files[1])
covplot(peak, weightCol="X41")

# annotation
peakAnno <- annotatePeak(files[3], TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene )
plotAnnoPie(peakAnno)
