###############################
# Test of trackViewer library #
###############################
library(Gviz)
library(rtracklayer)
library(trackViewer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

# MAPT test
mapt_entrez <- get("MAPT", org.Hs.egSYMBOL2EG)
theTrack <- geneTrack(mapt_entrez,TxDb.Hsapiens.UCSC.hg38.knownGene)[[1]]
