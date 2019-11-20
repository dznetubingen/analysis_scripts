######################
# CAGE-seq analysis with CAGEfightR
#########################
library(CAGEfightR)
library(rtracklayer)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
setwd("/data/kevin")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Create file lists
bw_plus_files = list.files("bigwig/", pattern = "*.plus.bw", full.names = T)
bw_minus_files = list.files("bigwig/", pattern = "*.minus.bw", full.names = T)
bw_plus = BigWigFileList(bw_plus_files)
bw_minus = BigWigFileList(bw_minus_files)
names(bw_plus) <- gsub(".plus.bw", "", basename(bw_plus_files))
names(bw_minus) <- gsub(".minus.bw", "", basename(bw_minus_files))
samples <- gsub(".plus.bw", "", basename(bw_plus_files))

# To be on safe side, make sure both lists have the same ordering
order(samples, names(bw_plus))
order(samples, names(bw_minus))

# Genome information 
gi = rtracklayer::SeqinfoForUCSCGenome("hg38")

# Load Data
ctss <- quantifyCTSSs(plusStrand = bw_plus,
                      minusStrand = bw_minus,
                      genome = gi)

###
# CTSS level analysis
###

# Calc TPM
ctss <- calcTPM(ctss, inputAssay = "counts", outputAssay = "TPM")

# Calc pooled CTSS
ctss <- calcPooled(ctss, inputAssay = "TPM")

# Calc CTSS support, subset data and recalculate TPM (uncommented currently)
ctss <- calcSupport(ctss, inputAssay = "counts")
#ctss <- subset(ctss, rowRanges(ctss)$support > 1)
#ctss <- calcTPM(ctss)
#ctss <- calcPooled(ctss)

# Define TSS clusters
cage_TCs <- clusterUnidirectionally(ctss, pooledCutoff=0.1, mergeDist=20) 

# Find Enhancer -> bidirectional clustering
BCs <- clusterBidirectionally(ctss, balanceThreshold=0.95)

BCs <- calcBidirectionality(BCs, samples=ctss)
enhancers <- subset(BCs, BCs$bidirectionality > 0)

##
# Cluster level analysis
##
quant_tss <- quantifyClusters(ctss, clusters = cage_TCs, inputAssay = "counts")
quant_enhancers <- quantifyClusters(ctss, clusters = enhancers, inputAssay = "counts")

# subset clusters by support
quant_tss <- subsetBySupport(quant_tss, inputAssay = "counts", unexpressed = 0, minSamples = 1)
quant_enhancers <- subsetBySupport(quant_enhancers, inputAssay = "counts", unexpressed = 0, minSamples = 1)

###
# Annotate with Gene expression
###
quant_tss <- assignTxID(quant_tss, txModels = txdb, outputColumn="txID")
quant_enhancers <- assignTxID(quant_enhancers, txModels = txdb, outputColumn="txID")

quant_tss <- assignTxType(quant_tss, txModels = txdb, outputColumn="txTyp")
quant_enhancers <- assignTxType(quant_enhancers, txModels = txdb, outputColumn="txTyp")

# Only keep intronic and intergenic enhancers
txType = rowRanges(quant_enhancers)$txTyp
quant_enhancers <- subset(quant_enhancers, txType %in% c("intron", "intergenic"))


###
# Find TSS-enhancer pairs
###
# mark the cluster type
rowRanges(quant_tss)$clusterType <- "TSS"
rowRanges(quant_enhancers)$clusterType = "enhancer"

# Merge ranges
colData(quant_enhancers) <- colData(quant_tss)
SE <- combineClusters(object1 = quant_tss,
                      object2 = quant_enhancers,
                      removeIfOverlapping="object1")

# Mark TSS as reference points
rowRanges(SE)$clusterType <- factor(rowRanges(SE)$clusterType, levels = c("TSS", "enhancer"))

# Recalculate TPM values
SE <- calcTPM(SE, totalTags = "totalTags")

# Find links between the groups
TCBC_pairs <- findLinks(SE,
                        inputAssay="TPM",
                        maxDist=10000,
                        directional="clusterType",
                        method="kendall")

en_tss_pairs <- data.frame(TCBC_pairs)

# subset for significant associations
en_tss_pairs <- en_tss_pairs[en_tss_pairs$p.value <= 0.01,]

# Only keep pairs with strong support
en_tss_pairs <- en_tss_pairs[en_tss_pairs$support1 > 5, ]
write.table(en_tss_pairs, "Frontal_enhancer_tss_pairs.txt", sep="\t", quote=F)
