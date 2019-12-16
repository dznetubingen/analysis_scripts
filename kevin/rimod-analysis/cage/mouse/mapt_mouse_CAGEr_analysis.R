######
# MAPT mouse CAGE-seq analysis (CAGEr)
#####
library(CAGEr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(stringr)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/mouse/mapt_mice/")

#ctss_files <- list.files("results/ctss", full.names = T)
#sample_labels <- str_split(basename(ctss_files), pattern="[.]", simplify = T)[,1]

bam_files <- list.files("results/STAR", full.names = T, pattern="*.bam")
sample_labels <- str_split(basename(bam_files), pattern="[.]", simplify = T)[,1]


cageset <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10", inputFiles = bam_files, 
               inputFilesType = "bam", sampleLabels = sample_labels)

getCTSS(cageset, removeFirstG = F, correctSystematicG = F)
ctss <- CTSStagCount(cageset)

save(cageset, file="cageset_object.RData")
#load("cageset_object.RData")
#ctss <- CTSStagCount(cageset)

# Remove sites only available in few samples
tmp <- ctss[, c(-1, -2, -3)]
rs <- apply(tmp, 1, sum)
keep <- rs > 1

ctss <- ctss[keep,]
cageset <- cageset[keep,]

# plot reverse cumulative
plotReverseCumulatives(cageset, fitInRange = c(5, 1000), onePlot = T)

# normalize (don't actually normalize)
normalizeTagCount(cageset, method="none")

# cluster CTSS
clusterCTSS(cageset, threshold = 1, thresholdIsTpm = F, nrPassThreshold = 10, method = "distclu", maxDist = 20, removeSingletons = T,
            keepSingletonsAbove = 5)
tc <- tagClusters(cageset, sample = sample_labels[1])

# aggregate clusters
aggregateTagClusters(cageset)

# make matrix of all consensus clusters
ccl <- consensusClusters(cageset)
ccl$id <- paste(ccl$chr, ccl$start, ccl$end, ccl$strand, sep="_")


ccl.df <- cageset@consensusClustersTpmMatrix
rownames(ccl.df) <- ccl$id
write.table(as.data.frame(ccl.df), "mouse_mapt_CAGE_clusters.txt", sep="\t", quote=F, col.names = NA)


###
# ChIP-seeker analysis
###
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# make peak file
peakFile <- data.frame(chr = ccl$chr, start = ccl$start, end = ccl$end)
peakFile <- cbind(peakFile, ccl.df)
write.table(peakFile, "mouse_mapt_CAGE_cluster_peaks.bed", sep="\t", quote=F, row.names = F, col.names = F)

# generate peak object
peaks <- readPeakFile("mouse_mapt_CAGE_cluster_peaks.bed")
peaks_gr <- GRanges(peakFile, strand = Rle(ccl$strand))


# define promoters around 3K;
promoter_3k <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000, by="transcript")
length(promoter_3k)


# Get a tag-matrix around that promoter region;
peaks_gr_tagmat_3k <- getTagMatrix(peaks_gr, windows=promoter_3k)

# Gene Annotation within the range of +/- 3kbases;
print("Peak annotation ...")
peaks_gr_3k <- annotatePeak(peaks_gr, tssRegion=c(-3000, 3000),
                            TxDb=txdb,
                            level = "gene",
                            assignGenomicAnnotation = TRUE,
                            genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"),
                            annoDb="org.Mm.eg.db")
plotAnnoPie(peaks_gr_3k)
plotAnnoBar(peaks_gr_3k)

upsetplot(peaks_gr_3k)

## Try to do peak annotation for all the ctss
?an

# Create gene-wise count table
counts <- as.data.frame(peaks_gr_3k)
genes <- counts$ENSEMBL
counts <- counts[, grepl("sample_", colnames(counts))]

counts <- aggregate(counts, by=list(genes), sum)
rownames(counts) <- counts$Group.1
counts <- counts[, -1]

#####
# PCA Plotting 
#####
md <- read.csv("rimod_mouse_mapt_meteadata.csv", sep=";")
md$sample <- paste0("sample_", md$sample_id)

# Perform PCA
library(edgeR)
library(ggfortify)
cpm.df <- cpm(counts)

pca_df <- prcomp(t(cpm.df))
autoplot(pca_df, shape=F)


# remove outliers
outs <- c("sample_16196", "sample_16195", "sample_16369", "sample_17754")
cpm.df.no_out <- cpm.df[, !colnames(cpm.df) %in% outs]
pca_df <- prcomp(t(cpm.df.no_out))
autoplot(pca_df, shape=F)

df <- data.frame(PC1=pca_df$x[,1], PC2=pca_df$x[,2], sample=rownames(pca_df$x))
df <- merge(df, md, by.x="sample", by.y="sample")
df$age <- factor(df$age)

ggplot(df, aes(x=PC1, y=PC2, color=group, shape=group)) +
  geom_point(size=10)

ggplot(df, aes(x=PC1, y=PC2, color=genotype, shape=age)) +
  geom_point(size=10)


ggplot(df, aes(x=PC1, y=PC2, color=age)) +
  geom_point(size=10)

ggplot(df, aes(x=PC1, y=PC2, color=sex)) +
  geom_point(size=10)


plot(pca_df$x[,1], pca_df$x[,2])

test <- removeBatchEffect(cpm.df.no_out, batch = df$sex)
