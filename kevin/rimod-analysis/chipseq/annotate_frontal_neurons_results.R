#############
# Annotate results from neuronal ChIP-seq analysis
# 07.09.2020
#############
library(DESeq2)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(ChIPseeker)
library(org.Hs.eg.db)



### Load and transform ChIP-seq stuff
setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/chipseq/frontal_sakib//")

# Load data
cts <- read.table("results/bwa/mergedLibrary/macs/narrowPeak/consensus/H3K4me3/deseq2/Neuron_C9orf72vsNeuron_NDC/Neuron_C9orf72vsNeuron_NDC.mLb.deseq2.FDR0.05.results.bed", sep="\t")
colnames(cts) <- c("chr", "start", "end", "interval", "lfc", "strand")
cts <- cts[, -4]
strand = cts$strand


peakFile <- data.frame(cts$chr, cts$start, cts$end)
peakFile$lfc <- cts$lfc

# Creat TxDB genome
gencode_TxDb <- makeTxDbFromGFF("~/resources/genomes/gencode.v34.annotation.gtf", format = "gtf", dataSource = "GENCODE",  organism = "Homo sapiens")


# Transform to GenomicRange object
peaks_gr <- GRanges(peakFile, strand = Rle(strand))

# define promoters around 3K;
promoter_3k <- getPromoters(TxDb=gencode_TxDb, upstream=2000, downstream=2000, by="gene")
length(promoter_3k)

# Get a tag-matrix around that promoter region;
peaks_gr_tagmat_3k <- getTagMatrix(peaks_gr, windows=promoter_3k)

print("Peak annotation ...")
peaks_gr_3k <- annotatePeak(peaks_gr, tssRegion=c(-2000, 2000),
                            TxDb=gencode_TxDb,
                            level = "gene",
                            assignGenomicAnnotation = TRUE,
                            genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"),
                            annoDb="org.Hs.eg.db",
                            addFlankGeneInfo = FALSE, flankDistance = 5000, verbose = TRUE)


pdf("annotation_gene_piechart.pdf")
plotAnnoPie(peaks_gr_3k)
dev.off()
pdf("annotation_gene_barplot.pdf")
plotAnnoBar(peaks_gr_3k)
dev.off()


# Save annotated DF
annoDF <- as.data.frame(peaks_gr_3k)
write.table(annoDF, "annoated_peaks.txt", sep="\t", row.names=F, quote = F)

genes <- as.character(na.omit(annoDF$SYMBOL))
genes <- genes[!duplicated(genes)]

write.table(genes, "C9orf72_chipseq_frontal_neurons_genes.txt", quote=F, row.names = F, col.names = F)
