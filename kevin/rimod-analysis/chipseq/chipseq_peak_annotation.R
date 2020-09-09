#############
# Correlate CAGE-seq, ChIP-seq, RNA-seq
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
setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/chipseq/brain_cemil/")

# Load data
cts <- read.table("results/bwa/mergedLibrary/macs/narrowPeak/consensus/H3K4me3/deseq2/H3K4me3.consensus_peaks.featureCounts.txt", header=T, sep="\t")
peak_id = paste(cts$Geneid, cts$Chr, cts$Start, cts$End, cts$Strand, sep="_")
rownames(cts) <- peak_id
cts <- cts[, c(-1, -2, -3, -4, -5, -6)]
colnames(cts) <- str_split(colnames(cts), pattern="[.]", simplify = T)[,1]
##
# Create Metadata file
##
chip.md <- read.csv("../design_chipseq_frontal_neuron_cemil.csv")
chip.md$sample <- str_split(chip.md$fastq_1, pattern="/", simplify = T)[,5]
chip.md <- chip.md[!grepl("Input", chip.md$group),]
group <- str_split(chip.md$sample, pattern="_", simplify = T)[,1]
chip.md$sample_id <- str_sub(group, start=1, end=5)
group[grepl("C9", group)] <- "C9"
group[grepl("GRN", group)] <- "GRN"
group[grepl("MAPT", group)] <- "MAPT"
group[grepl("NDC", group)] <- "NDC"
id <- str_split(chip.md$sample, pattern="_", simplify = T)[,7]
sample <- paste0("Brain_", group, "_R", id)
chip.md$sample <- sample


# rimod overall sample sheet
md <- read.csv("~/rimod/files/FTD_Brain_corrected.csv")
md$SAMPLEID <- str_pad(md$SAMPLEID, width = 5, side = "left", pad = "0")
md <- md[md$SAMPLEID %in% chip.md$sample_id,]
md <- md[md$REGION == "frontal",]
md <- data.frame(sample=md$SAMPLEID, age=md$AGE, sex=md$GENDER, mutation=md$GENE)

md <- merge(chip.md, md, by.x="sample_id", by.y="sample")

# match files
cts <- cts[,colnames(cts) %in% md$sample]
cts <- cts[, match(md$sample, colnames(cts))]

#========= end ChIP-seq loading ============#

# Generate a peak file
bed.info <- rownames(cts)
chr <- as.character(sapply(bed.info, function(x) {strsplit(x, split="_")[[1]][[3]]}))
start <- as.numeric(as.character(sapply(bed.info, function(x) {strsplit(x, split="_")[[1]][[4]]})))
end <- as.numeric(as.character(sapply(bed.info, function(x) {strsplit(x, split="_")[[1]][[5]]})))
strand <- as.character(sapply(bed.info, function(x) {strsplit(x, split="_")[[1]][[6]]}))

# remove them NAs
chr <- chr[!is.na(start)]
end <- end[!is.na(start)]
strand <- strand[!is.na(start)]
cts <- cts[!is.na(start),]
start <- start[!is.na(start)]

peakFile <- data.frame(chr, start, end)
peakFile <- cbind(peakFile, cts)

# Creat TxDB genome
gencode_TxDb <- makeTxDbFromGFF("~/resources/genomes/gencode.v34.annotation.gtf", format = "gtf", dataSource = "GENCODE",  organism = "Homo sapiens")


# Transform to GenomicRange object
peaks_gr <- GRanges(peakFile, strand = Rle(strand))

# define promoters around 3K;
promoter_3k <- getPromoters(TxDb=gencode_TxDb, upstream=2000, downstream=2000, by="transcript")
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

# Get annotates dataframe
annoDF <- as.data.frame(peaks_gr_3k)
# get correct colnames back
colnames(annoDF)[6:20] <- paste0("sample_", md$sample_id)


# test if order is correct
count_test <- as.numeric(cts[1,])
anno_test <- as.numeric(annoDF[1,])[6:20]
all(count_test == anno_test)
# [1] TRUE
# --> order is still correct and we can happily save and have our data

write.table(annoDF, "RiMod_genewise_annotated_ChIPseq_peaks.txt", sep="\t", quote=F)

# make the count table
count_table <- annoDF
genes <- count_table$geneId
count_table <- count_table[, 6:20]

test <- aggregate(count_table, by=list(genes), sum)
count_table <- test
rownames(count_table) <- count_table$Group.1
count_table <- count_table[,-1]
write.table(count_table, "RiMod_ChIPseq_annotated_count_table.txt", sep="\t", quote=F)

#================ end annotation chipseq ==============#

