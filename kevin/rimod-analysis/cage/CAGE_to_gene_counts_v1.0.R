########################################
## CAGE data analysis script ###########
#
# Script to generate Transcript and Genewise count tables given BAM files from CAGE-seq data
# The script is mostly following the same steps as the CAGE pipeline (Rscripts 4 and 5)
# NOTE: this sript uses the BAM files directly, and can eat a lot of memory. Be aware of this.
#
# Required input:
# - the bam files of interest
# - a GTF file for Gene annotation
# - directory where to perform the analysis
# The genome is hardcoded to HG39 (BSgenome.Hsapiens.UCSC.hg38)
#
# Author: Kevin Menden

#======================#
# Argument parsing
# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 4) {
  args <- c("--help")
}
## Help section
if("-h" %in% args ||"--help" %in% args) {
  cat("
        Please pass 4 arguments to the script in the following order:
        <bam_directory> the directory containin the BAM files for analysis
        <gtf_file> the GTF file
        <analysis_dir> the directory to store the results in (a subdirectory will be created)
        <tss_range> the range (in bp) for gene assignment around each TSS")
  stop()
}

#==== Hard coded section =======================================#
script_name = "CAGE_to_gene_counts_v1.0.R"
date = Sys.Date()
current_time = gsub(":", ".", gsub(" ", "_", Sys.time()))
# The genome
bsgenome = "BSgenome.Hsapiens.UCSC.hg38"
#================================================================#

# Load libraries
library(CAGEr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(ChIPseeker)
library(org.Hs.eg.db)
options(stringsAsFactors = F)

####### For testing purposes hard-coded argus #############
# Parse arguments
bam_directory = args[1]
gtf_file = args[2]
analysis_dir = args[3]
tss_range = as.numeric(args[4])


setwd(analysis_dir)
# Create sub-folder for current analysis
dir.create(paste("CAGE_gene_counts", current_time, sep="_"))
setwd(paste("CAGE_gene_counts", current_time, sep="_"))

# Save parameters in config file
params <- c(current_time, script_name, bam_directory, gtf_file, analysis_dir, as.character(tss_range))
param.names <- c("Time", "Script_name", "BAM_directory", "GTF_fiile", "Analysis_dir", "TSS_range")
params <- data.frame(param.name = param.names, param = params)
write.table(params, paste("config_file", current_time, sep="_"), quote = F, row.names = F, sep="\t")


### Data preparation
### List all bam files and extract sample labels
paths_to_bams<- list.files(bam_directory, full.names=TRUE, pattern = "*.bam$")
sample_labels <- list.files(bam_directory, pattern = "'.bam$")
sample_labels <- gsub(".bam", "", sample_labels)
sample_labels <- paste("sample",sample_labels, sep="_")

#==== Data loading and analysis ==================#

## Create CAGEset object
print("Generating CAGEset")
myCAGEset <- new("CAGEset",
                 genomeName = bsgenome,
                 inputFiles = paths_to_bams,
                 inputFilesType = "bam",
                 sampleLabels = sample_labels)



### Create CTSS files 
print("Creating CTSS files")
getCTSS(myCAGEset, removeFirstG = TRUE, correctSystematicG = TRUE)

## create ctss file
ctss <- CTSStagCount(myCAGEset)

## force the coordinate column to be integer, otherwise it creates
## problems in the step that creates the count table to be fed to DE
ctss[,2] <- as.integer(ctss[,2])

## save the ctss to a file with appropriate name
write.table(ctss, file=paste("cage_ctss_counts_", current_time, ".ctss", sep=""), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)


### Normalization
# no normalization --> keep raw counts so normalization can be done later
print("Normalize counts")
normalizeTagCount(myCAGEset, method="none")

### Clustering
# create Tag clusters
print("Clustering...")
clusterCTSS(object = myCAGEset, threshold = 10, thresholdIsTpm = FALSE,
            nrPassThreshold = 2, method = "distclu", maxDist =20,
            removeSingletons = TRUE)

### Cluster Aggregation
# aggregate cluster to consensus clusters
print("Aggregating cluters ...")
aggregateTagClusters(myCAGEset, tpmThreshold = 0, qLow = NULL, qUp = NULL, maxDist = 100)

countTable <- myCAGEset@consensusClustersTpmMatrix

for_deseq_consensus_cluster <- consensusClusters(myCAGEset)

for_deseq_consensus_cluster$id <- paste(for_deseq_consensus_cluster$chr,
                                        for_deseq_consensus_cluster$start, for_deseq_consensus_cluster$end,
                                        for_deseq_consensus_cluster$strand, sep="_")

rownames(countTable) <- for_deseq_consensus_cluster$id
print("Writing final table ...")
write.table(countTable, file=paste("cage_clustered_tss_counts_", current_time, ".txt", sep=""), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)

##===============================================##
# Annotation of TSS sites as transcripts          #

# Bring in bed-file format
countBed <- myCAGEset@consensusClustersTpmMatrix
consClus <- consensusClusters(myCAGEset)
chr <- consClus$chr
start <- consClus$start
end <- consClus$end
df <- data.frame(chr, start, end)
countBed <- cbind(df, countBed)
peak_file <- paste("cage_consClusters_peak_file_",current_time, ".bed", sep="")
write.table(countBed, peak_file, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)

# Generate peak object
print("Peak calling ...")
peaks <- readPeakFile(peak_file)
# Transform to GenomicRange object
peaks_gr <- GRanges(peaks, strand = Rle(strand(consClus$strand)))

# Creat TxDB genome
gencode_TxDb <- makeTxDbFromGFF(gtf_file, format = "gtf", dataSource = "GENCODE",  organism = "Homo sapiens")

# define promoters around 3K;
promoter_3k <- getPromoters(TxDb=gencode_TxDb, upstream=tss_range, downstream=tss_range, by="transcript")
length(promoter_3k)

# Get a tag-matrix around that promoter region;
peaks_gr_tagmat_3k <- getTagMatrix(peaks_gr, windows=promoter_3k)

# Plots Tags with heatmap
print("Plotting tag matrix ...")
pdf(paste("tag_heatmap_peaks_3k_", current_time,".pdf", sep=""))
tagHeatmap(peaks_gr_tagmat_3k,  xlim=c(-tss_range, tss_range), xlab="Genomic range", ylab="Transcripts", color="red")
dev.off()
# Plot average profile heatmaps
pdf(paste("tag_avgprf_heatmap_peaks_3k_", current_time,".pdf", sep=""))
plotAvgProf(peaks_gr_tagmat_3k,  xlim=c(-tss_range, tss_range), xlab="Genomic range", ylab="Read Count Frequency")
dev.off()

## Peak annotation

#=============== Genewise =======================#
# Gene Annotation within the range of +/- 3kbases;
print("Peak annotation ...")
peaks_gr_3k <- annotatePeak(peaks_gr, tssRegion=c(-tss_range, tss_range), 
                            TxDb=gencode_TxDb, 
                            level = "gene", 
                            assignGenomicAnnotation = TRUE, 
                            genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"), 
                            annoDb="org.Hs.eg.db", 
                            addFlankGeneInfo = FALSE, flankDistance = 5000, verbose = TRUE)


### Plot the results of the annotation
pdf(paste("annotation_gene_piechart_", current_time,".pdf", sep=""))
plotAnnoPie(peaks_gr_3k)
dev.off()
pdf(paste("annotation_gene_barplot_", current_time,".pdf", sep=""))
plotAnnoBar(peaks_gr_3k)
dev.off()
pdf(paste("annotation_gene_upsetplot_", current_time,".pdf", sep=""))
upsetplot(peaks_gr_3k)
dev.off()
pdf(paste("annotation_gene_vennpie_", current_time,".pdf", sep=""))
vennpie(peaks_gr_3k)
dev.off()
pdf(paste("annotation_gene_disttoTSS_", current_time,".pdf", sep=""))
plotDistToTSS(peaks_gr_3k, title="Distribution of CAGE clusters relative to TSS")
dev.off()

print("Sample name assignment ...")

### Sample name assignment for Gene count table
# Convert annotated object to dataframe
peaks_gr_3k_df <- as.data.frame(peaks_gr_3k)
# Generate genewise count table
peak_count_table <- data.frame(geneId = peaks_gr_3k_df$geneId,
                               ENTREZ = peaks_gr_3k_df$ENTREZID,
                               SYMBOL = peaks_gr_3k_df$SYMBOL)
# Get number of samples
nSamples <- ncol(countBed) - 3
counts <- peaks_gr_3k_df[,6:(5+nSamples)]
peak_count_table <- cbind(peak_count_table, counts)

## Assign sample columns back
for (i in 1:nSamples) {
  col_peak <- peak_count_table[,i+3]
  col_bed <- countBed[,i+3]
  if (all(col_peak == col_bed)){
    colnames(peak_count_table)[i+3] <- colnames(countBed[i+3])
  }
}

#============== Transcript wise ===================================#
# Transcript annotation
print("Peak annotation ...")
peaks_gr_3k_transcript <- annotatePeak(peaks_gr, tssRegion=c(-tss_range, tss_range), 
                                       TxDb=gencode_TxDb, 
                                       level = "transcript", 
                                       assignGenomicAnnotation = TRUE, 
                                       genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"), 
                                       annoDb="org.Hs.eg.db", 
                                       addFlankGeneInfo = FALSE, flankDistance = 5000, verbose = TRUE)
### Sample name assignment for Transcript count table
# Convert annotated object to dataframe
peaks_gr_3k_transcript_df <- as.data.frame(peaks_gr_3k_transcript)
# Generate genewise count table
peak_count_table_transcript <- data.frame(geneId = peaks_gr_3k_transcript_df$geneId,
                               ENTREZ = peaks_gr_3k_transcript_df$ENTREZID,
                               SYMBOL = peaks_gr_3k_transcript_df$SYMBOL)
# Get number of samples
nSamples <- ncol(countBed) - 3
counts_transcript <- peaks_gr_3k_transcript_df[,6:(5+nSamples)]
peak_count_table_transcript <- cbind(peak_count_table_transcript, counts_transcript)

## Assign sample columns back
for (i in 1:nSamples) {
  col_peak <- peak_count_table_transcript[,i+3]
  col_bed <- countBed[,i+3]
  if (all(col_peak == col_bed)){
    colnames(peak_count_table_transcript)[i+3] <- colnames(countBed[i+3])
  }
}

# Saving results
write.table(peak_count_table, paste("cage_gene_counts_", current_time, ".txt", sep=""), sep="\t", quote=F, row.names=F)
write.table(peak_count_table_transcript, paste("cage_transcript_counts_", current_time, ".txt", sep=""), sep="\t", quote=F, row.names=F)

print("Script finished")
