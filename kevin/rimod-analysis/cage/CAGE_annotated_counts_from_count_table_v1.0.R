########################################
## CAGE data analysis script ###########
#
# Script to generate Transcript and Genewise count tables given CTSS files from CAGE-seq data
# The script is mostly following the same steps as the CAGE pipeline (Rscripts 4 and 5)

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
      <count_table> the count table created from the CTSS file 
      <gtf_file> the GTF file
      <analysis_dir> the directory to store the results in (a subdirectory will be created)
      <tss_range> the range (in bp) for gene assignment around each TSS
      <annot_level> one of transcript | gene")
  stop()
}

#==== Hard coded section =======================================#
script_name = "CAGE_annotated_counts_from_count_table_v1.0.R"
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
count_table = args[1]
gtf_file = args[2]
analysis_dir = args[3]
tss_range = as.numeric(args[4])
annot_level = "transcript"



setwd(analysis_dir)
# Create sub-folder for current analysis
dir.create(paste("CAGE_annotated_counts_from_ct", current_time, sep="_"))
setwd(paste("CAGE_annotated_counts_from_ct", current_time, sep="_"))
dir.create("plots")

# Save parameters in config file
params <- c(current_time, script_name, count_table, gtf_file, analysis_dir, as.character(tss_range), annot_level)
param.names <- c("Time", "Script_name", "CTSS_directory", "GTF_fiile", "Analysis_dir", "TSS_range", "Annotation_level")
params <- data.frame(param.name = param.names, param = params)
write.table(params, paste("config_file", current_time, sep="_"), quote = F, row.names = F, sep="\t")


### Data preparation

# Read count table and transform to bed peak file
countTable <- read.table(count_table, sep="\t", header = T)
bed.info <- rownames(countTable)
chr <- as.character(sapply(bed.info, function(x) {strsplit(x, split="_")[[1]][[1]]}))
start <- as.character(sapply(bed.info, function(x) {strsplit(x, split="_")[[1]][[2]]}))
end <- as.character(sapply(bed.info, function(x) {strsplit(x, split="_")[[1]][[3]]}))
strand <- as.character(sapply(bed.info, function(x) {strsplit(x, split="_")[[1]][[4]]}))

peakFile <- data.frame(chr, start, end)
peakFile <- cbind(peakFile, countTable)

# Creat TxDB genome
gencode_TxDb <- makeTxDbFromGFF(gtf_file, format = "gtf", dataSource = "GENCODE",  organism = "Homo sapiens")



#==== Data loading and analysis ==================#
# Save the peak file in bed format
peak_file <- paste("cage_consClusters_peak_file_",current_time, ".bed", sep="")
write.table(peakFile, peak_file, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)

# Generate peak object
print("Peak calling ...")
peaks <- readPeakFile(peak_file)
# Transform to GenomicRange object
peaks_gr <- GRanges(peakFile, strand = Rle(strand))

# define promoters around 3K;
promoter_3k <- getPromoters(TxDb=gencode_TxDb, upstream=tss_range, downstream=tss_range, by="transcript")
length(promoter_3k)

# Get a tag-matrix around that promoter region;
peaks_gr_tagmat_3k <- getTagMatrix(peaks_gr, windows=promoter_3k)

# Plots Tags with heatmap
print("Plotting tag matrix ...")
pdf(paste("plots/tag_heatmap_peaks_3k_", current_time,".pdf", sep=""))
tagHeatmap(peaks_gr_tagmat_3k,  xlim=c(-tss_range, tss_range), xlab="Genomic range", ylab="Transcripts", color="red")
dev.off()
# Plot average profile heatmaps
pdf(paste("plots/tag_avgprf_heatmap_peaks_3k_", current_time,".pdf", sep=""))
plotAvgProf(peaks_gr_tagmat_3k,  xlim=c(-tss_range, tss_range), xlab="Genomic range", ylab="Read Count Frequency")
dev.off()

  
#=============== Annotation =======================#
# Gene Annotation within the range of +/- 3kbases;
print("Peak annotation ...")
peaks_gr_3k <- annotatePeak(peaks_gr, tssRegion=c(-tss_range, tss_range),
                              TxDb=gencode_TxDb,
                              level = annot_level,
                              assignGenomicAnnotation = TRUE,
                              genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"),
                              annoDb="org.Hs.eg.db",
                              addFlankGeneInfo = FALSE, flankDistance = 5000, verbose = TRUE)
  
  
### Plot the results of the annotation
pdf(paste("plots/annotation_gene_piechart_", current_time,".pdf", sep=""))
plotAnnoPie(peaks_gr_3k)
dev.off()
pdf(paste("plots/annotation_gene_barplot_", current_time,".pdf", sep=""))
plotAnnoBar(peaks_gr_3k)
dev.off()
pdf(paste("plots/annotation_gene_upsetplot_", current_time,".pdf", sep=""))
upsetplot(peaks_gr_3k)
dev.off()
pdf(paste("plots/annotation_gene_vennpie_", current_time,".pdf", sep=""))
vennpie(peaks_gr_3k)
dev.off()
pdf(paste("plots/annotation_gene_disttoTSS_", current_time,".pdf", sep=""))
plotDistToTSS(peaks_gr_3k, title="Distribution of CAGE clusters relative to TSS")
dev.off()

# Save annotated DF
annoDF <- as.data.frame(peaks_gr_3k)
write.table(annoDF, paste("annoated_peaks_", current_time, ".txt", sep=""), sep="\t", row.names=F, quote = F)


  
## Generate gene-level and trascript-level count table - Aggregate counts for genes
# remove chrM
nSamples <- ncol(countTable)
annoDF.noM <- annoDF[!annoDF$seqnames == "chrM",]

gene.counts <- data.frame(geneId = annoDF.noM$geneId)
gene.counts <- cbind(gene.counts, annoDF.noM[,6:(5+nSamples)])


# Merge genes
all.genes <- gene.counts$geneId
all.genes <- all.genes[!duplicated(all.genes)]
acc.df <- gene.counts[1,]
for (g in all.genes) {
  gdf <- gene.counts[gene.counts$geneId == g,]
  gdf.cts <- gdf[,2:ncol(gdf)]
  g.counts <- apply(gdf.cts, 2, sum)
  acc.df <- rbind(acc.df, gdf[1,])
}
acc.df <- acc.df[-1,]

write.table(acc.df, paste("merged_annot_gene_counts_", current_time, ".txt", sep=""), sep="\t", row.names = F, quote=F)

# Merge transcripts
transcript.counts <- data.frame(transcriptId = annoDF.noM$transcriptId)
transcript.counts <- cbind(transcript.counts, annoDF.noM[,6:(5+nSamples)])
all.transcripts <- transcript.counts$transcriptId
all.transcripts <- all.transcripts[!duplicated(all.transcripts)]
# aggreagte
acc.df <- transcript.counts[1,]
for (g in all.transcripts) {
  gdf <- transcript.counts[transcript.counts$transcriptId == g,]
  gdf.cts <- gdf[,2:ncol(gdf)]
  g.counts <- apply(gdf.cts, 2, sum)
  acc.df <- rbind(acc.df, gdf[1,])
}
acc.df <- acc.df[-1,]
write.table(acc.df, paste("merged_annot_transcript_counts_", current_time, ".txt", sep=""), sep="\t", row.names = F, quote=F)

print("Script finished")
