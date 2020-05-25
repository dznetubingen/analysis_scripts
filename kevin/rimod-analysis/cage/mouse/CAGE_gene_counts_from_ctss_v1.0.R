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
if(length(args) != 6) {
  args <- c("--help")
}
## Help section
if("-h" %in% args ||"--help" %in% args) {
  cat("
      Please pass 4 arguments to the script in the following order:
      <ctss_directory> the directory containin the ctss files for analysis
      <gtf_file> the GTF file
      <analysis_dir> the directory to store the results in (a subdirectory will be created)
      <tss_range> the range (in bp) for gene assignment around each TSS
      <annot_level> one of transcript | gene
      <keep_singleton> TRUE or FALSE, whether singletons should be kept")
  stop()
}

#==== Hard coded section =======================================#
script_name = "CAGE_gene_counts_from_ctss_v1.0.R"
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
# args <- c("~/rimod/CAGE/cager_testing/CAGE_ctss_files__2018-03-02_10.46.45/",
#           "~/resources/genomes/GRCh38_v27_gencode/gencode.v27.annotation.gtf",
#           "~/rimod/CAGE/cager_testing/",
#           "3000",
#           "transcript",
#           "TRUE")
# Parse arguments
ctss_directory = args[1]
gtf_file = args[2]
analysis_dir = args[3]
tss_range = as.numeric(args[4])
annot_level = args[5]
keep_singletons = args[6] == "TRUE"


setwd(analysis_dir)
# Create sub-folder for current analysis
dir.create(paste("CAGE_gene_counts_from_ctss", current_time, sep="_"))
setwd(paste("CAGE_gene_counts_from_ctss", current_time, sep="_"))
dir.create("plots")

# Save parameters in config file
params <- c(current_time, script_name, ctss_directory, gtf_file, analysis_dir, as.character(tss_range), annot_level, keep_singletons)
param.names <- c("Time", "Script_name", "CTSS_directory", "GTF_fiile", "Analysis_dir", "TSS_range", "Annotation_level", "Keep_singletons")
params <- data.frame(param.name = param.names, param = params)
write.table(params, paste("config_file", current_time, sep="_"), quote = F, row.names = F, sep="\t")


### Data preparation
### List all bam files and extract sample labels
paths_to_ctss<- list.files(ctss_directory, full.names=TRUE, pattern = "*.ctss$")
sample_labels <- list.files(ctss_directory, pattern = "*ctss$")
sample_labels <- gsub(".ctss", "", sample_labels)

#==== Data loading and analysis ==================#

# Creat TxDB genome
gencode_TxDb <- makeTxDbFromGFF(gtf_file, format = "gtf", dataSource = "GENCODE",  organism = "Homo sapiens")

annoPeakList <- list()

### Loop over every single CTSS file
# parse it
# generate peak bed file
# call the peaks
# annotate the peaks

for (i in 1:length(paths_to_ctss)){
  ctss_path = paths_to_ctss[i]
  sample_lbl <- sample_labels[i]
  print(paste("#### PROCESSING ", sample_lbl, " ###", sep=""))
  
  myCAGEset <- new("CAGEset",
                   genomeName = bsgenome,
                   inputFiles = ctss_path,
                   inputFilesType = "ctss",
                   sampleLabels = sample_lbl)
  
  getCTSS(myCAGEset, removeFirstG = TRUE, correctSystematicG = TRUE)
  ## create ctss file
  ctss <- CTSStagCount(myCAGEset)
  ## force the coordinate column to be integer, otherwise it creates
  ## problems in the step that creates the count table to be fed to DE
  ctss[,2] <- as.integer(ctss[,2])
  ## save the ctss to a file with appropriate name
  write.table(ctss, file=paste(sample_lbl,"_cage_ctss_counts_", current_time, ".ctss", sep=""), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
  ### Normalization
  # no normalization --> keep raw counts so normalization can be done later
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
  write.table(countTable, file=paste(sample_lbl,"_cage_clustered_tss_counts_", current_time, ".txt", sep=""), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
  
  
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
  peak_file <- paste(sample_lbl, "_cage_consClusters_peak_file_",current_time, ".bed", sep="")
  write.table(countBed, peak_file, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  # Generate peak object
  print("Peak calling ...")
  peaks <- readPeakFile(peak_file)
  # Transform to GenomicRange object
  peaks_gr <- GRanges(peaks, strand = Rle(strand(consClus$strand)))
  
  # define promoters around 3K;
  promoter_3k <- getPromoters(TxDb=gencode_TxDb, upstream=tss_range, downstream=tss_range, by="transcript")
  length(promoter_3k)
  
  # Get a tag-matrix around that promoter region;
  peaks_gr_tagmat_3k <- getTagMatrix(peaks_gr, windows=promoter_3k)
  
  # Plots Tags with heatmap
  print("Plotting tag matrix ...")
  pdf(paste("plots/",sample_lbl, "_tag_heatmap_peaks_3k_", current_time,".pdf", sep=""))
  tagHeatmap(peaks_gr_tagmat_3k,  xlim=c(-tss_range, tss_range), xlab="Genomic range", ylab="Transcripts", color="red")
  dev.off()
  # Plot average profile heatmaps
  pdf(paste("plots/",sample_lbl, "_tag_avgprf_heatmap_peaks_3k_", current_time,".pdf", sep=""))
  plotAvgProf(peaks_gr_tagmat_3k,  xlim=c(-tss_range, tss_range), xlab="Genomic range", ylab="Read Count Frequency")
  dev.off()
  
  ## Peak annotation
  
  #=============== Genewise =======================#
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
  pdf(paste("plots/",sample_lbl, "_annotation_gene_piechart_", current_time,".pdf", sep=""))
  plotAnnoPie(peaks_gr_3k)
  dev.off()
  pdf(paste("plots/",sample_lbl, "_annotation_gene_barplot_", current_time,".pdf", sep=""))
  plotAnnoBar(peaks_gr_3k)
  dev.off()
  pdf(paste("plots/",sample_lbl, "_annotation_gene_upsetplot_", current_time,".pdf", sep=""))
  upsetplot(peaks_gr_3k)
  dev.off()
  pdf(paste("plots/",sample_lbl, "_annotation_gene_vennpie_", current_time,".pdf", sep=""))
  vennpie(peaks_gr_3k)
  dev.off()
  pdf(paste("plots/",sample_lbl, "_annotation_gene_disttoTSS_", current_time,".pdf", sep=""))
  plotDistToTSS(peaks_gr_3k, title="Distribution of CAGE clusters relative to TSS")
  dev.off()
  
  # Append to list
  annoPeakList$tmp <- peaks_gr_3k
  names(annoPeakList)[i] <- sample_lbl
}


# Create plots together
pdf(paste("plots/all_samples_annotation_gene_barplot_", current_time,".pdf", sep=""))
plotAnnoBar(annoPeakList)
dev.off()
pdf(paste("plots/all_samples_annotation_gene_disttoTSS_", current_time,".pdf", sep=""))
plotDistToTSS(annoPeakList)
dev.off()

### Save all peak files
for (i in 1:length(annoPeakList)){
  sample_name <- names(annoPeakList)[i]
  annoPeak <- as.data.frame(annoPeakList[[i]])
  write.table(annoPeak, paste(sample_name, "_annotated_peaks_", current_time,".txt", sep=""), sep="\t", quote=F, row.names = F)
}


### Collapse peaks to genewise counts
count_tables <- list()
noSamples <- length(annoPeakList)
print("Removing chrM and aggregating gene peaks ...")
for (i in 1:noSamples){
  print(paste("processing ", i, " of ", noSamples, sep=""))
  annoPeak <- as.data.frame(annoPeakList[[i]])
  sample_name <- names(annoPeakList)[i]
  
  # Remove chr M
  annoPeak <- annoPeak[!annoPeak$seqnames == "chrM",]
  
  genes <- as.character(annoPeak$geneId)
  genes <- genes[!duplicated(genes)]
  
  sample.df <- data.frame(geneId = "decoy", count = 1)
  for (gene in genes) {
    gdf <- annoPeak[annoPeak$geneId == gene,]
    counts <- sum(gdf$V4)
    gene.df <- data.frame(geneId = gene, count = counts)
    sample.df <- rbind(sample.df, gene.df) 
  }
  sample.df <- sample.df[-1,] # remove decoy row
  count_tables$tmp <- sample.df
  names(count_tables)[i] <- sample_name
}

### Merge count tables into one
# Get all genes that have been called in any samples
all_genes <- c()
for (i in 1:noSamples){
  ct <- count_tables[[i]]
  tmp <- as.character(ct$geneId)
  all_genes <- c(all_genes, tmp)
}
all_genes <- all_genes[!duplicated(all_genes)]

# Merging
print("Merging count tables ...")
final_count_table <- data.frame(geneId = all_genes, decoy = rep(0, length(all_genes)))
for (i in 1:noSamples) {
  print(paste("Merging ", i, " ...", sep=""))
  ct <- count_tables[[i]]
  final_count_table <- merge(final_count_table, ct, by.x="geneId", by.y="geneId", all.x = T)
  final_count_table$count[is.na(final_count_table$count)] <- 0 # set unavailable values to zero
  # Adjust colname to sample name
  colnames(final_count_table)[i+2] <- names(count_tables)[i]
}
final_count_table <- final_count_table[,-2]

print("Writing merged count table ...")

write.table(final_count_table, paste("merged_genewise_count_table_", current_time, ".txt", sep=""), sep="\t", row.names=F, quote=F)

print("Script finished")
