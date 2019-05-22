########################################
## CAGE data analysis script ###########
#
# Script to generate Transcript and Genewise count tables given BAM files from CAGE-seq data
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
if(length(args) != 2) {
  args <- c("--help")
}
## Help section
if("-h" %in% args ||"--help" %in% args) {
  cat("
      Please pass 2 arguments to the script in the following order:
      <bam_directory> the directory containin the BAM files for analysis
      <analysis_dir> the directory to store the results in (a subdirectory will be created)")
  stop()
}

#==== Hard coded section =======================================#
script_name = "CAGE_create_ctss_files_from_bams_v1.0.R"
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
# args <- c("~/rimod/CAGE/frontal_bam_files/files_sub_4/",
#           "~/rimod/CAGE/cager_testing/")
bam_directory = args[1]
analysis_dir = args[2]


setwd(analysis_dir)
# Create sub-folder for current analysis
dir.create(paste("CAGE_ctss_files_", current_time, sep="_"))
setwd(paste("CAGE_ctss_files_", current_time, sep="_"))

# Save parameters in config file
params <- c(current_time, script_name, bam_directory, analysis_dir)
param.names <- c("Time","Script_name","BAM_directory", "Analysis_dir")
params <- data.frame(param.name = param.names, param = params)
write.table(params, paste("config_file", current_time, sep="_"), quote = F, row.names = F, sep="\t")


### Data preparation
### List all bam files and extract sample labels
paths_to_bams<- list.files(bam_directory, full.names=TRUE, pattern = "*.bam$")
sample_labels <- list.files(bam_directory, pattern = ".bam$")
sample_labels <- gsub(".bam", "", sample_labels)
sample_labels <- paste("sample",sample_labels, sep="_")

#==== Data loading and analysis ==================#
## Load every bam file separately and create a CTSS file 

for (i in 1:length(paths_to_bams)){
  bam <- paths_to_bams[i]
  sample <- sample_labels[i]
  print(paste("Processing ", sample, sep=""))
  cageset <- new("CAGEset",
                 genomeName = bsgenome,
                 inputFiles = bam,
                 inputFilesType = "bam",
                 sampleLabels = sample)
  
  getCTSS(cageset, removeFirstG = TRUE, correctSystematicG = TRUE)
  
  ## create ctss file
  ctss <- CTSStagCount(cageset)
  # Transform coordinate column to integer
  ctss[,2] <- as.integer(ctss[,2])
  ## save the ctss to a file with appropriate name
  write.table(ctss, file=paste(sample, ".ctss", sep=""), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
}

print("Script finished")
