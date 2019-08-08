##################
# Create BigWig files from bed files for use with CAGEfightR
##################
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/CAGEseq/")


#' Convers CTSS bed file to BigWig
#'
#' @param ctss_file path to ctss.bed file
#' @param genomeInfo genome info file created e.g. by "rtracklayer::SeqinfoForUCSCGenome("hg38")"
#' @return The function writes two files, one for each strand, (".plus.bw" and ".minus.bw"), in the same directory as the original file
#' @importFrom GenomicRanges GRanges
#' @importFrom rtracklayer import export
#' @export
convert_ctss_to_bigwig = function( ctss_file, genomeInfo, fname ) {
  tpm_bed = rtracklayer::import(ctss_file)
  tpm_bed = GenomicRanges::GRanges(tpm_bed , seqinfo = genomeInfo)
  tpm_bed_plus = tpm_bed[ tpm_bed@strand == "+", ]
  tpm_bed_minus = tpm_bed[ tpm_bed@strand == "-", ]
  rtracklayer::export(object = tpm_bed_plus , paste(fname , ".plus.bw" , sep = "") , format = "BigWig" )
  rtracklayer::export(object = tpm_bed_minus, paste(fname , ".minus.bw" , sep = "") , format = "BigWig" )
}

# genome info 
gi = rtracklayer::SeqinfoForUCSCGenome("hg38")

# Load data
inputFiles = list.files(path = "bed_files/", pattern = ".bed$", full.names=TRUE)

# Only perform analysis of frontal samples
inputFiles <- inputFiles[grepl("_fro_", inputFiles)]


for (f in inputFiles){
  fname = gsub("_no_chrM.ctss.bed", "", basename(f))
  
  # Load Bed file back in and transform to BigWig Format
  convert_ctss_to_bigwig(f, genomeInfo = gi, fname = fname)
}









