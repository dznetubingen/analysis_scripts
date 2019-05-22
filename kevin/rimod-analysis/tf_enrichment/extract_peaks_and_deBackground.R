#################################################################################################
## Script to extract DE CAGE peaks in bed format
## and extrant non-DE CAGE peaks in bed format to
## be used as background
##
## MODIFIED VERSION
## This version extracts the peaks that are DE in the opposite direction as background
##
## Usage: Rscript extract_peaks_and_background.R <deseq_deg_result.txt> <pval> <lfc>
#################################################################################################

args <- commandArgs(trailingOnly = TRUE)

deg_file <- as.character(args[1])
pval <- as.numeric(args[2])
lfc <- as.numeric(args[3])
print(deg_file)
print(paste("LFC: ", as.character(lfc), sep=""))
print(paste("P-val: ", as.character(pval), sep=""))
if (lfc > 0){
  print("Extracing up-regulated peaks")
} else {
  print("Extracting down-regulated peaks")
}


deg <- read.table(deg_file, sep="\t", row.names = 1, header = T)
deg.sig <- deg[deg$padj <= pval,]
if (lfc > 0){
  deg.bg <- deg.sig[deg.sig$log2FoldChange < -lfc,]
  deg.sig <- deg.sig[deg.sig$log2FoldChange >= lfc,]
  
} else{
  deg.bg <- deg.sig[deg.sig$log2FoldChange > abs(lfc),]
  deg.sig <- deg.sig[deg.sig$log2FoldChange <= lfc,]
}
dim(deg.sig)

createBedFile = function(df){
  genes <- rownames(df)
  chr <- as.character(sapply(genes, function(x){strsplit(x, split="_")[[1]][[1]]}))
  start <- as.character(sapply(genes, function(x){strsplit(x, split="_")[[1]][[2]]}))
  end <- as.character(sapply(genes, function(x){strsplit(x, split="_")[[1]][[3]]}))
  strand <- as.character(sapply(genes, function(x){strsplit(x, split="_")[[1]][[4]]}))
  new_df = data.frame(chr, start, end, strand)
  return(new_df)
}

sig.bed <- createBedFile(deg.sig)
bg.bed <- createBedFile(deg.bg)

write.table(sig.bed, "sig_genes.bed", sep="\t", quote=F, col.names = F, row.names = F)
write.table(bg.bed, "background.bed", sep="\t", quote=F, col.names = F, row.names = F)

