##############
# methylation analysis for plotting
##############
library(stringr)
setwd("~/rimod/paper/figures/figure4/")


# Read data
mapt <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_mapt.ndc_quant.txt", sep="\t", header=T)
grn <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_grn.ndc_quant.txt", sep="\t", header=T)
c9 <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_c9orf72.ndc_quant.txt", sep="\t", header=T)

mapt <- mapt[mapt$adj.P.Val <= 0.05,]
grn <- grn[grn$adj.P.Val <= 0.05,]

# order by fold change
mapt <- mapt[order(abs(mapt$logFC), decreasing = T),]
mapt <- mapt[order(mapt$adj.P.Val),]
grn <- grn[order(grn$adj.P.Val),]

mapt.up <- mapt[mapt$logFC > 0,]
mapt.down <- mapt[mapt$logFC < 0,]

genes <- as.character(mapt.up$GencodeBasicV12_NAME)
group <- as.character(mapt.up$GencodeBasicV12_Group)
# subset

# Filter the gene list
filterGenes <- function(genes, group){
  genes <- str_split(genes, pattern=";", simplify = T)
  group <- str_split(group, pattern=";", simplify = T)

  # filter exon CpGs
  noexon <- str_detect(group, "Exon", negate = T)
  genes <- genes[noexon]
  group <- group[noexon]
  
  # filter 3'UTR CpGs
  no3p <- str_detect(group, "3'UTR", negate = T)
  genes <- genes[no3p]
  group <- group[no3p]
  
  # filter out duplicates and empty genes
  genes <- genes[!genes == ""]
  genes <- genes[!duplicated(genes)]
  
  return(genes)
}





