##############
# methylation analysis for plotting
##############
library(stringr)
setwd("/Users/kevin/dzne/rimod_package/frontal_methylation_0818/")


# Read data
mapt <- read.table("/Users/kevin/dzne/rimod_package/frontal_methylation_0818/DMPs_mapt.ndc_quant.txt", sep="\t", header=T)
grn <- read.table("/Users/kevin/dzne/rimod_package/frontal_methylation_0818/DMPs_grn.ndc_quant.txt", sep="\t", header=T)
c9 <- read.table("/Users/kevin/dzne/rimod_package/frontal_methylation_0818/DMPs_c9orf72.ndc_quant.txt", sep="\t", header=T)


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


# MAPT
mapt <- mapt[mapt$adj.P.Val <= 0.05,]
# up CpGs
mapt.up <- mapt[mapt$logFC > 0,]
genes <- as.character(mapt.up$GencodeBasicV12_NAME)
group <- as.character(mapt.up$GencodeBasicV12_Group)
genes <- filterGenes(genes, group)
write.table(genes, "MAPT_UP_CpGs_genes.txt", col.names = F, row.names = F, quote=F)

# down CpGs
mapt.down <- mapt[mapt$logFC < 0,]
genes <- as.character(mapt.down$GencodeBasicV12_NAME)
group <- as.character(mapt.down$GencodeBasicV12_Group)
genes <- filterGenes(genes, group)
write.table(genes, "MAPT_DOWN_CpGs_genes.txt", col.names = F, row.names = F, quote=F)
#==================##

# GRN
grn <- grn[grn$adj.P.Val <= 0.05,]
# up CpGs
grn.up <- grn[grn$logFC > 0,]
genes <- as.character(grn.up$GencodeBasicV12_NAME)
group <- as.character(grn.up$GencodeBasicV12_Group)
genes <- filterGenes(genes, group)
write.table(genes, "GRN_UP_CpGs_genes.txt", col.names = F, row.names = F, quote=F)

# down CpGs
grn.down <- grn[grn$logFC < 0,]
genes <- as.character(grn.down$GencodeBasicV12_NAME)
group <- as.character(grn.down$GencodeBasicV12_Group)
genes <- filterGenes(genes, group)
write.table(genes, "GRN_DOWN_CpGs_genes.txt", col.names = F, row.names = F, quote=F)
#==================##


# C9orf72
grn <- c9[c9$adj.P.Val <= 0.05,]
# up CpGs
grn.up <- grn[grn$logFC > 0,]
genes <- as.character(grn.up$GencodeBasicV12_NAME)
group <- as.character(grn.up$GencodeBasicV12_Group)
genes <- filterGenes(genes, group)
write.table(genes, "C9orf72_UP_CpGs_genes.txt", col.names = F, row.names = F, quote=F)

# down CpGs
grn.down <- grn[grn$logFC < 0,]
genes <- as.character(grn.down$GencodeBasicV12_NAME)
group <- as.character(grn.down$GencodeBasicV12_Group)
genes <- filterGenes(genes, group)
write.table(genes, "C9orf72_DOWN_CpGs_genes.txt", col.names = F, row.names = F, quote=F)
#==================##

