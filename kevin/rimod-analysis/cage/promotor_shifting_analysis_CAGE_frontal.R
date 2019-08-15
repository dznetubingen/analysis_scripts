####################################
# Shifting promotor analysis (frontal)
#####################################
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
setwd("~/rimod/CAGE/cage_analysis/promotor_shifting/frontal/")


## Create GRanges object from shifting result
createGR <- function(mat){
  gr <- GRanges(seqnames = mat$chr,
                strand = mat$strand,
                ranges = IRanges(start = mat$start,
                                 end = mat$end))
  return(gr)
}

# Define cutoffs
fdr <- 0.05
score_cutoff <- 0

# Define TxDB
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Load results
mapt <- read.table("mapt_shifting_promotors.txt", sep="\t", header=T)
grn <- read.table("grn_shifting_promotors.txt", sep="\t", header=T)
c9 <- read.table("c9_shifting_promotors.txt", sep="\t", header=T)

# Apply cutoffs
mapt <- mapt[mapt$fdr.KS <= fdr,]
mapt <- mapt[mapt$shifting.score > score_cutoff,]
grn <- grn[grn$fdr.KS <= fdr,]
grn <- grn[grn$shifting.score > score_cutoff,]
c9 <- c9[c9$fdr.KS <= fdr,]
c9 <- c9[c9$shifting.score > score_cutoff,]

# Create GR objects
mapt.gr <- createGR(mapt)
grn.gr <- createGR(grn)
c9.gr <- createGR(c9)

# Annotate shifting promotors
anno.mapt <- annotatePeak(mapt.gr, TxDb = txdb, annoDb = 'org.Hs.eg.db')
anno.grn <- annotatePeak(grn.gr, TxDb = txdb, annoDb = 'org.Hs.eg.db')
anno.c9 <- annotatePeak(c9.gr, TxDb = txdb, annoDb = 'org.Hs.eg.db')

mapt.df <- as.data.frame(anno.mapt)
grn.df <- as.data.frame(anno.grn)
c9.df <- as.data.frame(anno.c9)

# Save results
write.table(mapt.df$ENSEMBL, "mapt_promotor_shifting_genes_fro.txt", quote=F, row.names=F)
write.table(grn.df$ENSEMBL, "grn_promotor_shifting_genes_fro.txt", quote=F, row.names=F)
write.table(c9.df$ENSEMBL, "c9_promotor_shifting_genes_fro.txt", quote=F, row.names=F)


