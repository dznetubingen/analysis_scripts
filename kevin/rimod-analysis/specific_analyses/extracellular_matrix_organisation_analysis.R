##################
# Analysis of the Extracellular Matrix Organization (REACTOME)
#################
library(pheatmap)
library(stringr)
library(viridis)
library(biomaRt)


ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

setwd("~/rimod/integrative_analysis/extracellular_matrix/")

# load expression
mat <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-08-12_07.58.35/deseq_vst_values_2019-08-12_07.58.35.txt",
                  sep="\t", header=T, row.names=1, check.names = F)
rownames(mat) <- str_split(rownames(mat), pattern="[.]", simplify = T)[,1]
colnames(mat) <- gsub("X", "", colnames(mat))

# load md
md <- read.table("~/rimod/RNAseq/rnaseq_frontal_md.txt", header=T)

# Load ubiquting pathway
ubi <- read.table("ECM_reactome.tsv", sep="\t", header=T)
genes <- str_split(ubi$MoleculeName, pattern=" ", simplify = T)[,2]
ubi$name <- genes
# Get ensembl IDs
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="hgnc_symbol", values=genes, mart=ensembl)
ubi <- merge(ubi, bm, by.x="name", by.y="hgnc_symbol")
ubi <- ubi[!duplicated(ubi$ensembl_gene_id),]
rownames(ubi) <- ubi$ensembl_gene_id


# Subset and order the expression matrix
mat <- mat[rownames(mat) %in% ubi$ensembl_gene_id,]
mat <- mat[,md$ids]
ubi <- ubi[rownames(mat),]
# change ensemble ID to hgnc
rownames(mat) <- ubi$name

# Generate heatmap for DEGs
anno_col <- data.frame(dc = md$mutated_gene)
rownames(anno_col) <- colnames(mat)

# GRN
deg.grn <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/GRN_cell_composition_filterd_DEGs.txt", sep="\t", header=T)
mat.grn <- mat[rownames(mat) %in% deg.grn$hgnc_symbol,]
pheatmap(mat.grn, scale="row", color = viridis(200), annotation_col = anno_col)

# MAPT
deg.mapt <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/MAPT_cell_composition_filterd_DEGs.txt", sep="\t", header=T)
mat.mapt <- mat[rownames(mat) %in% deg.mapt$hgnc_symbol,]
pheatmap(mat.mapt, scale="row", color = viridis(200), annotation_col = anno_col)

# C9orf72 (no overlap)
deg.c9 <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/C9orf72_cell_composition_filterd_DEGs.txt", sep="\t", header=T)
mat.c9 <- mat[rownames(mat) %in% deg.c9$hgnc_symbol,]
pheatmap(mat.c9, scale="row", color = viridis(200), annotation_col = anno_col)

###
# Check for alternative splicing
###

as.grn <- read.table("~/rimod/RNAseq/as_analysis/majiq/grn_AS_genes_dPSI_0.2.txt", sep="\t", header = T)
as.grn <- as.grn[as.grn$x %in% ubi$ensembl_gene_id,]

as.mapt <- read.table("~/rimod/RNAseq/as_analysis/majiq/mapt_AS_genes_dPSI_0.2.txt", sep="\t", header = T)
as.mapt <- as.character(as.mapt[as.mapt$x %in% ubi$ensembl_gene_id,])

#===================================#

###
# Check for promotor shifting
###
# grn
ps.grn <- read.table("~/rimod/CAGE/cage_analysis/promotor_shifting/frontal/grn_promotor_shifting_genes_fro.txt", header=T)
ps.grn <- as.character(ps.grn[ps.grn$x %in% ubi$ensembl_gene_id,])
ps.grn <- ubi[ps.grn,]
# mapt
ps.mapt <- read.table("~/rimod/CAGE/cage_analysis/promotor_shifting/frontal/mapt_promotor_shifting_genes_fro.txt", header=T)
ps.mapt <- as.character(ps.mapt[ps.mapt$x %in% ubi$ensembl_gene_id,])
ps.mapt <- ubi[ps.mapt,]


#==========================================#

###
# Check for overlap with miRNA targets
###

mir.grn <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/GRN_miRNA_target_edge_table.txt", 
                      sep="\t", header=T, stringsAsFactors = F)
mir.grn <- mir.grn[mir.grn$targets %in% ubi$name,]
table(mir.grn$targets)

# miR-4492 has lots of targets in pathway --> further check its function
mir.grn.all <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/GRN_miRNA_target_edge_table.txt", 
                          sep="\t", header=T, stringsAsFactors = F)
mir.4492 <- mir.grn.all[mir.grn.all$mirna == "hsa-miR-4492",]
write.table(mir.4492$targets, "GRN_miR.4492_targets.txt", quote=F, row.names=F, col.names = F)


mir.mapt <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/MAPT_miRNA_target_edge_table.txt", 
                       sep="\t", header=T, stringsAsFactors = F)
mir.mapt <- mir.mapt[mir.mapt$targets %in% ubi$name,]
table(mir.mapt$targets)


#============================================#


###
# Check for differential methylation at this locus
###

# Extrac gene accessions
getGenes <- function(x){
  genes <- as.character(x$GencodeCompV12_NAME)
  genes <- genes[!genes == ""]
  genes <- as.character(sapply(genes, function(y){strsplit(y, split="[;]")[[1]][[1]]}))
  genes <- genes[!duplicated(genes)]
  return(genes)
}

# MAPT
mapt.met <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_mapt.ndc_quant.txt", sep="\t")
mapt.met <- mapt.met[mapt.met$adj.P.Val <= 0.01,]
mapt.met.genes <- getGenes(mapt.met)
mapt.met.genes <- mapt.met.genes[mapt.met.genes %in% ubi$name]

# GRN
grn.met <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_grn.ndc_quant.txt", sep="\t")
grn.met <- grn.met[grn.met$P.Value <= 0.01,]
grn.met.genes <- getGenes(grn.met)
grn.met.genes <- grn.met.genes[grn.met.genes %in% ubi$name]

# C9orf72
c9.met <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_c9orf72.ndc_quant.txt", sep="\t")
c9.met <- c9.met[c9.met$P.Value <= 0.01,]
c9.met.genes <- getGenes(c9.met)
c9.met.genes <- c9.met.genes[c9.met.genes %in% ubi$name]

# Further look into some genes
mapt.pdp <- mapt.met[grepl("PDP2", mapt.met$GencodeBasicV12_NAME),]
grn.pdp <- grn.met[grepl("PDP2", grn.met$GencodeBasicV12_NAME),]
c9.pdp <- c9.met[grepl("PDP2", c9.met$GencodeBasicV12_NAME),]

#====================================================================#

###
# ChIP-seq
###

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(ChIPseeker)
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

# MAPT
chip.mapt <- read.table("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/chipseq/results_mapt_narrow/bwa/mergedLibrary/macs/narrowPeak/consensus/H3K4me3/deseq2/H3K4me3.consensus_peaks.results.txt",
                        sep="\t", header=T)
chip.mapt <- chip.mapt[chip.mapt$Brain_MAPTvsBrain_NDC.pvalue <= 0.05,]
df <- data.frame(chr=chip.mapt$Chr, start=chip.mapt$Start, end=chip.mapt$End, strand=chip.mapt$Strand,
                 pval=chip.mapt$Brain_MAPTvsBrain_NDC.pvalue, lfc=chip.mapt$Brain_MAPTvsBrain_NDC.log2FoldChange)
df <- na.omit(df)
gr <- makeGRangesFromDataFrame(df)
mapt.anno <- as.data.frame(annotatePeak(gr, TxDb=txdb, annoDb = "org.Hs.eg.db"))
mapt.anno <- mapt.anno[mapt.anno$SYMBOL %in% ubi$name,]

# GRN
chip.grn <- read.table("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/chipseq/results_grn_narrow/bwa/mergedLibrary/macs/narrowPeak/consensus/H3K4me3/deseq2/H3K4me3.consensus_peaks.results.txt",
                       sep="\t", header=T)
chip.grn <- chip.grn[chip.grn$Brain_GRNvsBrain_NDC.pvalue <= 0.05,]
df <- data.frame(chr=chip.grn$Chr, start=chip.grn$Start, end=chip.grn$End, strand=chip.grn$Strand,
                 pval=chip.grn$Brain_GRNvsBrain_NDC.pvalue, lfc=chip.grn$Brain_GRNvsBrain_NDC.log2FoldChange)
df <- na.omit(df)
gr <- makeGRangesFromDataFrame(df)
grn.anno <- as.data.frame(annotatePeak(gr, TxDb=txdb, annoDb = "org.Hs.eg.db"))
grn.anno <- grn.anno[grn.anno$SYMBOL %in% ubi$name,]

# C9or72
chip.c9 <- read.table("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/chipseq/results_c9orf72_narrow/bwa/mergedLibrary/macs/narrowPeak/consensus/H3K4me3/deseq2/H3K4me3.consensus_peaks.results.txt",
                      sep="\t", header=T)
chip.c9 <- chip.c9[chip.c9$Brain_C9vsBrain_NDC.pvalue <= 0.05,]
df <- data.frame(chr=chip.c9$Chr, start=chip.c9$Start, end=chip.c9$End, strand=chip.c9$Strand,
                 pval=chip.c9$Brain_C9vsBrain_NDC.pvalue, lfc=chip.c9$Brain_C9vsBrain_NDC.log2FoldChange)
df <- na.omit(df)
gr <- makeGRangesFromDataFrame(df)
c9.anno <- as.data.frame(annotatePeak(gr, TxDb=txdb, annoDb = "org.Hs.eg.db"))
c9.anno <- c9.anno[c9.anno$SYMBOL %in% ubi$name,]

#==================================#

####
# Check for TFs
####
# mapt
tf.mapt <- read.table("~/rimod/CAGE/cage_analysis/chea_tf_regulators/frontal/Integrated_meanRank_MAPT_up.tsv", sep="\t", header=T)
tf.mapt <- tf.mapt[1:50,]
gene_overlaps <- c()
for (i in 1:nrow(tf.mapt)) {
  ovl <- as.character(tf.mapt$Overlapping_Genes[i])
  ovl <- str_split(ovl, pattern=",")[[1]]
  ovl <- ovl[ovl %in% ubi$name]
  gene_overlaps <- c(gene_overlaps, length(ovl))
}
tf.mapt$wnt_overlap <- gene_overlaps
tf.mapt <- tf.mapt[order(tf.mapt$wnt_overlap, decreasing = T),]
#tf.mapt <- tf.mapt[1:20,]
tf.mapt <- tf.mapt[tf.mapt$TF %in% deg.mapt$hgnc_symbol,]

# grn
tf.grn <- read.table("~/rimod/CAGE/cage_analysis/chea_tf_regulators/frontal/Integrated_meanRank_GRN_up.tsv", sep="\t", header=T)
tf.grn <- tf.grn[1:50,]
gene_overlaps <- c()
for (i in 1:nrow(tf.grn)) {
  ovl <- as.character(tf.grn$Overlapping_Genes[i])
  ovl <- str_split(ovl, pattern=",")[[1]]
  ovl <- ovl[ovl %in% ubi$name]
  gene_overlaps <- c(gene_overlaps, length(ovl))
}
tf.grn$wnt_overlap <- gene_overlaps
tf.grn.wnt <- tf.grn[tf.grn$TF %in% ubi$name,]
tf.grn <- tf.grn[order(tf.grn$wnt_overlap, decreasing = T),]
#tf.grn <- tf.grn[1:20,]
tf.grn <- tf.grn[tf.grn$TF %in% deg.grn$hgnc_symbol,]

####
# Build iGraph
####
library(igraph)

####
# Build network for FTD-GRN
####

### grn
genes <- as.character(ubi$name)
deg <- deg.grn[deg.grn$hgnc_symbol %in% genes,]
# get all fold change values as well
exp.grn <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-08-12_07.58.35/deseq_result_grn.ndc_fro_2019-08-12_07.58.35.txt", sep="\t", header=T)
exp.grn <- exp.grn[exp.grn$X %in% ubi$ensembl_gene_id,]
exp.grn <- merge(exp.grn, ubi, by.x="X", by.y="ensembl_gene_id")
exp.grn <- exp.grn[, -ncol(exp.grn)]
# get miRNA fold changes
exp.mir <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_result_grn.ndc_frontal_smRNAseq.txt", sep="\t", header=T)

# get PPI connections
ppi <- read.table("ecm_string.tsv", sep="\t", header=T)
ppi <- ppi[, c(1,2,15)]

# Create edge list from PPIs
edges <- c()
for (i in 1:nrow(ppi)){
  e <- c(as.character(ppi[i,1]), as.character(ppi[i,2]))
  edges <- c(edges, e)
}

# Add miRNA edges
for (i in 1:nrow(mir.grn)){
  e <- c(as.character(mir.grn[i,1]), as.character(mir.grn[i,2]))
  edges <- c(edges, e)
}


# Add TF connections to network
for (i in 1:nrow(tf.grn)) {
  tf <- as.character(tf.grn$TF[i])
  targets <- str_split(tf.grn$Overlapping_Genes[i], pattern = ",")[[1]]
  targets <- targets[targets %in% ubi$name]
  for (tg in targets){
    e <- c(tf, as.character(tg))
    edges <- c(edges, e)
  }
}

g <- graph(edges=edges)


# Assign type
vnames <- V(g)$name
types = c()
vtypes <- c()
for (v in vnames) {
  type = ""
  if (grepl("hsa-", v)){
    type = "miRNA"
  }
  else if (v %in% tf.grn$TF){
    type = "TF"
  }
  else {
    type = "Gene"
  }
  vtypes <- c(vtypes, type)
}
V(g)$type <- vtypes

# Assign LFC
lfcs <- c()
for (v in vnames) {
  fc = 0
  if (v %in% exp.mir$X){
    tmp <- exp.mir[exp.mir$X == v,]
    fc <- as.numeric(tmp$log2FoldChange)
  }
  else if (v %in% exp.grn$name){
    tmp <- exp.grn[exp.grn$name == v,]
    fc <- as.numeric(tmp$log2FoldChange)
  }
  else if (v %in% deg.grn$hgnc_symbol){
    tmp <- deg.grn[deg.grn$hgnc_symbol == v,]
    fc <- as.numeric(tmp$log2FoldChange)
  }
  lfcs <- c(lfcs, fc)
}
V(g)$lfc <- lfcs

# Save file
file_name = "GRN_ECM_network.gml"
write_graph(g, file=file_name, format="gml")
print("Network creation succesfull.")

