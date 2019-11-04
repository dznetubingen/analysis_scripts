#################################
# Membrane trafficking sub-module investigation
# Check how different sub-modules from membrane trafficking are regulated
#################################
library(pheatmap)
library(stringr)
library(viridis)
library(biomaRt)

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

setwd("~/rimod/integrative_analysis/membrane_trafficking//")

# load expression
mat <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_vst_values_2019-10-23_13.33.11.txt",
                  sep="\t", header=T, row.names=1, check.names = F)
rownames(mat) <- str_split(rownames(mat), pattern="[.]", simplify = T)[,1]
colnames(mat) <- gsub("X", "", colnames(mat))

# load md
md <- read.table("~/rimod/RNAseq/rnaseq_frontal_md.txt", header=T)


# Function to load a pathway
loadPathway <- function(path){
  pw <- read.table(path, sep="\t", header=T)
  genes <- str_split(pw$MoleculeName, pattern=" ", simplify = T)[,2]
  pw$name <- genes
  # Get ensembl IDs
  bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="hgnc_symbol", values=genes, mart=ensembl)
  pw <- merge(pw, bm, by.x="name", by.y="hgnc_symbol")
  pw <- pw[!duplicated(pw$ensembl_gene_id),]
  rownames(pw) <- pw$ensembl_gene_id
  return(pw)
}


# load pathway
modules <- c("clathrin_mediated_endocytosis.tsv", "er_to_golgi_anterograde_transport.tsv", "escrt_pathway.tsv",
             "gap_junction_trafficking_and_regulation.tsv", "intra_golgi_and_retrograde_golgi_to_er_traffick.tsv",
             "rab_regulation_of_trafficking.tsv", "trans_golgi_network_vesicle_budding.tsv",
             "translocation_of_glut4_to_the_plasma_membrane.tsv")

mt <- loadPathway("membrane_trafficking_reactome.tsv")
mod1 <- loadPathway(modules[1])
mod2 <- loadPathway(modules[2])
mod3 <- loadPathway(modules[3])
mod4 <- loadPathway(modules[4])
mod5 <- loadPathway(modules[5])
mod6 <- loadPathway(modules[6])
mod7 <- loadPathway(modules[7])
mod8 <- loadPathway(modules[8])




# Subset and order the expression matrix
mat <- mat[rownames(mat) %in% mt$ensembl_gene_id,]
mat <- mat[,md$ids]
mt <- mt[rownames(mat),]
# change ensemble ID to hgnc
rownames(mat) <- mt$name

# Assign module IDs to pathway
mt$module <- rep("mod1", nrow(mt))
for (i in 1:nrow(mt)) {
  gene <- rownames(mt)[i]
  if (gene %in% rownames(mod1)){
    mt$module[i] <- "mod1"
  }
  else if (gene %in% rownames(mod2)){
    mt$module[i] <- "mod2"
  }
  else if (gene %in% rownames(mod3)){
    mt$module[i] <- "mod3"
  }
  else if (gene %in% rownames(mod4)){
    mt$module[i] <- "mod4"
  }
  else if (gene %in% rownames(mod5)){
    mt$module[i] <- "mod5"
  }
  else if (gene %in% rownames(mod6)){
    mt$module[i] <- "mod6"
  }
  else if (gene %in% rownames(mod7)){
    mt$module[i] <- "mod7"
  }
  else if (gene %in% rownames(mod8)){
    mt$module[i] <- "mod8"
  }
}


# MAPT
deg.mapt <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_result_mapt.ndc_fro_2019-10-23_13.33.11.txt",
                       sep="\t", header=T, row.names=1)
deg.mapt <- deg.mapt[deg.mapt$padj <= 0.05,]
deg.mapt <- deg.mapt[rownames(deg.mapt) %in% rownames(mt),]
deg.mapt <- deg.mapt[rownames(mt),]
mt.mapt <- mt
mt.mapt$lfc <- deg.mapt$log2FoldChange
mt.mapt$padj <- deg.mapt$padj
mt.mapt <- mt.mapt[, c(-2, -3, -4, -5)]

# GRN
deg.grn <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_result_grn.ndc_fro_2019-10-23_13.33.11.txt",
                       sep="\t", header=T, row.names=1)
deg.mapt <- deg.mapt[deg.mapt$padj <= 0.05,]
deg.grn <- deg.grn[rownames(deg.grn) %in% rownames(mt),]
deg.grn <- deg.grn[rownames(mt),]
mt.grn <- mt
mt.grn$lfc <- deg.grn$log2FoldChange
mt.grn$padj <- deg.grn$padj
mt.grn <- mt.grn[, c(-2, -3, -4, -5)]



###
# Check for overlap with miRNA targets
###

# MAPT
mir.mapt <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/MAPT_miRNA_target_edge_table.txt", 
                       sep="\t", header=T, stringsAsFactors = F)
mir.mapt <- mir.mapt[mir.mapt$targets %in% mt$name,]
table(mir.mapt$mirna)
# get miRNA fold changes
mapt.mir <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_result_mapt.ndc_frontal_smRNAseq.txt", sep="\t", header=T)

mt.mapt$mir_target <- rep("no", nrow(mt.mapt))
mt.mapt$mir_target[mt.mapt$name %in% mir.mapt$targets] <- "yes"

## GRN
mir.grn <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/GRN_miRNA_target_edge_table.txt", 
                       sep="\t", header=T, stringsAsFactors = F)
mir.grn <- mir.grn[mir.grn$targets %in% mt$name,]
# get miRNA fold changes
grn.mir <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_result_grn.ndc_frontal_smRNAseq.txt", sep="\t", header=T)

mt.grn$mir_target <- rep("no", nrow(mt.grn))
mt.grn$mir_target[mt.grn$name %in% mir.grn$targets] <- "yes"

#==========================================#


###
# Check for alternative splicing
###

## MAPT
as.mapt <- read.table("~/rimod/RNAseq/as_analysis/majiq/mapt_AS_genes_dPSI_0.2.txt", sep="\t", header = T)
as.mapt <- mt.mapt[mt.mapt$ensembl_gene_id %in% as.mapt$x,]

mt.mapt$alt_splice <- rep("no", nrow(mt.mapt))
mt.mapt$alt_splice[mt.mapt$name %in% as.mapt$name] <- "yes"

## GRN
as.grn <- read.table("~/rimod/RNAseq/as_analysis/majiq/mapt_AS_genes_dPSI_0.2.txt", sep="\t", header = T)
as.grn <- mt.grn[mt.grn$ensembl_gene_id %in% as.grn$x,]

mt.grn$alt_splice <- rep("no", nrow(mt.grn))
mt.grn$alt_splice[mt.grn$name %in% as.grn$name] <- "yes"


#===================================#

###
# Check for promotor shifting
###

# MAPT
ps.mapt <- read.table("~/rimod/CAGE/cage_analysis/promotor_shifting/frontal/mapt_promotor_shifting_genes_fro.txt", header=T)
ps.mapt <- mt.mapt[mt.mapt$ensembl_gene_id %in% ps.mapt$x,]

mt.mapt$prom_shift <- rep("no", nrow(mt.mapt))
mt.mapt$prom_shift[mt.mapt$name %in% ps.mapt$name] <- "yes"

# GRN
ps.grn <- read.table("~/rimod/CAGE/cage_analysis/promotor_shifting/frontal/grn_promotor_shifting_genes_fro.txt", header=T)
ps.grn <- mt.grn[mt.grn$ensembl_gene_id %in% ps.grn$x,]

mt.grn$prom_shift <- rep("no", nrow(mt.grn))
mt.grn$prom_shift[mt.grn$name %in% ps.grn$name] <- "yes"


#==========================================#

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
#mapt.met <- mapt.met[mapt.met$adj.P.Val <= 0.01,]
mapt.met.genes <- getGenes(mapt.met)

met_lfc <- c()
for (i in 1:nrow(mt.mapt)) {
  res <- 0
  g <- mt.mapt$name[i]
  if (g %in% mapt.met.genes){
    tmp <- mapt.met[grepl(g, mapt.met$GencodeCompV12_NAME),]
    res <- mean(tmp$logFC)
  }
  met_lfc <- c(met_lfc, res)
}
mt.mapt$meth_lfc <- met_lfc

# GRN
grn.met <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_grn.ndc_quant.txt", sep="\t")
#mapt.met <- mapt.met[mapt.met$adj.P.Val <= 0.01,]
grn.met.genes <- getGenes(grn.met)

met_lfc <- c()
for (i in 1:nrow(mt.grn)) {
  res <- 0
  g <- mt.grn$name[i]
  if (g %in% grn.met.genes){
    tmp <- grn.met[grepl(g, grn.met$GencodeCompV12_NAME),]
    res <- mean(tmp$logFC)
  }
  met_lfc <- c(met_lfc, res)
}
mt.grn$meth_lfc <- met_lfc

















####
# Build iGraph
####
library(igraph)

####
# Build network for FTD-MAPT
####
# get PPI connections
ppi <- read.table("membrane_trafficking_STRING.tsv", sep="\t", header=F)
ppi <- ppi[, c(1,2,15)]

# Create edge list from PPIs
edges <- c()
for (i in 1:nrow(ppi)){
  e <- c(as.character(ppi[i,1]), as.character(ppi[i,2]))
  edges <- c(edges, e)
}

# Add miRNA edges
for (i in 1:nrow(mir.mapt)){
  e <- c(as.character(mir.mapt[i,1]), as.character(mir.mapt[i,2]))
  edges <- c(edges, e)
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
  if (v %in% mapt.mir$X){
    tmp <- mapt.mir[mapt.mir$X == v,]
    fc <- as.numeric(tmp$log2FoldChange)
  }
  else if (v %in% mt.mapt$name){
    tmp <- mt.mapt[mt.mapt$name == v,]
    fc <- as.numeric(tmp$lfc)
  }
  lfcs <- c(lfcs, fc)
}
V(g)$lfc <- lfcs

# Assign Pvalues
pvals <- c()
for (v in vnames) {
  fc = 0
  if (v %in% mapt.mir$X){
    tmp <- mapt.mir[mapt.mir$X == v,]
    p <- as.numeric(tmp$padj)
  }
  else if (v %in% mt.mapt$name){
    tmp <- mt.mapt[mt.mapt$name == v,]
    p <- as.numeric(tmp$padj)
  }
  pvals <- c(pvals, p)
}
V(g)$sig <- pvals <= 0.05


# Assign modules
mods <- c()
for (v in vnames) {
  m = "miRNA"
  if (v %in% mt.mapt$name){
    tmp <- mt.mapt[mt.mapt$name == v,]
    m <- as.character(tmp$module)
  }
  mods <- c(mods, m)
}
V(g)$module <- mods

# Save file
file_name = "MAPT_MembraneTraficking_signaling_network.gml"
write_graph(g, file=file_name, format="gml")
print("Network creation succesfull.")

