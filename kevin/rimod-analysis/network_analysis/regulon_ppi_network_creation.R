########################################
## TF-regulon PPI network integration ##
########################################
library(igraph)
library(biomaRt)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
reg_file <- as.character(args[1])
deg_file <- as.character(args[2])
pval_cut <- as.numeric(args[3])
lfc_cut <- as.numeric(args[4])
outname <- as.character(args[5])


## for testing
setwd("~/rimod/integrative_analysis/tf_regulon_ppi_integration/")
reg_file <- "~/rimod/CAGE/cage_analysis/regulon_analysis_degbg/c9orf72/c9orf72_regulons.rds"
deg_file <- "~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2018-04-26_genewise/deseq_result_c9.ndc_2018-04-26_14.22.04.txt"
pval_cut <- 0.05
lfc_cut <- 0.5
outname = "c9orf72"


#================= Load Data ================#
# Read TF-gene mappings
tfgene <- read.table("~/rimod/CAGE/cage_analysis/regulon_activity_enrichment_analysis_300418/TF_gene_targets_030518.txt", sep="\t", header = T)
# Load IntAct PPIs
ppi <- read.table("~/resources/ppi/IntAct/intact/hs.intact.sub.ensembl.txt", sep="\t", header = T, fill = T)
ppi <- ppi[,7:ncol(ppi)]

# Load regulons
regs <- readRDS(reg_file)

# Load DEGs
deg <- read.table(deg_file, sep="\t", header = T, row.names=1)
deg <- deg[deg$padj <= pval_cut,]
deg <- deg[abs(deg$log2FoldChange) >= lfc_cut,]

# Get ensembl gene IDs for TFs
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
up.tfs <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "hgnc_symbol", values = names(regs), mart = ensembl)
# Only consider sig. up-regulated TF-regulons
deg.up <- deg[deg$log2FoldChange >= lfc_cut,]
up.tfs <- up.tfs[up.tfs$ensembl_gene_id %in% rownames(deg.up),]
regs <- regs[names(regs) %in% up.tfs$hgnc_symbol]

#======== Network Creation ================#

#### Create PPI network
gran <- "ENSG00000030582" ## add GRN manually
mapt <- "ENSG00000186868" ## add Mapt manually
c9orf <- "ENSG00000147894"
degs <- rownames(deg)
degs <- c(degs, gran, mapt, c9orf)
ppi <- ppi[ppi$ensgA %in% degs,]
ppi <- ppi[ppi$ensgB %in% degs,]
ppi <- ppi[!duplicated(ppi),]
edges_ppi <- c()
for (i in 1:nrow(ppi)) {
  e <- c(as.character(ppi[i,1]), as.character(ppi[i,2]))
  edges_ppi <- c(edges_ppi, e)
}

#### Create TF-TF network
tf_names <- up.tfs$ensembl_gene_id
up.net<- tfgene[tfgene$gene_id %in% up.tfs$ensembl_gene_id,]
up.net <- up.net[up.net$tf_gene_id %in% rownames(deg[deg$log2FoldChange > 0,]),]
up.net <- up.net[!up.net$tf_symbol == "",]
up.net <- up.net[!up.net$gene_symbol == "",]
edges_tf <- c()
for (i in 1:nrow(up.net)){
  e <- c(as.character(up.net[i,2]), as.character(up.net[i,4]))
  edges_tf <- c(edges_tf,e)
}

#### Create TF-gene network
edges_tfgene <- c()
# Up-regulons
for (i in 1:length(regs)) {
  tf <- up.tfs[up.tfs$hgnc_symbol == names(regs)[i],]$ensembl_gene_id
  if (length(tf) > 1) tf = tf[1]
  targets <- as.character(unlist(regs[i]))
  for (tgt in targets) {
    e <- c(as.character(tf), as.character(tgt))
    if(length(e) != 2){
    }
    edges_tfgene <- c(edges_tfgene, e)
  }
}

### Merge networks
edges <- c(edges_ppi, edges_tf, edges_tfgene)
g <- graph(edges=edges)

genes <- V(g)$name
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = genes, mart = ensembl)
bm <- bm[match(genes, bm$ensembl_gene_id),]
bm$hgnc_symbol[bm$hgnc_symbol == ""] <- "no_hgnc"
V(g)$symbol <- as.character(bm$hgnc_symbol)

### Assign node and edges properties
V(g)$type <- rep("gene", length(V(g)$name))
V(g)$type[V(g)$name %in% levels(factor(tfgene$tf_gene_id))] <- "TF"
V(g)$type[V(g)$name %in% tf_names] <- "TF_sig"
V(g)$state <- rep("uknown", length(V(g)$name))
deg.down <- as.character(rownames(deg[deg$log2FoldChange < 0,]))
V(g)$state[V(g)$name %in% deg.down] <- "down"
V(g)$state[V(g)$name %in%  as.character(rownames(deg[deg$log2FoldChange > 0,]))] <- "up"
ppi_edge <- rep("ppi", length(edges_ppi)/2)
tf_edge <- rep("tf_tf", length(edges_tf)/2)
tfg_edge <- rep("tf_gene", length(edges_tfgene)/2)
E(g)$type <- c(ppi_edge, tf_edge, tfg_edge)
write_graph(g, file=paste(outname, "_tf_ppi_gene_network.gml", sep=""), format="gml")



