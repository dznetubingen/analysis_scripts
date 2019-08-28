#########################################################
# Generate igraph network for export to Cytoscape
# Input: miRNA-target edges, TF-target edges
######################################################
library(igraph)
library(stringr)
library(biomaRt)

setwd("~/rimod/integrative_analysis/gene_regulatory_network/")

###
# MAPT up-regulated
###
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# load DEG list
deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-08-12_07.58.35/deseq_result_mapt.ndc_fro_2019-08-12_07.58.35.txt",
                  sep="\t", header=T)
deg <- deg[deg$padj <= 0.05,]
deg <- deg[deg$log2FoldChange > 0,]

# get HGNC symbols
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=as.character(deg$X), mart=ensembl)
deg <- merge(deg, bm, by.x = "X", by.y="ensembl_gene_id")

# load miRNA DEGs
mirs <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_result_mapt.ndc_frontal_smRNAseq.txt", sep="\t", header=T)
mirs <- mirs[mirs$padj <= 0.05,]
mirs <- mirs[mirs$log2FoldChange < -0.8,]


# load miRNA-target interactions
mir.edge <- read.table("resources/MAPT_miRNA_target_edge_table.txt", sep="\t", header=T)
mir.edge$action <- rep(-1, nrow(mir.edge))
mir.edge <- mir.edge[mir.edge$mirna %in% mirs$X,]


# load TFs of interest and generate tf edge list
tf <- read.table("~/rimod/CAGE/cage_analysis/chea_tf_regulators/frontal/Integrated_meanRank_MAPT_up.tsv", sep="\t", header=T)
tf <- tf[1:30,]
tfs <- c()
tf_targets <- c()
for (i in 1:nrow(tf)){
  n <- as.character(tf[i,3])
  n_targets <- str_split(tf[i,6], pattern=",")[[1]]
  n <- rep(n, length(n_targets))
  
  tfs <- c(tfs, n)
  tf_targets <- c(tf_targets, n_targets)
}
tf.edge <- data.frame(tf = tfs, targets = tf_targets)

# subset with DE genes
tf.edge <- tf.edge[tf.edge$targets %in% deg$hgnc_symbol,]


##
# Load PPIs
##
score_threshold <- 500
ppi <- read.table("resources/9606.protein.links.v11.0.txt", sep=" ", header = T)
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_peptide_id"), filters="hgnc_symbol", values=deg$hgnc_symbol, mart=ensembl)
ppi$protein1 <- gsub("9606.", "", ppi$protein1)
ppi$protein2 <- gsub("9606.", "", ppi$protein2)
ppi <- ppi[ppi$protein1 %in% bm$ensembl_peptide_id,]
ppi <- ppi[ppi$protein2 %in% bm$ensembl_peptide_id,]
ppi <- ppi[ppi$combined_score >= score_threshold,]
# convert to HGNC symbols
ppi <- merge(ppi, bm, by.x = "protein1", by.y = "ensembl_peptide_id")
ppi$gene1 <- ppi$hgnc_symbol
ppi <- merge(ppi, bm, by.x= "protein2", by.y="ensembl_peptide_id")
ppi$gene2 <- ppi$hgnc_symbol.y
ppi <- ppi[,c("gene1", "gene2", "combined_score")]


#=====================================================================#

#####
# Create iGraph network
######
edges <- c()
# generate mir edges
for (i in 1:nrow(mir.edge)){
  e <- c(as.character(mir.edge[i,1]), as.character(mir.edge[i,2]))
  edges <- c(edges, e)
}

# generate tf-gene edges
for (i in 1:nrow(tf.edge)){
  e <- c(as.character(tf.edge[i,1]), as.character(tf.edge[i,2]))
  edges <- c(edges, e)
}

# generate PPI edges

edges <- make.names(edges)
for (i in 1:nrow(ppi)){
  e <- c(as.character(ppi[i,1]), as.character(ppi[i,2]))
  edges <- c(edges, e)
}


g <- graph(edges=edges, directed=T)
l <- layout_with_fr(g)
plot(g, vertex.color = "gold", vertex.size = 5, edge.arrow.size = .5, layout = l, 
     vertex.label.cex = 0.6, vertex.label.color = "black")

###
# Collect information on vertices and add to graph
###

# Assign log fold changes
vnames <- V(g)$name
net.lfcs <- c()
for (i in 1:length(vnames)){
  v <- vnames[i]
  lfc <- 0
  if (v %in% deg$hgnc_symbol){
    lfc <- as.numeric(deg[deg$hgnc_symbol == v,]$log2FoldChange)
  }
  else if (v %in% make.names(mirs$X)){
    lfc <- as.numeric(mirs[make.names(mirs$X) == v,]$log2FoldChange)
  }
  net.lfcs <- c(net.lfcs, lfc)
}
V(g)$lfc <- net.lfcs

# Assign types
typ <- c()
for (v in vnames){
  tp <- ""
  if (grepl("hsa.", v)){
    tp <- "miRNA"
  }
  else if (v %in% tf.edge$tf){
    tp <- "TF"
  }
  else {
    tp <- "Gene"
  }
  typ <- c(typ, tp)
}
V(g)$type <- typ



write.graph(g, file = "MAPT_UPdegs_miRNA_TF_network.gml", format="gml")



