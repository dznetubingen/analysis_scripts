#########################
# MAPT HumanBase Module M4-down in-depth analysis
#########################
library(stringr)
library(igraph)
setwd("~/rimod/integrative_analysis/mitochondrion_module_analysis/")


# Load module
mod <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_mapt_filtered_down_modules.txt", sep="\t", header=T)
mod <- mod[mod$CLUSTER_NAME == "M1",]
mod <- as.character(str_split(mod$CLUSTER_GENES, pattern=",", simplify = T))

mod.grn <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_grn_filtered_down_modules.txt", sep="\t", header=T)
mod.grn <- mod.grn[mod.grn$CLUSTER_NAME == "M1",]
mod.grn <- as.character(str_split(mod.grn$CLUSTER_GENES, pattern=",", simplify = T))


# load module connections
con <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_mapt_filtered_down_connections.txt", sep="\t", header=T)
con <- con[con$SOURCE %in% mod,]
con <- con[con$TARGET %in% mod,]
con <- con[, c(1,2,3)]

con.grn <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_grn_filtered_down_connections.txt", sep="\t", header=T)
con.grn <- con.grn[con.grn$SOURCE %in% mod,]
con.grn <- con.grn[con.grn$TARGET %in% mod,]
con.grn <- con.grn[, c(1,2,3)]


# Load miRNA connections
mir <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/MAPT_miRNA_target_edge_table.txt", sep="\t", header=T)
mir <- mir[mir$targets %in% mod,]

mir.grn <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/GRN_miRNA_target_edge_table.txt", sep="\t",  header=T)
mir.grn <- mir.grn[mir.grn$targets %in% mod.grn,]

# Load TF connections
tf <- read.table("~/rimod/integrative_analysis/tf_cage_rnaseq/MAPT_common_TFs_down.txt", sep="\t", header=T, stringsAsFactors = F)
for (i in 1:nrow(tf)) {
  tmp <- intersect(mod, as.character(str_split(tf$Overlapping_Genes[i], pattern=",", simplify = T)))
  if (length(tmp) > 0) {
    tf$Overlapping_Genes[i] <- paste(tmp, collapse=",")
  }
  else {
    tf$Overlapping_Genes[i] <- "none"
  }
}
# Make TF target table
tf.df <- data.frame(TF = "dummy", Target = "dummy")
for (i in 1:nrow(tf)) {
  s <- tf$TF[i]
  targets <- as.character(str_split(tf$Overlapping_Genes[i], pattern=",", simplify = T))
  tmp <- data.frame(TF = rep(s, length(targets)), Target = targets)
  tf.df <- rbind(tf.df, tmp)
}
tf.df <- tf.df[-1,]



# Load methylation data
met <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_mapt.ndc_quant.txt", sep="\t", header = T, stringsAsFactors = F)
met <- met[!met$GencodeBasicV12_NAME == "",]
met.up <- met[met$logFC > 0,]
met.down <- met[met$logFC < 0,]

met.grn <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_grn.ndc_quant.txt", sep="\t", header = T, stringsAsFactors = F)
met.grn <- met.grn[!met.grn$GencodeBasicV12_NAME == "",]
met.grn.up <- met.grn[met.grn$logFC > 0,]
met.grn.down <- met.grn[met.grn$logFC < 0,]

###
# Get genes associated with a certian annotation group from a DMP file
# E.g. get all Genes that have DMPs in their TSS
getDMgenes <- function(m, groupHandle = "TSS"){
  tmp <- c()
  for (i in 1:nrow(m)) {
    group <- as.character(str_split(m$GencodeBasicV12_Group[i], pattern=";", simplify = T))
    genes <- as.character(str_split(m$GencodeBasicV12_NAME[i], pattern=";", simplify = T))
    genes <- genes[grepl(groupHandle, group)]
    tmp <- c(tmp, genes)
  }
  return(tmp)
}

# mapt
met.up.genes <- getDMgenes(met.up)
met.down.genes <- getDMgenes(met.down)
met.up.genes <- intersect(met.up.genes, mod)
met.down.genes <- intersect(met.down.genes, mod)

# grn
met.grn.up.genes <- getDMgenes(met.grn.up)
met.grn.down.genes <- getDMgenes(met.grn.down)
met.grn.up.genes <- intersect(met.grn.up.genes, mod.grn)
met.grn.down.genes <- intersect(met.grn.down.genes, mod.grn)


# Load splicing data
as <- read.table("~/rimod/RNAseq/as_analysis/majiq/mapt_AS_genes_dPSI_0.2_CCF_HGNC.txt", sep="\t")
as <- as.character(as$V1)
as <- intersect(as, mod)

as.grn <- read.table("~/rimod/RNAseq/as_analysis/majiq/grn_AS_genes_dPSI_0.2_CCF_HGNC.txt", sep="\t")
as.grn <- as.character(as.grn$V1)
as.grn <- intersect(as.grn, mod)


#=========================================================#



#####
# Make iGraph network
#####

# gene-gene
edges <- c()
for (i in 1:nrow(con)) {
  e <- c(as.character(con[i,1]), as.character(con[i,2]))
  edges <- c(edges, e)
}

# mirna-gene
for (i in 1:nrow(mir)) {
  e <- c(as.character(mir[i,1]), as.character(mir[i,2]))
  edges <- c(edges, e)
}

# tf-gene
for (i in 1:nrow(tf.df)) {
  e <- c(as.character(tf.df[i,1]), as.character(tf.df[i,2]))
  edges <- c(edges, e)
}

# methylation-gene
for (m in met.up.genes) {
  e <- c("HyperMethylation", m)
  edges <- c(edges, e)
}
for (m in met.down.genes) {
  e <- c("HypoMethylation", m)
  edges <- c(edges, e)
}

# alternative splicing
for (a in as) {
  e <- c("AS", a)
  edges <- c(edges, e)
}

g <- graph(edges)

# assign types to nodes and edges
verts <- V(g)$name
types <- rep("Gene", length(verts))
types[verts %in% mir$mirna] <- "miRNA"
types[verts %in% tf.df$TF] <- "TF"
types[verts %in% c("HyperMethylation", "HypoMethylation")] <- "CpG"
types[verts == "AS"] <- "AlternativeSplicing"
V(g)$type <- types

# save the graph
write_graph(g, file = "MAPT_M1down_network.gml", format = "gml")
write.table(mod, "mapt_module_m1down_genes.txt", row.names=F, col.names=F, quote=F)

#===============================#


#####
# Network analysis in iGraph
#####

###
# MAPT
###
# genes only
edges <- c()
for (i in 1:nrow(con)) {
  e <- c(as.character(con[i,1]), as.character(con[i,2]))
  edges <- c(edges, e)
}

net <- graph(edges, directed = F)
deg <- degree(net, mode="all")
plot(net, vertex.size = deg*2)
deg <- deg[order(deg, decreasing = T)]
deg.mapt <- deg

###
# GRN
###
# genes only
edges <- c()
for (i in 1:nrow(con.grn)) {
  e <- c(as.character(con.grn[i,1]), as.character(con.grn[i,2]))
  edges <- c(edges, e)
}

net <- graph(edges, directed = F)
deg <- degree(net, mode="all")
plot(net, vertex.size = deg*2)
deg <- deg[order(deg, decreasing = T)]
deg.grn <- deg
