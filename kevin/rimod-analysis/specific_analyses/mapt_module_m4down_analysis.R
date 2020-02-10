#########################
# MAPT HumanBase Module M4-down in-depth analysis
#########################
library(stringr)
library(igraph)
setwd("~/rimod/integrative_analysis/mapt_m4_down/")


# Load module
mod <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_mapt_filtered_down_modules.txt", sep="\t", header=T)
mod <- mod[mod$CLUSTER_NAME == "M4",]
mod <- as.character(str_split(mod$CLUSTER_GENES, pattern=",", simplify = T))

# load module connections
con <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_mapt_filtered_down_connections.txt", sep="\t", header=T)
con <- con[con$SOURCE %in% mod,]
con <- con[con$TARGET %in% mod,]
con <- con[, c(1,2,3)]


# Load miRNA connections
mir <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/MAPT_miRNA_target_edge_table.txt", sep="\t", header=T)
mir <- mir[mir$targets %in% mod,]

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
#met <- met[abs(met$logFC) > 0.6,]
met.up <- met[met$logFC > 0,]
met.down <- met[met$logFC < 0,]

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

met.up.genes <- getDMgenes(met.up)
met.down.genes <- getDMgenes(met.down)
met.up.genes <- intersect(met.up.genes, mod)
met.down.genes <- intersect(met.down.genes, mod)

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
#for (i in 1:nrow(tf.df)) {
#  e <- c(as.character(tf.df[i,1]), as.character(tf.df[i,2]))
#  edges <- c(edges, e)
#}

# methylation-gene
for (m in met.up.genes) {
  e <- c("HyperMethylation", m)
  edges <- c(edges, e)
}
for (m in met.down.genes) {
  e <- c("HypoMethylation", m)
  edges <- c(edges, e)
}

g <- graph(edges)

# assign types to nodes and edges
verts <- V(g)$name
types <- rep("Gene", length(verts))
types[verts %in% mir$mirna] <- "miRNA"
types[verts %in% tf.df$TF] <- "TF"
types[verts %in% c("HyperMethylation", "HypoMethylation")] <- "CpG"
V(g)$type <- types

# save the graph
write_graph(g, file = "MAPT_M4down_network.gml", format = "gml")
