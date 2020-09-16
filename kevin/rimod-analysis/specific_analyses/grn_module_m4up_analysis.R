#########################
# GRN HumanBase Module M4-up in-depth analysis
#########################
library(stringr)
library(igraph)
setwd("~/dzne/rimod_analysis/grn_m4_up_integration//")


# Load module
mod <- read.table("~/dzne/rimod_package/analysis/human_base/rnaseq_grn_filtered_up_modules.txt", sep="\t", header=T)
mod <- mod[mod$CLUSTER_NAME == "M4",]
mod <- as.character(str_split(mod$CLUSTER_GENES, pattern=",", simplify = T))

# load module connections
con <- read.table("~/dzne/rimod_package/analysis/human_base/rnaseq_grn_filtered_up_connections.txt", sep="\t", header=T)
con <- con[con$SOURCE %in% mod,]
con <- con[con$TARGET %in% mod,]
con <- con[, c(1,2,3)]


# Load miRNA connections
mir <- read.table("~/dzne/rimod_package/smRNAseq/analysis/target_mrna_correlation_analysis_0819/GRN_miRNA_target_edge_table.txt", sep="\t", header=T)
mir <- mir[mir$targets %in% mod,]

# Load TF connections
tf <- read.table("~/dzne/rimod_package/integrative_analysis/tf_cage_rnaseq/GRN_common_TFs_up.txt", sep="\t", header=T, stringsAsFactors = F)
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
met <- read.table("~/dzne/rimod_package/frontal_methylation_0818/DMPs_grn.ndc_quant.txt", sep="\t", header = T, stringsAsFactors = F)
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

# Load splicing data
as <- read.table("~/dzne/rimod_package/as_analysis/grn_AS_genes_dPSI_0.2_CCF_HGNC.txt", sep="\t")

as <- as.character(as$V1)
as <- intersect(as, mod)
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
write_graph(g, file = "GRN_M4up_network.gml", format = "gml")
write.table(mod, "grn_module_m4up_genes.txt", row.names=F, col.names=F, quote=F)

#===============================#


#####
# Network analysis in iGraph
#####

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


hs = hub_score(net, weights=NA)$vector
plot(net, vertex.size = hs*100)


ceb <- cluster_fast_greedy(net)
plot(ceb, net)


#
library(ggplot2)
n_genes = 20

# degree analysis
deg <- degree(g)
deg <- deg[types == "Gene"]
deg <- deg[order(deg, decreasing = T)]
df <- as.data.frame(deg)
df$Gene <- rownames(df)
colnames(df) <- c("Degree", "Gene")
df <- df[1:n_genes,]
df$Gene <- factor(as.character(df$Gene), levels = as.character(df$Gene))

p <- ggplot(df, aes(x=Gene, y=Degree)) + 
  geom_bar(stat="identity") + 
  theme_minimal() + 
  coord_flip() + ggtitle("Degree")
p

ggsave("degree_plot.png", width=2, height=3)



# betweennes
bet <- betweenness(g)
bet <- bet[types == "Gene"]
bet <- bet[order(bet, decreasing = T)]
df <- as.data.frame(bet)
df$Gene <- rownames(df)
colnames(df) <- c("Betweenness", "Gene")
df <- df[1:n_genes,]
df$Gene <- factor(as.character(df$Gene), levels = as.character(df$Gene))

p <- ggplot(df, aes(x=Gene, y=Betweenness)) + 
  geom_bar(stat="identity") + 
  theme_minimal() + 
  coord_flip() + ggtitle("Betweenness")
p
ggsave("betweenness_plot.png", width=2, height=3)


# for string
genes <- V(g)$name
write.table(genes, "module_genes.txt", quote=F, col.names = F, row.names = F)
