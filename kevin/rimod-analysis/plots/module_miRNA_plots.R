#############
# Figure 4
# Integration of miRNAs and modules 
# to generate plots similar to the regulon analysis plots
#############
library(stringr)
library(igraph)

setwd("~/rimod/paper/figures/figure4/")

###
# Parameters
###
jac_cutoff = 0.02

#===========#

jaccard <- function(a, b){
  jac <- length(intersect(a, b)) / length(union(a, b))
  return(jac)
}


####
# MAPT
####
# miRNA
mir <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/MAPT_miRNA_target_edge_table.txt", sep="\t", header=T, stringsAsFactors = F)
mir.deg <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_result_mapt.ndc_frontal_smRNAseq.txt", sep="\t", header=T, stringsAsFactors = F)

# DEG
deg <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/MAPT_cell_composition_filterd_DEGs.txt", sep="\t", header=T, stringsAsFactors = F)

# Modules
mod.up <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_mapt_filtered_up_modules.txt", sep="\t", header=T)
mod.down <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_mapt_filtered_down_modules.txt", sep="\t", header=T)




####
# Down-modules
down.mods <- list()
for (m in as.character(mod.down$CLUSTER_NAME)) {
  tmp <- mod.down[mod.down$CLUSTER_NAME == m,]
  genes <- as.character(tmp$CLUSTER_GENES)
  genes <- str_split(genes, pattern=",")[[1]]
  down.mods$tmp <- genes
  names(down.mods)[length(down.mods)] <- m
}

# Graph edges
edges <- c()

# Make TF edges
jac.df <- data.frame(reg = "dummy", module = "dummy", jaccard = 0)


# Get DE miRNAs
mirs <- as.character(mir$mirna)
mirs <- mirs[!duplicated(mirs)]
for (m in mirs) {
  tmp <- mir[mir$mirna == m,]
  targets <- as.character(tmp$targets)
  
  for (j in 1:length(down.mods)){
    module = names(down.mods)[j]
    ovl <- intersect(targets, down.mods[[j]])
    jac <- jaccard(targets, down.mods[[j]])
    
    
    if (jac > jac_cutoff){
      edges <- c(edges, c(m, module))
      jac.df <- rbind(jac.df, data.frame(reg = m, module = module, jaccard = jac))
    }
  }
}

g <- graph(edges=edges)

# Add sizes
names <- V(g)$name
sizes <- rep(1, length(names))
for (i in 1:length(names)) {
  n <- names[i]
  if (n %in% mod.down$CLUSTER_NAME){
    sizes[i] <- length(down.mods[names(down.mods) == n][[1]])
  }
}
V(g)$size <- sizes

# Add jaccard scores
jac.df <- jac.df[-1,]
E(g)$jaccard <- jac.df$jaccard

# Make Graph

write_graph(g, file = "MAPT_downModules_miRNA.gml", format = "gml")
plot(g)
#=================================================#

####
# Up-modules
## Up-modules
up.mods <- list()
for (m in as.character(mod.up$CLUSTER_NAME)) {
  tmp <- mod.up[mod.up$CLUSTER_NAME == m,]
  genes <- as.character(tmp$CLUSTER_GENES)
  genes <- str_split(genes, pattern=",")[[1]]
  up.mods$tmp <- genes
  names(up.mods)[length(up.mods)] <- m
}

# Graph edges
edges <- c()

# Get DE miRNAs
mirs <- as.character(mir$mirna)
mirs <- mirs[!duplicated(mirs)]
for (m in mirs) {
  tmp <- mir[mir$mirna == m,]
  targets <- as.character(tmp$targets)
  
  for (j in 1:length(up.mods)){
    module = names(up.mods)[j]
    ovl <- intersect(targets, up.mods[[j]])
    jac <- jaccard(targets, up.mods[[j]])
    
    
    if (jac > jac_cutoff){
      edges <- c(edges, c(m, module))
      jac.df <- rbind(jac.df, data.frame(reg = m, module = module, jaccard = jac))
    }
  }
}

g <- graph(edges=edges)

# Add sizes
names <- V(g)$name
sizes <- rep(1, length(names))
for (i in 1:length(names)) {
  n <- names[i]
  if (n %in% mod.up$CLUSTER_NAME){
    sizes[i] <- length(up.mods[names(up.mods) == n][[1]])
  }
}
V(g)$size <- sizes

# Add jaccard scores
jac.df <- jac.df[-1,]
E(g)$jaccard <- jac.df$jaccard

# Make Graph

write_graph(g, file = "MAPT_upModules_miRNA.gml", format = "gml")
plot(g)
#=========================================================================#


#####
# FTD-GRN
#####

# miRNA
mir <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/GRN_miRNA_target_edge_table.txt", sep="\t", header=T, stringsAsFactors = F)
mir.deg <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_result_grn.ndc_frontal_smRNAseq.txt", sep="\t", header=T, stringsAsFactors = F)

# DEG
deg <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/GRN_cell_composition_filterd_DEGs.txt", sep="\t", header=T, stringsAsFactors = F)

# Modules
mod.up <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_grn_filtered_up_modules.txt", sep="\t", header=T)
mod.down <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_grn_filtered_down_modules.txt", sep="\t", header=T)




####
# Down-modules
down.mods <- list()
for (m in as.character(mod.down$CLUSTER_NAME)) {
  tmp <- mod.down[mod.down$CLUSTER_NAME == m,]
  genes <- as.character(tmp$CLUSTER_GENES)
  genes <- str_split(genes, pattern=",")[[1]]
  down.mods$tmp <- genes
  names(down.mods)[length(down.mods)] <- m
}

# Graph edges
edges <- c()

# Get DE miRNAs
mirs <- as.character(mir$mirna)
mirs <- mirs[!duplicated(mirs)]
for (m in mirs) {
  tmp <- mir[mir$mirna == m,]
  targets <- as.character(tmp$targets)
  
  for (j in 1:length(down.mods)){
    module = names(down.mods)[j]
    ovl <- intersect(targets, down.mods[[j]])
    jac <- jaccard(targets, down.mods[[j]])
    
    
    if (jac > jac_cutoff){
      edges <- c(edges, c(m, module))
      jac.df <- rbind(jac.df, data.frame(reg = m, module = module, jaccard = jac))
    }
  }
}

g <- graph(edges=edges)

# Add sizes
names <- V(g)$name
sizes <- rep(1, length(names))
for (i in 1:length(names)) {
  n <- names[i]
  if (n %in% mod.down$CLUSTER_NAME){
    sizes[i] <- length(down.mods[names(down.mods) == n][[1]])
  }
}
V(g)$size <- sizes

# Add jaccard scores
jac.df <- jac.df[-1,]
E(g)$jaccard <- jac.df$jaccard

# Make Graph

write_graph(g, file = "GRN_downModules_miRNA.gml", format = "gml")
plot(g)
#=================================================#

####
# Up-modules
## Up-modules
up.mods <- list()
for (m in as.character(mod.up$CLUSTER_NAME)) {
  tmp <- mod.up[mod.up$CLUSTER_NAME == m,]
  genes <- as.character(tmp$CLUSTER_GENES)
  genes <- str_split(genes, pattern=",")[[1]]
  up.mods$tmp <- genes
  names(up.mods)[length(up.mods)] <- m
}

# Graph edges
edges <- c()



# Get DE miRNAs
mirs <- as.character(mir$mirna)
mirs <- mirs[!duplicated(mirs)]
for (m in mirs) {
  tmp <- mir[mir$mirna == m,]
  targets <- as.character(tmp$targets)
  
  for (j in 1:length(up.mods)){
    module = names(up.mods)[j]
    ovl <- intersect(targets, up.mods[[j]])
    jac <- jaccard(targets, up.mods[[j]])
    
    
    if (jac > jac_cutoff){
      edges <- c(edges, c(m, module))
      jac.df <- rbind(jac.df, data.frame(reg = m, module = module, jaccard = jac))
    }
  }
}

g <- graph(edges=edges)

# Add sizes
names <- V(g)$name
sizes <- rep(1, length(names))
for (i in 1:length(names)) {
  n <- names[i]
  if (n %in% mod.up$CLUSTER_NAME){
    sizes[i] <- length(up.mods[names(up.mods) == n][[1]])
  }
}
V(g)$size <- sizes

# Add jaccard scores
jac.df <- jac.df[-1,]
E(g)$jaccard <- jac.df$jaccard

# Make Graph

write_graph(g, file = "GRN_upModules_miRNA.gml", format = "gml")
plot(g)
#=========================================================================#