###############################
# Figure 3: only TF-module networks
################################
library(stringr)
library(igraph)

setwd("~/rimod/paper/figures/figure3/")

###
# Parameters
###
tf_cutoff = 50
jac_cutoff = 0.01

#===========#

jaccard <- function(a, b){
  jac <- length(intersect(a, b)) / length(union(a, b))
  return(jac)
}


####
# MAPT
####

# TF
tf.up <- read.table("~/rimod/integrative_analysis/tf_cage_rnaseq/MAPT_common_TFs_up.txt", sep="\t", header=T, stringsAsFactors = F)
tf.down <- read.table("~/rimod/integrative_analysis/tf_cage_rnaseq/MAPT_common_TFs_down.txt", sep="\t", header=T, stringsAsFactors = F)
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


# Format Down-TF
tf.down <- tf.down[1:tf_cutoff,]
tf.down <- tf.down[, c(3, 6)]
tf.down <- tf.down[tf.down$TF %in% deg$hgnc_symbol,]

# Make TF edges
jac.df <- data.frame(reg = "dummy", module = "dummy", jaccard = 0)
for (i in 1:nrow(tf.down)) {
  tf <- tf.down[i,1]
  targets <- str_split(tf.down[i,2], pattern=",")[[1]]
  
  for (j in 1:length(down.mods)) {
    module = names(down.mods)[j]
    ovl <- intersect(targets, down.mods[[j]])
    jac <- jaccard(targets, down.mods[[j]])
    
    if (jac > jac_cutoff){
      edges <- c(edges, c(tf, module))
      jac.df <- rbind(jac.df, data.frame(reg = tf, module = module, jaccard = jac))
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

write_graph(g, file = "MAPT_downModules_TF.gml", format = "gml")
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


# Format Down-TF
tf.up <- tf.up[1:tf_cutoff,]
tf.up <- tf.up[, c(3, 6)]
tf.up <- tf.up[tf.up$TF %in% deg$hgnc_symbol,]

# Make TF edges
jac.df <- data.frame(reg = "dummy", module = "dummy", jaccard = 0)
for (i in 1:nrow(tf.up)) {
  tf <- tf.up[i,1]
  targets <- str_split(tf.up[i,2], pattern=",")[[1]]
  
  for (j in 1:length(up.mods)) {
    module = names(up.mods)[j]
    ovl <- intersect(targets, up.mods[[j]])
    jac <- jaccard(targets, up.mods[[j]])
    
    if (jac > jac_cutoff){
      edges <- c(edges, c(tf, module))
      jac.df <- rbind(jac.df, data.frame(reg = tf, module = module, jaccard = jac))
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

write_graph(g, file = "MAPT_upModules_TF.gml", format = "gml")
plot(g)
#=========================================================================#


#####
# FTD-GRN
#####

# TF
tf.up <- read.table("~/rimod/integrative_analysis/tf_cage_rnaseq/GRN_common_TFs_up.txt", sep="\t", header=T, stringsAsFactors = F)
tf.down <- read.table("~/rimod/integrative_analysis/tf_cage_rnaseq/GRN_common_TFs_down.txt", sep="\t", header=T, stringsAsFactors = F)
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


# Format Down-TF
tf.down <- tf.down[1:tf_cutoff,]
tf.down <- tf.down[, c(3, 6)]
tf.down <- tf.down[tf.down$TF %in% deg$hgnc_symbol,]

# Make TF edges
jac.df <- data.frame(reg = "dummy", module = "dummy", jaccard = 0)
for (i in 1:nrow(tf.down)) {
  tf <- tf.down[i,1]
  targets <- str_split(tf.down[i,2], pattern=",")[[1]]
  
  for (j in 1:length(down.mods)) {
    module = names(down.mods)[j]
    ovl <- intersect(targets, down.mods[[j]])
    jac <- jaccard(targets, down.mods[[j]])
    
    if (jac > jac_cutoff){
      edges <- c(edges, c(tf, module))
      jac.df <- rbind(jac.df, data.frame(reg = tf, module = module, jaccard = jac))
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

write_graph(g, file = "GRN_downModules_TF.gml", format = "gml")
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


# Format Down-TF
tf.up <- tf.up[1:tf_cutoff,]
tf.up <- tf.up[, c(3, 6)]
tf.up <- tf.up[tf.up$TF %in% deg$hgnc_symbol,]

# Make TF edges
jac.df <- data.frame(reg = "dummy", module = "dummy", jaccard = 0)
for (i in 1:nrow(tf.up)) {
  tf <- tf.up[i,1]
  targets <- str_split(tf.up[i,2], pattern=",")[[1]]
  
  for (j in 1:length(up.mods)) {
    module = names(up.mods)[j]
    ovl <- intersect(targets, up.mods[[j]])
    jac <- jaccard(targets, up.mods[[j]])
    
    if (jac > jac_cutoff){
      edges <- c(edges, c(tf, module))
      jac.df <- rbind(jac.df, data.frame(reg = tf, module = module, jaccard = jac))
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

write_graph(g, file = "GRN_upModules_TF.gml", format = "gml")
plot(g)
#=========================================================================#