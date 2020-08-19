#############
# Figure 6
# Integration of miRNAs and modules 
# to generate plots similar to the regulon analysis plots
#############
library(stringr)
library(igraph)

setwd("~/rimod/paper_v2/figures/figure6//")

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
# GRN M3 down
####

mod <- read.table("string_interactions_GRN_M3down.tsv", sep="\t")
mod.genes <- union(mod$V1, mod$V2)
# miRNA
mir <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/GRN_miRNA_target_edge_table.txt", sep="\t", header=T, stringsAsFactors = F)
mir <- mir[mir$targets %in% mod.genes,]

# make them edges
edges <- c()
for (i in 1:nrow(mod)) {
  e <- c(as.character(mod[i,]$V1), as.character(mod[i,]$V2))
  edges <- c(edges, e)
}
for (i in 1:nrow(mir)) {
  e <- c(as.character(mir[i,]$mirna), as.character(mir[i,]$targets))
  edges <- c(edges, e)
}

# Make Graph

g <- graph(edges=edges)
write_graph(g, file = "GRN_M3down_miRNA_targets.gml", format = "gml")
plot(g)


####
# MAPT
####
# miRNA
mir <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/MAPT_miRNA_target_edge_table.txt", sep="\t", header=T, stringsAsFactors = F)
mir <- mir[mir$mirna %in% c("hsa-miR-150-5p", "hsa-miR-193a-3p", "hsa-miR-142-3p", "hsa-miR-148a-3p", "hsa-miR-363-3p"),]

write.table(mir$targets, "MAPT_mir_targetes.txt", quote=F, row.names = F, col.names = F)

# load STRING network
mod <- read.table("string_interactions_MAPT_common_mir_targets.tsv", sep="\t")
mod.genes <- union(mod$V1, mod$V2)
mir <- mir[mir$targets %in% mod.genes,]


# make them edges
edges <- c()
for (i in 1:nrow(mod)) {
  e <- c(as.character(mod[i,]$V1), as.character(mod[i,]$V2))
  edges <- c(edges, e)
}
for (i in 1:nrow(mir)) {
  e <- c(as.character(mir[i,]$mirna), as.character(mir[i,]$targets))
  edges <- c(edges, e)
}

# Make Graph

g <- graph(edges=edges)
write_graph(g, file = "MAPT_common_miRNA_targets.gml", format = "gml")
plot(g)


