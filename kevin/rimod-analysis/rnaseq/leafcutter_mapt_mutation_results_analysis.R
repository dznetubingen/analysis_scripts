##################################
# Analysis of LeafCutter results
# for MAPT Mutations
##################################
library(stringr)
setwd("~/rimod/RNAseq/as_analysis/leafcutter/mapt_mutation_analysis/")

#===========================================#
# Utility functions
getClusters <- function(x){
  res <- as.character(str_split(x, pattern=":", simplify = TRUE)[,2])
  return(res)
}

getESClusters <- function(x){
  res <- as.character(str_split(x, pattern=":", simplify = TRUE)[,4])
  return(res)
}

unlistGenes <- function(x){
  genes <- c()
  for (elem in x){
    elems <- strsplit(elem, split=",")
    for (e in elems){
      genes <- c(genes, e)
    }
  }
  genes <- genes[!duplicated(genes)]
  return(genes)
}
#============================================#

#====================================#
# Cutoffs
q_cutoff = 0.05
dpsi_cutoff = 0.1
print(paste("Cutoff:", q_cutoff))
print(paste("dPSI:", dpsi_cutoff))
#=====================================#

# P301L
p301l <- read.table("p301l_results/leafcutter_ds_cluster_significance.txt", sep="\t", header=T, stringsAsFactors = F)
p301l_es <- read.table("p301l_results/leafcutter_ds_effect_sizes.txt", sep="\t", header=T, stringsAsFactors = F)
p301l <- p301l[!is.na(p301l$p),]
p301l <- p301l[p301l$p.adjust <= q_cutoff,]
p301l$clusterID <- getClusters(p301l$cluster)
p301l_es$clusterID <- getESClusters(p301l_es$intron)
p301l_es <- p301l_es[p301l_es$clusterID %in% p301l$clusterID,]
# filter for dPSI
p301l_es <- p301l_es[abs(p301l_es$deltapsi) >= dpsi_cutoff,]
p301l <- p301l[p301l$clusterID %in% p301l_es$clusterID,]
p301l_genes <- unlistGenes(p301l$genes)
write.table(p301l_genes, paste0("p301l_results/p301l_ds_genes_Q", q_cutoff, "_dPSI", dpsi_cutoff, ".txt"), sep="\t", quote=F, row.names = F)

print(paste("p301l:", dim(p301l)[1]))

# G272V
g272v <- read.table("g272v_results/leafcutter_ds_cluster_significance.txt", sep="\t", header=T, stringsAsFactors = F)
g272v_es <- read.table("g272v_results//leafcutter_ds_effect_sizes.txt", sep="\t", header=T, stringsAsFactors = F)
g272v <- g272v[!is.na(g272v$p),]
g272v <- g272v[g272v$p.adjust <= q_cutoff,]
g272v$clusterID <- getClusters(g272v$cluster)
g272v_es$clusterID <- getESClusters(g272v_es$intron)
g272v_es <- g272v_es[g272v_es$clusterID %in% g272v$clusterID,]
# filter for dPSI
g272v_es <- g272v_es[abs(g272v_es$deltapsi) >= dpsi_cutoff,]
g272v <- g272v[g272v$clusterID %in% g272v_es$clusterID,]
g272v_genes <- unlistGenes(g272v$genes)
write.table(g272v_genes, paste0("g272v_results/g272v_ds_genes_Q", q_cutoff, "_dPSI", dpsi_cutoff, ".txt"), sep="\t", quote=F, row.names = F)

print(paste("G272V:", dim(g272v)[1]))

# r406v
r406v <- read.table("r406v_results/leafcutter_ds_cluster_significance.txt", sep="\t", header=T, stringsAsFactors = F)
r406v_es <- read.table("r406v_results/leafcutter_ds_effect_sizes.txt", sep="\t", header=T, stringsAsFactors = F)
r406v <- r406v[!is.na(r406v$p),]
r406v <- r406v[r406v$p.adjust <= q_cutoff,]
r406v$clusterID <- getClusters(r406v$cluster)
r406v_es$clusterID <- getESClusters(r406v_es$intron)
r406v_es <- r406v_es[r406v_es$clusterID %in% r406v$clusterID,]
# filter for dPSI
r406v_es <- r406v_es[abs(r406v_es$deltapsi) >= dpsi_cutoff,]
r406v <- r406v[r406v$clusterID %in% r406v_es$clusterID,]
r406v_genes <- unlistGenes(r406v$genes)
write.table(r406v_genes, paste0("r406v_results/r406v_ds_genes_Q", q_cutoff, "_dPSI", dpsi_cutoff, ".txt"), sep="\t", quote=F, row.names = F)

print(paste("r406v:", dim(r406v)[1]))

# l315r
l315r <- read.table("l315r_results/leafcutter_ds_cluster_significance.txt", sep="\t", header=T, stringsAsFactors = F)
l315r_es <- read.table("l315r_results/leafcutter_ds_effect_sizes.txt", sep="\t", header=T, stringsAsFactors = F)
l315r <- l315r[!is.na(l315r$p),]
l315r <- l315r[l315r$p.adjust <= q_cutoff,]
l315r$clusterID <- getClusters(l315r$cluster)
l315r_es$clusterID <- getESClusters(l315r_es$intron)
l315r_es <- l315r_es[l315r_es$clusterID %in% l315r$clusterID,]
# filter for dPSI
l315r_es <- l315r_es[abs(l315r_es$deltapsi) >= dpsi_cutoff,]
l315r <- l315r[l315r$clusterID %in% l315r_es$clusterID,]
l315r_genes <- unlistGenes(l315r$genes)
write.table(l315r_genes, paste0("l315r_results/l315r_ds_genes_Q", q_cutoff, "_dPSI", dpsi_cutoff, ".txt"), sep="\t", quote=F, row.names = F)

print(paste("l315r:", dim(l315r)[1]))


