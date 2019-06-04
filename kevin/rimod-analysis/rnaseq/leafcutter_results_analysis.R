##################################
# Analysis of LeafCutter results
##################################
library(stringr)
setwd("~/rimod/RNAseq/as_analysis/leafcutter/")

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


# Cutoffs
q_cutoff = 0.05
dpsi_cutoff = 0.1
print(paste("Cutoff:", q_cutoff))
print(paste("dPSI:", dpsi_cutoff))
# MAPT
mapt <- read.table("mapt_results/leafcutter_ds_cluster_significance.txt", sep="\t", header=T, stringsAsFactors = F)
mapt_es <- read.table("mapt_results/leafcutter_ds_effect_sizes.txt", sep="\t", header=T, stringsAsFactors = F)
mapt <- mapt[!is.na(mapt$p),]
mapt <- mapt[mapt$p.adjust <= q_cutoff,]
mapt$clusterID <- getClusters(mapt$cluster)
mapt_es$clusterID <- getESClusters(mapt_es$intron)
mapt_es <- mapt_es[mapt_es$clusterID %in% mapt$clusterID,]
# filter for dPSI
mapt_es <- mapt_es[abs(mapt_es$deltapsi) >= dpsi_cutoff,]
mapt <- mapt[mapt$clusterID %in% mapt_es$clusterID,]
mapt_genes <- unlistGenes(mapt$genes)
write.table(mapt_genes, paste0("mapt_results/mapt_ds_genes_Q", q_cutoff, "_dPSI", dpsi_cutoff, ".txt"), sep="\t", quote=F, row.names = F)

print(paste("MAPT:", dim(mapt)[1]))

# GRN
grn <- read.table("grn_results/leafcutter_ds_cluster_significance.txt", sep="\t", header=T, stringsAsFactors = F)
grn_es <- read.table("grn_results/leafcutter_ds_effect_sizes.txt", sep="\t", header=T, stringsAsFactors = F)
grn <- grn[!is.na(grn$p),]
grn <- grn[grn$p.adjust <= q_cutoff,]
grn$clusterID <- getClusters(grn$cluster)
grn_es$clusterID <- getESClusters(grn_es$intron)
grn_es <- grn_es[grn_es$clusterID %in% grn$clusterID,]
# filter for dPSI
grn_es <- grn_es[abs(grn_es$deltapsi) >= dpsi_cutoff,]
grn <- grn[grn$clusterID %in% grn_es$clusterID,]
grn_genes <- unlistGenes(grn$genes)
write.table(grn_genes, paste0("grn_results/grn_ds_genes_Q", q_cutoff, "_dPSI", dpsi_cutoff, ".txt"), sep="\t", quote=F, row.names = F)

print(paste("GRN:", dim(grn)[1]))

# C9orf72
c9orf <- read.table("c9or7f2_results/leafcutter_ds_cluster_significance.txt", sep="\t", header=T, stringsAsFactors = F)
c9orf_es <- read.table("c9or7f2_results/leafcutter_ds_effect_sizes.txt", sep="\t", header=T, stringsAsFactors = F)
c9orf <- c9orf[!is.na(c9orf$p),]
c9orf <- c9orf[c9orf$p.adjust <= q_cutoff,]
c9orf$clusterID <- getClusters(c9orf$cluster)
c9orf_es$clusterID <- getESClusters(c9orf_es$intron)
c9orf_es <- c9orf_es[c9orf_es$clusterID %in% c9orf$clusterID,]
# filter for dPSI
c9orf_es <- c9orf_es[abs(c9orf_es$deltapsi) >= dpsi_cutoff,]
c9orf <- c9orf[c9orf$clusterID %in% c9orf_es$clusterID,]
c9orf_genes <- unlistGenes(c9orf$genes)
write.table(c9orf_genes, paste0("c9or7f2_results/c9orf72_ds_genes_Q", q_cutoff, "_dPSI", dpsi_cutoff, ".txt"), sep="\t", quote=F, row.names = F)

print(paste("C9ORF72:", dim(c9orf)[1]))

