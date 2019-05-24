##################################
# Analysis of LeafCutter results
##################################
setwd("~/rimod/RNAseq/as_analysis/leafcutter/")

# MAPT
mapt <- read.table("mapt_results/leafcutter_ds_cluster_significance.txt", sep="\t", header=T, stringsAsFactors = F)
mapt_es <- read.table("mapt_results/leafcutter_ds_effect_sizes.txt", sep="\t", header=T, stringsAsFactors = F)

print(dim(mapt))
mapt <- mapt[!is.na(mapt$p),]
print(dim(mapt))

test <- mapt[mapt$genes == 'MAPT',]
test <- test[!is.na(test$p),]

# Get MAPT cluster ids
clusters <- test$cluster
clusters <- as.character(sapply(clusters, function(x){strsplit(x, split=":")[[1]][[2]]}))

mapt_es$clusters <- as.character(sapply(mapt_es$intron, function(x){strsplit(x, split=":")[[1]][[4]]}))

es <- mapt_es[mapt_es$clusters %in% clusters,]

majiq <- read.table("../majiq/majiq_psi/mapt_mapt_psi.tsv", sep="\t", header=T, check.names = F)
majiq <- majiq[, c(-2, -3, -6, -9, -10, -11, -12, -13)]
dpsi <- as.character(majiq$`E(dPSI) per LSV junction`)
dpsi <- sapply(dpsi, function(x){strsplit(x, split=";")})

keep_vec = rep(FALSE, length(dpsi))
count = 0
for (elem in dpsi) {
  count = count + 1
  tmp = c()
  for (e in elem){
    tmp = c(tmp, as.numeric(e))
  }
  if (any(abs(tmp) > 0.05)){
    keep_vec[count] <- TRUE
  }
}

majiq <- majiq[keep_vec,]

