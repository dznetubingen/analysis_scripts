####################
# Format CAGE clusters for frontal and temporal samples
####################
setwd("~/rimod/CAGE/cage_analysis/tf_enrichment_analysis/")

# Load cluster matrix
mat <- read.table("Raw_All7RegionSamps_3kbGR_DF.txt", sep="\t", header=T)
rownames(mat) <- mat$Clus_ID

# Split in frontal and temporal
fro <- mat[, grepl("_fro_", colnames(mat))]
tem <- mat[, grepl("_tem_", colnames(mat))]

write.table(fro, "frontal_clusters_CAGE.txt", sep="\t", quote=F, col.names = NA)
write.table(tem, "temporal_clusters_CAGE.txt", sep="\t", quote=F, col.names = NA)