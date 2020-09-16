############################################################################
#
# Aggregation of cluster-based CAGE counts to gene-level counts
#
# 
############################################################################
library(stringr)
# load data
setwd("~/rimod/CAGE/results_annotation/")
cage <- read.csv("RiMod_genewise_annotated_CAGE_clusters.txt", sep="\t")

# Divide table in annotation and counts
cage <- cage[, c(-1, -2, -3, -4, -5)] # remove locations and stuff
annot <- cage[,249:ncol(cage)]
counts <- cage[,1:248]
#counts$geneId <- cage$geneId
# Get all genes
#genes <- as.character(levels(factor(cage$geneId)))
genes <- as.character(cage$geneId)

new_counts <- aggregate(counts, by=list(genes), FUN=sum)
genes <- str_split(new_counts$Group.1, pattern="[.]", simplify = T)[,1]
rownames(new_counts) <- genes
new_counts <- new_counts[, -1]


write.table(new_counts, "RiMod_aggrGeneCounts_CAGEseq_all.txt", quote=F, sep="\t")

# Normalize counts as CPM
library(edgeR)
y <- DGEList(new_counts)
cpms <- cpm(y)
write.table(cpms, "RiMod_aggrGeneCPM_CAGEseq_all.txt", sep="\t", quote=F)

# split counts in frontal and temporal
fro <- new_counts[, grepl("fro", colnames(new_counts))]
tem <- new_counts[, grepl("tem", colnames(new_counts))]
fro_cpm <- cpms[, grepl("fro", colnames(cpms))]
tem_cpm <- cpms[, grepl("tem", colnames(cpms))]

write.table(fro, "RiMod_aggrGeneCounts_CAGEseq_fro.txt", quote=F, sep="\t")
write.table(tem, "RiMod_aggrGeneCounts_CAGEseq_tem.txt", quote=F, sep="\t")

write.table(fro_cpm, "RiMod_aggrGeneCPM_CAGEseq_fro.txt", quote=F, sep="\t")
write.table(tem_cpm, "RiMod_aggrGeneCPM_CAGEseq_tem.txt", quote=F, sep="\t")
