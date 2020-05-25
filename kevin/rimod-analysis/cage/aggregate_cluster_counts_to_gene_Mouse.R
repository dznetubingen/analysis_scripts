############################################################################
#
# Aggregation of cluster-based CAGE counts to gene-level counts
#
# 
############################################################################
library(stringr)
# load data
setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/mouse/mapt_mice//")
cage <- read.csv("RiMod_MouseMAPT_genewise_annotated_CAGE_clusters.txt", sep="\t")

# Divide table in annotation and counts
cage <- cage[, c(-1, -2, -3, -4, -5)] # remove locations and stuff

# Keep only things close to promoter and 5'UTR
keep_regions = c("Promoter (1-2kb)", "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")
cage <- cage[cage$annotation %in% keep_regions,]

annot <- cage[,(ncol(cage)-7):ncol(cage)]
counts <- cage[,1:(ncol(cage)-8)]



#counts$geneId <- cage$geneId
# Get all genes
#genes <- as.character(levels(factor(cage$geneId)))
genes <- as.character(cage$geneId)

new_counts <- aggregate(counts, by=list(genes), FUN=sum)
genes <- str_split(new_counts$Group.1, pattern="[.]", simplify = T)[,1]
rownames(new_counts) <- genes
new_counts <- new_counts[, -1]


write.table(new_counts, "RiMod_MouseMAPT_aggrGeneCounts_CAGEseq_all.txt", quote=F, sep="\t")

# Normalize counts as CPM
library(edgeR)
y <- DGEList(new_counts)
cpms <- cpm(y)
write.table(cpms, "RiMod_MouseMAPT_aggrGeneCPM_CAGEseq_all.txt", sep="\t", quote=F)
