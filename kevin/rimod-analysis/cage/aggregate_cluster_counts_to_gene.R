############################################################################
#
# Aggregation of cluster-based CAGE counts to gene-level counts
#
# 
############################################################################

# load data
setwd("~/rimod/CAGE/cage_data//")
cage <- read.table("~/rimod/CAGE/cage_data/Raw_All7RegionSamps_3kbGR_DF.txt", sep="\t", header =T)

# Divide table in annotation and counts
annot <- cage[,1:14]
counts <- cage[,15:ncol(cage)]
counts$geneId <- cage$geneId
# Get all genes
genes <- as.character(levels(factor(cage$geneId)))

# Aggregate genes
rnames <- c()
df <- counts[1,]
df <- df[,-(ncol(df))]
for (g in genes) {
  sub <- counts[counts$geneId == g,]
  sub <- sub[,-(ncol(sub))]
  sub <- apply(sub, 2, sum)
  df <- rbind(df, sub)
  rnames <- c(rnames, g)
}
# Remove decoy column and write table
df <- df[-1,]
rownames(df) <- rnames
write.table(df, "cage_all7regions_3kbgr_aggr.txt", sep="\t", quote=F, col.names = NA)

# Normalize counts as CPM
library(edgeR)
y <- DGEList(df)
cpms <- cpm(y)
write.table(cpms, "cage_all7regions_3kbr_aggr_CPM.txt", sep="\t", quote=F, col.names=NA)
