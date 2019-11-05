######
# RNA-seq WGCNA analysis
######
library(WGCNA)
library(stringr)
library(dynamicTreeCut)
library(GO.db)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(viridis)
setwd("~/rimod/RNAseq/analysis/")

# multithreading
enableWGCNAThreads()

# Load exression data
mat <- read.table("RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_vst_values_2019-10-23_13.33.11.txt", sep="\t", row.names=1, header=T)
rownames(mat) <- str_split(rownames(mat), patter="[.]", simplify = T)[,1]
datExpr <- mat

# Load metadata
md <- read.csv("~/rimod/files/FTD_Brain.csv", stringsAsFactors = FALSE)
md$SAMPLEID <- as.character(sapply(md$SAMPLEID, function(x){strsplit(x, split="_")[[1]][[1]]}))
md$SAMPLEID <- str_pad(md$SAMPLEID, width = 5, side = "left", pad = "0") # fill sample ids to 5 digits

# bring counts and md in similar format
rna.samples <- as.character(sapply(colnames(mat), function(x){strsplit(x, split="_")[[1]][[1]]}))
rna.samples <- str_pad(gsub("X", "", rna.samples), width=5, side='left', pad='0')
md <- md[md$SAMPLEID %in% rna.samples,]
md <- md[match(rna.samples, md$SAMPLEID),]

# Load deconvolution results
fracs <- read.table("deconvolution/cdn_predictions.txt", sep="\t", header = T, row.names = 1)
rownames(fracs) <- str_split(rownames(fracs), pattern="_", simplify = T)[,1]
rownames(fracs) <- str_pad(gsub("X", "", rownames(fracs)), width=5, side="left", pad='0')
fracs <- fracs[match(rna.samples, rownames(fracs)),]
md <- cbind(md, fracs)

# keep only interesteing metadata
md <- md[, c(-1, -2, -3, -4, -5, -6, -9, -11, -12, -13, -14, -15, -16, -17, -20 ,-21, -23, -24, -25, -26, -27)]
md <- md[, c(-7, -8, -9, -10)]

# encode some values as binary
md$group = rep(0, nrow(md))
md$group[!md$CASE.CONTROL == "control"] <- 1
md$mapt <- rep(0, nrow(md))
md$mapt[md$DISEASE.CODE == "FTD-MAPT"] <- 1
md$grn <- rep(0, nrow(md))
md$grn[md$DISEASE.CODE == "FTD-GRN"] <- 1
md$c9 <- rep(0, nrow(md))
md$c9[md$DISEASE.CODE == "FTD-C9"] <- 1
md$gender <- rep(0, nrow(md))
md$gender[md$GENDER == "F"] <- 1

md <- md[, c(-1,-2,-3, -5)]
md <- md[, -2]

# Remove lowly-expressed genes
keep <- rowSums(as.matrix(mat) >= 5) >= 10
mat <- mat[keep,]
# transpose matrix
mat <- t(mat)


#####
# Choose soft-thresholding power
#####
sft <- pickSoftThreshold(mat)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab = "SoftThreshold (power)",ylab="Scale Free Topolog Model fit, signed Rr",
     main = "Scale independence")
abline(h=0.9, col='red')

plot(sft$fitIndices[,1], sft$fitIndices[,5], ylab="Mean connectivity", xlab="Soft threshold")

# Based on analysis, choose softPower = 5
softPower = 5
#========================================#


###
# Calculate adjacencies
###

adjacency = adjacency(mat, power = softPower)

#=====================#

##
# Topological Overlap Matrix (TOM) and clustering
##
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

# Clustering with TOM
geneTree = hclust(as.dist(dissTOM), method="average")
plot(geneTree, xlab="", ylab="", main="Gene clustering on TOM-based dissimilarity", labels=FALSE, hang=0.04)

#===============#

###
# Module seletion by dynamic tree cut
###
minModuleSize = 25
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
table(dynamicMods)

# minSize 25 yields 41 modules (smallest has 46 genes in it)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, dendroLabels = FALSE)
#=====================#

###
# Merge similar modules
###
# Caclulate eigengenes
MEList = moduleEigengenes(mat, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot result
plot(METree, main = "Clustering of module eigengenes", xlab="", sub="")

# Choose cutoff 
MEDissTres = 0.1 # 0.1 corresponds to 0.9 correlation
abline(h=MEDissTres, col='red')

# do the merging
merge = mergeCloseModules(mat, dynamicColors, cutHeight = MEDissTres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

# check result in plot
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = F)

# Save the results
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
save(MEs, moduleLabels, moduleColors, geneTree, file = "Merged_Modules_Tresh0.1_WGCNA.RData")
#===========================================================#
# Load again if necessary
load("Merged_Modules_Tresh0.1_WGCNA.RData")

####
# Quantification of module-trait associations
####
# number of genes and samples
nGenes <- ncol(mat)
nSamples <- nrow(mat)

# Recalculate MEs with color lables
MEsO = moduleEigengenes(mat, moduleColors)$eigengenes
MEs = orderMEs(MEsO)

# calculate correlation
moduleTraitCor = cor(MEs, md, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Plotting
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep="")
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(md),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = viridis(250),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5, zlim = c(-1, 1), main = "Module-trait relationships")

#=================================#

####
# Gene relationship to trait and important modules
####
group = as.data.frame(md$group)
names(group) <- "group"
modNames = substring(names(MEs), 3)
# gene Module membership
geneModuleMembership <- as.data.frame(cor(md, MEs, use="p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(mat, group, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(group), sep="")
names(GSPvalue) = paste("p.GS", names(group), sep="")

## Create Dataframe with information
geneInfo0 = data.frame(geneID = colnames(mat),
                       moduleColors = moduleColors,
                       geneTraitSignificance[,1],
                       GSPvalue[,1])
modOrder = order(-abs(cor(MEs, group, use="p")))
geneInfo = geneInfo0[modOrder,]

#===============================================#

###
# Get entrez IDs
library(biomaRt)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
bm <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), filters="ensembl_gene_id", values=colnames(mat), mart = ensembl)

# make lists of entrez gene IDs and modules
df <- data.frame(genes = colnames(mat), module = moduleColors)
df <- merge(df, bm, by.x="genes", by.y="ensembl_gene_id")
df <- df[!duplicated(df$genes),]

# Perform enrichment
go.enr <- GOenrichmentAnalysis(df$module, df$entrezgene_id, organism = "human", nBestP = 10)
