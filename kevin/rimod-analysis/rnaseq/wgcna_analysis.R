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
library(pheatmap)
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

setwd("~/rimod/RNAseq/analysis/wgcna_modules/")
###
# own testing
cor.df <- moduleTraitCor[, c("group", "mapt", "grn", "c9")]
pheatmap(cor.df ,color = viridis(200), filename = "WGCNA_correlation_heatmap.png")

# Make Heatmap for Figure 3
hm.df <- moduleTraitCor
colnames(hm.df) <- c("Age", "InNeurons", "Oligodendrocytes", "Endothelial", "Microglia", "Asctrocytes", "OPC", "ExNeurons", "FTD",
                     "FTD-MAPT", "FTD-GRN", "FTD-C9orf72", "Sex")
rownames(hm.df) <- gsub("ME", "", rownames(hm.df))
pheatmap(t(hm.df), color = viridis(200), filename = "WGCNA_correlation_heatmap_all.png", angle_col = "315",
         treeheight_col = 0, treeheight_row = 0, width = 6.3, height=4, cluster_rows =F)

# pvalue df
pval.df <- moduleTraitPvalue[, c("group", "mapt", "grn", "c9")]
pval.df[pval.df > 0.05] <- 1

pheatmap(t(pval.df), color = viridis(200, option = "A"), filename = "WGCNA_pvalue_heatmap.png", cluster_rows = F, cluster_cols = F,
         width=5, height=5, legend = F)


# others
# save all modules
mods = levels(df$module)
for (m in mods) {
  tmp = df[df$module == m,]
  write.table(tmp$genes, paste(m, "module.txt", sep="_"), row.names=F, col.names=F, quote=F)
}


# Which modules are MAPT, GRN and C9orf72 in
mapt.gene = "ENSG00000186868"
grn.gene = "ENSG00000030582"
c9orf72.gene = "ENSG00000147894"

mapt.mod = df[df$genes == mapt.gene,]
grn.mod = df[df$genes == grn.gene,]
c9.mod = df[df$genes == c9orf72.gene,]
###


pval.df <- moduleTraitPvalue
pval.df[pval.df > 0.05] <- 1
pheatmap(pval.df, color = viridis(200, option = "A"), filename = "WGCNA_pvalue_heatmap_all.png", cluster_rows = F, cluster_cols = F,
         width = 5, height=5)
pheatmap(t(moduleTraitCor), viridis(200, option = "A"), filename = "WGCNA_cor_heatmap_all.png", cluster_rows = F, cluster_cols = F,
         width=5, height=5)



# save the data frame
write.table(df, "WGCNA_modules.txt", sep="\t", quote=F)
write.table(moduleTraitCor, "WGCNA_moduleTraitCor.txt", sep="\t", quote=F)
write.table(moduleTraitPvalue, "WGCNA_moduleTraitPvalule.txt", sep="\t", quote=F)

#================================================#

#####
# Closer look into specific modules
#####
library(ggplot2)
library(biomaRt)
library(ggrepel)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
geneModuleMembership <- as.data.frame(cor(mat, MEs, use = "p"))

# load DE results
grn.deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_result_grn.ndc_fro_2019-10-23_13.33.11.txt", sep="\t", header=T, row.names=1)
c9.deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_result_c9.ndc_fro_2019-10-23_13.33.11.txt", sep="\t", header=T, row.names=1)
mapt.deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_result_mapt.ndc_fro_2019-10-23_13.33.11.txt", sep="\t", header=T, row.names=1)

### ORANGE MODULE and C9orf72
module = "orange"
trait = md$c9

geneTraitSignificance <- as.data.frame(cor(mat, trait, use="p"))
moduleGenes = moduleColors==module
moduleGenesMM = abs(geneModuleMembership[moduleGenes, match(module, modNames)])
moduleGenesSig = abs(geneTraitSignificance[moduleGenes, 1])
moduleGenesEnsemble = colnames(mat)[moduleGenes]

# get bm
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=moduleGenesEnsemble, mart=ensembl)
tmp = data.frame(Membership = moduleGenesMM, TraitCor = moduleGenesSig, Gene = moduleGenesEnsemble)
tmp <- merge(tmp, bm, by.x="Gene", by.y="ensembl_gene_id")
tmp <- merge(tmp, c9.deg ,by.x="Gene", by.y="row.names")

# Make labels
labels = tmp$hgnc_symbol
labels[tmp$Membership < 0.4] <- ""
labels[tmp$TraitCor < 0.2] <- ""
labels[tmp$padj > 0.1] <- ""
table(labels)

p = ggplot(tmp, aes(x=Membership, y=TraitCor, color=padj)) +
  geom_point(size = 5, alpha=0.4) + 
  theme_minimal() +
  scale_color_gradient(low = "#f0650e", high = "#0091ff") +
  geom_label_repel(aes(label=labels), box.padding = 0.35, point.padding = 0.5, segment.colour = 'grey50')
p

ggsave(filename = "c9orf72_orange.png", width=6, height=3.5)

### DARKOLIVEGREEN MODULE and C9orf72
module = "darkolivegreen"
trait = md$c9

geneTraitSignificance <- as.data.frame(cor(mat, trait, use="p"))
moduleGenes = moduleColors==module
moduleGenesMM = abs(geneModuleMembership[moduleGenes, match(module, modNames)])
moduleGenesSig = abs(geneTraitSignificance[moduleGenes, 1])
moduleGenesEnsemble = colnames(mat)[moduleGenes]

# get bm
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=moduleGenesEnsemble, mart=ensembl)
tmp = data.frame(Membership = moduleGenesMM, TraitCor = moduleGenesSig, Gene = moduleGenesEnsemble)
tmp <- merge(tmp, bm, by.x="Gene", by.y="ensembl_gene_id")
tmp <- merge(tmp, c9.deg ,by.x="Gene", by.y="row.names")

# Make labels
labels = tmp$hgnc_symbol
labels[tmp$Membership < 0.8] <- ""
labels[tmp$TraitCor < 0.3] <- ""
#labels[tmp$padj > 0.05] <- ""

p = ggplot(tmp, aes(x=Membership, y=TraitCor, color=padj)) +
  geom_point(size = 5, alpha=0.4) + 
  theme_minimal() +
  scale_color_gradient(low = "#f0650e", high = "#0091ff") +
  geom_label_repel(aes(label=labels), box.padding = 0.35, point.padding = 0.5, segment.colour = 'grey50')
p

ggsave(filename = "c9orf72_darkolivegreen.png", width=6, height=3.5)


### LIGHTCYAN MODULE and GRN
module = "lightcyan"
trait = md$grn

geneTraitSignificance <- as.data.frame(cor(mat, trait, use="p"))
moduleGenes = moduleColors==module
moduleGenesMM = abs(geneModuleMembership[moduleGenes, match(module, modNames)])
moduleGenesSig = abs(geneTraitSignificance[moduleGenes, 1])
moduleGenesEnsemble = colnames(mat)[moduleGenes]

# get bm
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=moduleGenesEnsemble, mart=ensembl)
tmp = data.frame(Membership = moduleGenesMM, TraitCor = moduleGenesSig, Gene = moduleGenesEnsemble)
tmp <- merge(tmp, bm, by.x="Gene", by.y="ensembl_gene_id")
tmp <- merge(tmp, grn.deg ,by.x="Gene", by.y="row.names")

# make labels
labels = tmp$hgnc_symbol
labels[tmp$Membership < 0.8] <- ""
labels[tmp$TraitCor < 0.3] <- ""
labels[tmp$padj >= 0.01] <- ""

p = ggplot(tmp, aes(x=Membership, y=TraitCor, color=padj)) +
  geom_point(size = 5, alpha=0.4) + 
  theme_minimal() +
  scale_color_gradient(low = "#f0650e", high = "#0091ff") +
  geom_label_repel(aes(label=labels), box.padding = 0.35, point.padding = 0.5, segment.colour = 'grey50')
p
ggsave(filename="grn_lightcyan_module.png", width=6, height=3.5)


### PALETURQUOISE MODULE and GRN
module = "paleturquoise"
trait = md$grn

geneTraitSignificance <- as.data.frame(cor(mat, trait, use="p"))
moduleGenes = moduleColors==module
moduleGenesMM = abs(geneModuleMembership[moduleGenes, match(module, modNames)])
moduleGenesSig = abs(geneTraitSignificance[moduleGenes, 1])
moduleGenesEnsemble = colnames(mat)[moduleGenes]

# get bm
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=moduleGenesEnsemble, mart=ensembl)
tmp = data.frame(Membership = moduleGenesMM, TraitCor = moduleGenesSig, Gene = moduleGenesEnsemble)
tmp <- merge(tmp, bm, by.x="Gene", by.y="ensembl_gene_id")
tmp <- merge(tmp, grn.deg ,by.x="Gene", by.y="row.names")

# make labels
labels = tmp$hgnc_symbol
labels[tmp$Membership < 0.8] <- ""
labels[tmp$TraitCor < 0.3] <- ""
labels[tmp$padj >= 0.01] <- ""

p = ggplot(tmp, aes(x=Membership, y=TraitCor, color=padj)) +
  geom_point(size = 5, alpha=0.4) + 
  theme_minimal() +
  scale_color_gradient(low = "#f0650e", high = "#0091ff") +
  geom_label_repel(aes(label=labels), box.padding = 0.35, point.padding = 0.5, segment.colour = 'grey50')
p
ggsave(filename="grn_lpaleturquoise_module.png", width=6, height=3.5)

### PINK MODULE and MAPT
module = "pink"
trait = md$mapt

geneTraitSignificance <- as.data.frame(cor(mat, trait, use="p"))
moduleGenes = moduleColors==module
moduleGenesMM = abs(geneModuleMembership[moduleGenes, match(module, modNames)])
moduleGenesSig = abs(geneTraitSignificance[moduleGenes, 1])
moduleGenesEnsemble = colnames(mat)[moduleGenes]

# get bm
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=moduleGenesEnsemble, mart=ensembl)
tmp = data.frame(Membership = moduleGenesMM, TraitCor = moduleGenesSig, Gene = moduleGenesEnsemble)
tmp <- merge(tmp, bm, by.x="Gene", by.y="ensembl_gene_id")
tmp <- merge(tmp, mapt.deg ,by.x="Gene", by.y="row.names")

# make labels
labels = tmp$hgnc_symbol
labels[tmp$Membership < 0.8] <- ""
labels[tmp$TraitCor < 0.3] <- ""
labels[tmp$padj >= 0.01] <- ""

p = ggplot(tmp, aes(x=Membership, y=TraitCor, color=padj)) +
  geom_point(size = 5, alpha=0.4) + 
  theme_minimal() +
  scale_color_gradient(low = "#f0650e", high = "#0091ff") +
  geom_label_repel(aes(label=labels), box.padding = 0.35, point.padding = 0.5, segment.colour = 'grey50')
p

ggsave(filename="mapt_pink_module.png", width=6, height=3.5)

### BROWN MODULE and GRN
module = "brown"
trait = md$grn

geneTraitSignificance <- as.data.frame(cor(mat, trait, use="p"))
moduleGenes = moduleColors==module
moduleGenesMM = abs(geneModuleMembership[moduleGenes, match(module, modNames)])
moduleGenesSig = abs(geneTraitSignificance[moduleGenes, 1])
moduleGenesEnsemble = colnames(mat)[moduleGenes]

# get bm
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=moduleGenesEnsemble, mart=ensembl)
tmp = data.frame(Membership = moduleGenesMM, TraitCor = moduleGenesSig, Gene = moduleGenesEnsemble)
tmp <- merge(tmp, bm, by.x="Gene", by.y="ensembl_gene_id")
tmp <- merge(tmp, grn.deg ,by.x="Gene", by.y="row.names")

# make labels
labels = tmp$hgnc_symbol
labels[tmp$Membership < 0.8] <- ""
labels[tmp$TraitCor < 0.3] <- ""
labels[tmp$padj >= 0.001] <- ""

p = ggplot(tmp, aes(x=Membership, y=TraitCor, color=padj)) +
  geom_point(size = 5, alpha=0.4) + 
  theme_minimal() +
  scale_color_gradient(low = "#f0650e", high = "#0091ff") +
  geom_label_repel(aes(label=labels), box.padding = 0.35, point.padding = 0.5, segment.colour = 'grey50')
p

ggsave(filename="grn_brown_module.png", width=6, height=3.5)


### lightcyan MODULE and GRN
module = "lightcyan"
trait = md$grn

geneTraitSignificance <- as.data.frame(cor(mat, trait, use="p"))
moduleGenes = moduleColors==module
moduleGenesMM = abs(geneModuleMembership[moduleGenes, match(module, modNames)])
moduleGenesSig = abs(geneTraitSignificance[moduleGenes, 1])
moduleGenesEnsemble = colnames(mat)[moduleGenes]

# get bm
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=moduleGenesEnsemble, mart=ensembl)
tmp = data.frame(Membership = moduleGenesMM, TraitCor = moduleGenesSig, Gene = moduleGenesEnsemble)
tmp <- merge(tmp, bm, by.x="Gene", by.y="ensembl_gene_id")
tmp <- merge(tmp, grn.deg ,by.x="Gene", by.y="row.names")

# make labels
labels = tmp$hgnc_symbol
labels[tmp$Membership < 0.8] <- ""
labels[tmp$TraitCor < 0.3] <- ""
labels[tmp$padj >= 0.001] <- ""

p = ggplot(tmp, aes(x=Membership, y=TraitCor, color=padj)) +
  geom_point(size = 5, alpha=0.4) + 
  theme_minimal() +
  scale_color_gradient(low = "#f0650e", high = "#0091ff") +
  geom_label_repel(aes(label=labels), box.padding = 0.35, point.padding = 0.5, segment.colour = 'grey50')
p

ggsave(filename="grn_lightcyan_module.png", width=6, height=3.5)



