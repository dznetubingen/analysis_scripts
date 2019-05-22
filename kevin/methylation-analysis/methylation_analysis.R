# Methylation Data Testing
# Pipeline generation according to Maksimovic et al., 2017
# Note: This script is specific to the RIMOD methylation data set
# It is not a general purpose script and should not be used for ohter data as is but can be used as guidance
# Ask Kevin if you have question :)

# load packages required for analysis 
# TODO: ADJUST THIS 
library(limma)
library(RColorBrewer)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(missMethyl) 
library(matrixStats) 
library(minfiData) 
library(Gviz) 
library(DMRcate) 
library(stringr)
library(viridis)
library(plyr)
library(ggplot2)
library(quantro)


######
# Parsing and QC
#######

data.dir <- "/media/phdstudent/Data/rimod/Methylation/ON-2017_8527_ILLUMINA/"
setwd(data.dir)
list.files(data.dir)
design <- read.csv(paste0(data.dir,"rimod_design.csv"))
# Only include temporal samples
design <- design[design$REGION == "temporal",]

# Get annotation
annEpicObj <- getAnnotationObject(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
annEpic <- getAnnotation(annEpicObj)

# read sample sheet and manually adapt it
targets <- read.metharray.sheet(data.dir, pattern="samplesheet.csv$")
colnames(targets)[7] <- "Slide"
colnames(targets)[8] <- "Array"
# extend targets
colnames(design)[4] <- "Sample_Name"
mdata <- join(targets, design, by="Sample_Name", type="left")
targets <- mdata


# Construct basenames
targets$Basename <- file.path(data.dir, targets$Slide, paste0(targets$Slide, "_",targets$Array))
# read data
RGset <- read.metharray.exp(targets = targets)

# Calculate detection p-values
## Note: detection p-values are calculated by comparing the total intensity to the background (negative control probes)
## High p-values indicate poor quality
detP <- detectionP(RGset)
# Examine mean detection p-values across all samples
pal <- viridis(48)
png("Detection_Pvalues.png")
barplot(colMeans(detP),  col=pal[factor(targets$Sample_Name)], las=2, cex.names = 0.8, ylab = "Mean Detection p-values", main="Detection p-values")
abline(h=0.01, col="red")
dev.off()

# Generate minfi QC report
group <- as.character(targets$DISEASE.CODE)
qcReport(RGset, sampNames = targets$Sample_Name, sampGroups = group, pdf="QC_report.pdf")

# remove samples with bad detection p-values
keep <- colMeans(detP) < 0.05
RGset <- RGset[,keep]
#####
# quantro analysis
#####
# This analysis is done to guide in choice of the normalization
# mSet <- preprocessRaw(RGset)
# beta.vals <- getBeta(mSet)
# # Generate group factor
# group <- as.character(targets$DISEASE.CODE)
# for(i in 1:length(group)){
#   new_group <- group
#   if (is.na(group[i])){
#     new_group[i] <- "unknown"
#   }
#   group <- as.factor(new_group)
# }
# 
# 
# 
# matdensity(beta.vals, groupFactor = group, main="Beta Values")
# legend('top', levels(group))
# matboxplot(beta.vals, groupFactor = group, main="Beta Values")
# 
# # quantro function
# qtest <- quantro(object=mSet, groupFactor = group, B=100)


#####
# Normalization
#####
# the result from the quantro function can be used to choose the normalization procedure

# Apply preprocessQuantile function
mSetSq <- preprocessQuantile(RGset)
mSetRaw <- preprocessRaw(RGset)
mSetFn <- preprocessFunnorm(RGset)

# Create density plots of the Beta values to compare raw with normalized values
png("raw_density_plot.png")
densityPlot(getBeta(mSetRaw), main ="Raw")
dev.off()
png("normalized_density_plot.png")
densityPlot(getBeta(mSetSq), main="Normalized")
dev.off()
png("normalized_functional_density_plot", width=800, height=500)
densityPlot(getBeta(mSetFn), main="Funnorm normalized")
dev.off()


####
# MDS plot
####
# Find confounding factors using MDS analysis (basically PCA)
# Color by Disease Code
pal <- brewer.pal(4, "Dark2")
png("mds_disease_code12.png")
plotMDS(getM(mSetSq), top=1000, gene.selection="common" , labels= targets$DISEASE.CODE, col=pal[factor(targets$DISEASE.CODE)])
dev.off()
png("mds_disease_code13.png")
plotMDS(getM(mSetSq), dim=c(1,3), top=1000, gene.selection="common" , labels= targets$DISEASE.CODE, col=pal[factor(targets$DISEASE.CODE)])
dev.off()
png("mds_disease_code23.png")
plotMDS(getM(mSetSq), dim=c(2,3), top=1000, gene.selection="common" , labels= targets$DISEASE.CODE, col=pal[factor(targets$DISEASE.CODE)])
dev.off()

# Color by gender
pal <- brewer.pal(2, "Dark2")
png("mds_gender12.png")
plotMDS(getM(mSetSq), top=1000, gene.selection="common" , labels= targets$GENDER, col=pal[factor(targets$GENDER)])
dev.off()

##############
### FILTERING
##############
# Filtering can be applied using various metrics. Here SNPs, sex chromosomes and detP values are used
# Filter probes based on detection p-values
detP <- detP[match(featureNames(mSetSq), rownames(detP)),]
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
mSetSqFlt <- mSetSq[keep,]

## Ignore methylation on sex-chromosomes
keep <- !(featureNames(mSetSqFlt) %in% annEpic$Name[annEpic$chr %in% c("chrX", "chrY")])
mSetSqFlt <- mSetSqFlt[keep,]

# Remove probes with SNPs at CpG Site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)

# Potentially remove cross-reactive probes
# TODO download correct file for this
# edit: find file for EPIC array

# Redo the MDS plots to see what has changed
# Color by Disease Code
pal <- brewer.pal(4, "Dark2")
png("mds_disease_code12_filtered.png")
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common" , labels= targets$DISEASE.CODE, col=pal[factor(targets$DISEASE.CODE)])
dev.off()
png("mds_disease_code13_filtered.png")
plotMDS(getM(mSetSqFlt), dim=c(1,3), top=1000, gene.selection="common" , labels= targets$DISEASE.CODE, col=pal[factor(targets$DISEASE.CODE)])
dev.off()
png("mds_disease_code23_filtered.png")
plotMDS(getM(mSetSqFlt), dim=c(2,3), top=1000, gene.selection="common" , labels= targets$DISEASE.CODE, col=pal[factor(targets$DISEASE.CODE)])
dev.off()
# Color by gender
pal <- brewer.pal(2, "Dark2")
png("mds_gender12_filtered.png")
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common" , labels= targets$GENDER, col=pal[factor(targets$GENDER)])
dev.off()


##############
### Differential Methylation Analysis
##############
# Get M-values for analysis and Beta-values for visualization
mVals <- getM(mSetSqFlt)
betaVals <- getBeta(mSetSqFlt)

# Remove samples with no complete phenotype information
group <- as.character(targets$DISEASE.CODE)
keep <- !is.na(group)
group <- group[keep]
mVals <- mVals[,keep]
betaVals <- betaVals[,keep]
targets <- targets[keep,]

G <- factor(group)
design.matrix <- model.matrix(~ 0 + G)
colnames(design.matrix) <- make.names(levels(G))

# Create contrasts matrix
conts <- c("FTD.C9-control", "FTD.GRN-control", "FTD.MAPT-control")
contrast.matrix <- makeContrasts(contrasts = conts, levels = design.matrix)

fit <- lmFit(mVals, design = design.matrix)
cont.fit <- contrasts.fit(fit = fit, contrasts = contrast.matrix)
fit2 <- eBayes(cont.fit)
res <- decideTests(fit2)
summary(res)

######
# Explore methylation of contrasts
# Here: FTD-MAPT
# Get sub-annotation 
pvalCutoff <- 0.05
annEpic <- annEpic[match(rownames(mVals), annEpic$Name), ]
dmp.mapt <- topTable(fit2, num=Inf, coef="FTD.MAPT-control", genelist = annEpic)
dmp.mapt <- dmp.mapt[dmp.mapt$adj.P.Val <= pvalCutoff,] # only consider significant genes
# Genomic region of the CpG
mapt.region <- as.character(dmp.mapt$UCSC_RefGene_Group)
# Region is annotated for every transcript. To reduce complexity only consider first listed region annotation
new.region <- c()
for (i in 1:length(mapt.region)){
  tmp <- strsplit(mapt.region[i], split=";")[[1]]
  if (length(tmp) > 1){
    tmp <- tmp[1]
  }
  new.region <- c(new.region, tmp)
}
mapt.region <- new.region
barplot(table(mapt.region))

# Island region testing
mapt.island <- dmp.mapt$Relation_to_Island
barplot(table(mapt.island))

# HMM predicted island testing
mapt.predisl <- dmp.mapt$HMM_Island
tmp.isl <- c()
for (i in 1:length(mapt.predisl)){
  if (mapt.predisl[i] == ""){
    tmp <- "OpeanSea"
  } else {
    tmp <- "HMM_Island"
  }
  tmp.isl <- c(tmp.isl, tmp)
}
mapt.predisl <- tmp.isl
barplot(table(mapt.predisl))

###########
### Look at distribution of regions for whole EPIC array
##########
# Genomic region of the CpG
epic.region  <- as.character(annEpic$UCSC_RefGene_Group)
# Region is annotated for every transcript. To reduce complexity only consider first listed region annotation
new.region <- c()
for (i in 1:length(epic.region)){
  tmp <- strsplit(epic.region[i], split=";")[[1]]
  if (length(tmp) > 1){
    tmp <- tmp[1]
  }
  new.region <- c(new.region, tmp)
}
epic.region <- new.region
barplot(table(epic.region))

########
## Differentially Methylated Regions Analysis
########
# Create annotation objects
myAnnotC9 <- cpg.annotate(object = mVals, datatype = "array", what = "M", analysis.type = "differential", design = design.matrix,
                            contrasts = T, cont.matrix = contrast.matrix, coef = "FTD.C9-control", arraytype = "EPIC")
myAnnotGrn <- cpg.annotate(object = mVals, datatype = "array", what = "M", analysis.type = "differential", design = design.matrix,
                            contrasts = T, cont.matrix = contrast.matrix, coef = "FTD.GRN-control", arraytype = "EPIC")
myAnnotMapt <- cpg.annotate(object = mVals, datatype = "array", what = "M", analysis.type = "differential", design = design.matrix,
                        contrasts = T, cont.matrix = contrast.matrix, coef = "FTD.MAPT-control", arraytype = "EPIC")

# Calculate differentially methylated regions
dmr.c9 <- dmrcate(myAnnotC9, lambda = 1000, C=2)
dmr.grn <- dmrcate(myAnnotGrn, lambda = 1000, C=2)
dmr.mapt <- dmrcate(myAnnotMapt, lambda=1000, C=2)

# Further processing
results.ranges <- extractRanges(dmr.mapt, genome="hg19")
pal <- brewer.pal(4, "Dark2")
groups <- pal
names(groups) <- colnames(design.matrix)
cols <- groups[as.character(factor(colnames(design.matrix)))]
# DMR.plot is currently not working
# png("testplot.png")
# DMR.plot(ranges = results.ranges, dmr = 1, CpGs = betaVals, what="Beta", phen.col = cols,
#          arraytype = "EPIC", pch=16, toscale=TRUE, plotmedians=TRUE, genome="hg19", samps=c(1:43))
# dev.off()

# Try plotting using the Gviz package

genome <- "hg19"
dmrIdx = 1
# extract coordinates
coords <- strsplit2(dmr.mapt$results$coord[dmrIdx],":")
chrom <- coords[1]
start <- as.numeric(strsplit2(coords[2], "-")[1])
end <- as.numeric(strsplit2(coords[2], "-")[2])
# add extra space to plot
minbase = start - (0.25*(end-start))
maxbase = end + (0.25*(end-start))

# Create the tracks
iTrack <- IdeogramTrack(genome = genome, chromosome = chrom, name="")
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
rTrack <- UcscTrack(genome = genome, chromosome = chrom, track ="refGene", 
                    from = minbase, to = maxbase, trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
                    symbol="name2", transcript="name", strand="strand", fill="darkblue", stacking="squish", name="RefSeq",
                    showId = T, geneSymbol=T)

# Check ordering
annEpicSubOrd <- annEpicSub[order(annEpicSub$chr, annEpicSub$pos),]
bValsOrd <- betaVals[match(annEpicSubOrd$Name,rownames (betaVals)),]

# Create data tracks
cpgData <- GRanges(seqnames = Rle(annEpicSubOrd$chr),
                   ranges = IRanges(start=annEpicSubOrd$pos, end = annEpicSubOrd$pos),
                   strand = Rle(rep("*", nrow(annEpicSubOrd))),
                   betas <- bValsOrd)
# extract data on CpGs in DMR
cpgData <- subsetByOverlaps(cpgData, results.ranges[dmrIdx])
# methylation data track
methTrack <- DataTrack(range=cpgData, groups = group, genome = genome,
                       chromosome = chrom, ylim = c(-0.05, 1.05, col = pal),
                       type=c("a", "p"), name="DNA methy \n beta value", background.panel="white",
                       legend=T, cex.title = 0.8, cex.axis = 0.8, cex.legend = 0.8)
# DMR position track
dmrTrack <- AnnotationTrack(start = start, end = end, genome = genome, name = "DMR", chromosome = chrom, fill = "darkred")

tracks <- list(iTrack, gTrack, rTrack)
sizes <- c(2,2,5)
plotTracks(tracks, from = minbase, to = maxbase, showTitle=TRUE, add53=TRUE, add35=TRUE, grid=T, lty.grid=3, length=tracks)


