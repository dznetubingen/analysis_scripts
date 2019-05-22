#############################################
# Analysis of frontal FTD methylation data  #
#############################################

# Load packages
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
library(stringr)

# Get annotation
annEpicObj <- getAnnotationObject(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
annEpic <- getAnnotation(annEpicObj)

###############
# Read in Data
###############
# Data directory
data.dir <- "~/rimod/Methylation/frontal_methylation_0818"
setwd(data.dir)

#=== Read design file and format it ===#
design_file = paste0(data.dir, "/FTD_methylation_July2018.csv")
design <- read.csv(design_file)
colnames(design) <- c("SampleID", "Group")
# Change Sample A144/12
samples = as.character(design$SampleID)
sid = which(samples == "A144/12")
samples[sid] <- "14412"
design$SampleID = samples
# Merge with other design matrix
md <- read.csv("~/rimod/files/FTD_Brain.csv")
#md <- md[md$REGION == "frontal",]
sids <- as.character(md$SAMPLEID)
idx <- which(sids == "A144_12")
sids[idx] <- "14412"
sids <- str_pad(sids, 5, side="left", pad="0")
md$SampleID <- sids
md <- data.frame(age = md$AGE, gender = md$GENDER, sampleid = md$SampleID, gene=md$GENE)
mdata <- merge(design, md, by.x="SampleID", by.y="sampleid")
design <- mdata[!duplicated(mdata),]
# ========================================#

#== Read in the data ==#
# read sample sheet and manually adapt it
targets <- read.metharray.sheet(data.dir, pattern="2.csv$")
colnames(targets)[7] <- "Slide"
colnames(targets)[8] <- "Array"

# extend targets
sample_names = as.character(targets$Sample_Name)
samples <- substr(sample_names, 1, 5)
targets$Sample_Name2 <- samples
colnames(design)[1] <- "Sample_Name2"
mdata <- join(targets, design, by="Sample_Name2", type="left")
targets <- mdata

# Construct basenames
data_file_dir = paste0(data.dir, "/Data")
targets$Basename <- file.path(data_file_dir, targets$Slide, paste0(targets$Slide, "_",targets$Array))
# read data
RGset <- read.metharray.exp(targets = targets)
#===================================================#


#=== Quality Control ===#
detP <- detectionP(RGset)
# Examine mean detection p-values across all samples
pal <- viridis(48)
png("Detection_Pvalues.png")
barplot(colMeans(detP),  col=pal[factor(targets$Sample_Name)], las=2, cex.names = 0.8, ylab = "Mean Detection p-values", main="Detection p-values")
abline(h=0.01, col="red")
dev.off()

# Generate minfi QC report
group <- as.character(targets$Group)
qcReport(RGset, sampNames = targets$Sample_Name, sampGroups = group, pdf="QC_report.pdf")

# remove samples with bad detection p-values
keep <- colMeans(detP) < 0.01
RGset <- RGset[,keep]
targets <- targets[keep,]
#=====================================================#


#=== Normalization ===#
# Apply preprocessQuantile function as it is more suited for samples with largely different
# expression profiles
mSetFn <- preprocessFunnorm(RGset)
mSetRaw <- preprocessRaw(RGset) # for plotting

# Create density plots of the Beta values to compare raw with normalized values
png("raw_density_plot.png", width=800, height=500)
densityPlot(getBeta(mSetRaw), main ="Raw")
dev.off()
png("normalized_functional_density_plot", width=800, height=500)
densityPlot(getBeta(mSetFn), main="Funnorm normalized")
dev.off()
#======================================================#

#=== Filtering ===#
# Filtering can be applied using various metrics. Here SNPs, sex chromosomes and detP values are used

# Filter probes based on detection p-values
detP <- detP[match(featureNames(mSetFn), rownames(detP)),]
keep <- rowSums(detP < 0.01) == ncol(mSetFn)
mSetFnFlt <- mSetFn[keep,]

## Ignore methylation on sex-chromosomes
keep <- !(featureNames(mSetFnFlt) %in% annEpic$Name[annEpic$chr %in% c("chrX", "chrY")])
mSetFnFlt <- mSetFnFlt[keep,]

# Remove probes with SNPs at CpG Site
mSetFnFlt <- dropLociWithSnps(mSetFnFlt)


#=======================================================#


#=== MDS plotting ===#
# Find confounding factors using MDS analysis (basically PCA)
# Color by Disease Code
mds_dir <- "mds_plots/"
pal <- brewer.pal(4, "Dark2")
png("mds_plots/mds_disease_code12.png")
plotMDS(getM(mSetFnFlt), top=1000, gene.selection="common" , labels= targets$Group, col=pal[factor(targets$Group)])
dev.off()
png("mds_plots/mds_disease_code13.png")
plotMDS(getM(mSetFnFlt), dim=c(1,3), top=1000, gene.selection="common" , labels=targets$Group, col=pal[factor(targets$Group)])
dev.off()
png("mds_plots/mds_disease_code23.png")
plotMDS(getM(mSetFnFlt), dim=c(2,3), top=1000, gene.selection="common" , labels= targets$Group, col=pal[factor(targets$Group)])
dev.off()

# Color by Slide
pal <- brewer.pal(2, "Dark2")
png("mds_plots/mds_slide12.png")
plotMDS(getM(mSetFnFlt), top=1000, gene.selection="common" , labels= targets$Slide, col=pal[factor(targets$Slide)])
dev.off()

# Color by Array
pal <- brewer.pal(2, "Dark2")
png("mds_plots/mds_array12.png")
plotMDS(getM(mSetFnFlt), top=1000, gene.selection="common" , labels= targets$Array, col=pal[factor(targets$Array)])
dev.off()

# Color by Gender
pal <- brewer.pal(2, "Dark2")
png("mds_plots/mds_gender12.png")
plotMDS(getM(mSetFnFlt), top=1000, gene.selection="common" , labels= targets$gender, col=pal[factor(targets$gender)])
dev.off()

# Color by mutation
pal <- brewer.pal(10, "Dark2")
png("mds_plots/mds_disease_code12.png")
mutation <- as.character(targets$gene)
mutation[is.na(mutation)] <- "ndc"
plotMDS(getM(mSetFnFlt), top=1000, gene.selection="common" , labels= mutation, col=pal[factor(targets$Group)])
dev.off()
#======================================================#


#=== Differential Methylation Analysis ===#
mVals <- getM(mSetFnFlt)
betaVals <- getBeta(mSetFnFlt)

# Apply additional filtering using DMRcate
mVals <- rmSNPandCH(mVals, rmcrosshyb = TRUE, rmXY = TRUE)
betaVals <- rmSNPandCH(betaVals, rmcrosshyb = TRUE, rmXY = TRUE) 

# Create matrix for external use
cnames <- paste(targets$Sample_Name2, targets$Group, sep="_")
mvals_ext <- mVals
colnames(mvals_ext) <- cnames
write.table(mvals_ext, "mVals_matrix_frontal_methylation.txt", quote=F, sep="\t")

# Create Matrix
targets$Group[targets$Group == "FTD_MAPT"] <- "FTD-MAPT"
targets$Group <- make.names(targets$Group)
G <- factor(targets$Group)
A <- targets$age
design.matrix <- model.matrix(~ 0 + G + A)
colnames(design.matrix) <- c("FTD.C9", "FTD.GRN", "FTD.MAPT", "NDC", "AGE")

# Create contrasts matrix
conts <- c("FTD.C9-NDC", "FTD.GRN-NDC", "FTD.MAPT-NDC")
contrast.matrix <- makeContrasts(contrasts = conts, levels = design.matrix)
# Limma fitting
fit <- lmFit(mVals, design = design.matrix)
cont.fit <- contrasts.fit(fit = fit, contrasts = contrast.matrix)
fit2 <- eBayes(cont.fit)
res <- decideTests(fit2)
summary(res)

# Extract MAPT DMPs
annEpicSub <- annEpic[match(rownames(mVals),annEpic$Name), c(1:4, 12:19, 24:ncol(annEpic))]
dmps_mapt <- topTable(fit2, num=Inf, coef=3, genelist = annEpicSub)
write.table(dmps_mapt, "DMPs_mapt.ndc_funnNorm.txt")

# plot the top 4 most significantly differentially methylated CpGs 
par(mfrow=c(2,2))
sapply(rownames(dmps_mapt)[1:4], function(cpg){
  plotCpg(betaVals, cpg=cpg, pheno=targets$Group, ylab = "Beta values")
})
par(mfrow=c(1,1))

# Save mvals

write.table(mVals, "mVals_")
#============================================#

#=== Sub-classifying MAPT ===#
keep_groups <- c("NDC", "FTD.MAPT")
keep <- targets$Group %in% keep_groups
mapt.mvals <- mVals[,keep]
mapt.targets <- targets[keep,]
# Separate mutations
mapt.mut <- as.character(mapt.targets$gene)
mapt.mut[is.na(mapt.mut)] <- "ndc"
mapt.mut[!mapt.mut %in% c("P301L", "ndc")] <- "other"

# Create Matrix
G <- factor(mapt.mut)
A <- mapt.targets$age
design.matrix.mapt <- model.matrix(~ 0 + G + A)
colnames(design.matrix.mapt) <- c("ndc", "other", "P301L", "A")
# Create contrasts matrix
conts <- c("P301L-ndc", "other-ndc")
contrast.matrix.mapt <- makeContrasts(contrasts = conts, levels = design.matrix.mapt)
# Limma fitting
fit <- lmFit(mapt.mvals, design = design.matrix.mapt)
cont.fit <- contrasts.fit(fit = fit, contrasts = contrast.matrix.mapt)
fit2 <- eBayes(cont.fit)
res <- decideTests(fit2)
summary(res)

# Extract MAPT_P301l DMPs
annEpicSub <- annEpic[match(rownames(mapt.mvals),annEpic$Name), c(1:4, 12:19, 24:ncol(annEpic))]
dmps_mapt <- topTable(fit2, num=Inf, coef=1, genelist = annEpicSub)
write.table(dmps_mapt, "DMPs_maptP301L.ndc_funnNorm.txt")
#=============================================#


#=== Differentially methylated regions analysis ===#
annot.mapt = cpg.annotate(object = mapt.mvals, datatype = "array", what="M", analysis.type="differential",
                     design=design.matrix, contrasts=T, cont.matrix=contrast.matrix.mapt,
                     coef="P301L-ndc", arraytype="EPIC")

# Find DMRs for MAPT-NDC
dmrs_mapt <- dmrcate(annot.mapt, lambda=1000, C=2)

# convert to annotated genomic regions
data(dmrcatedata)
results.ranges <- extractRanges(dmrs_mapt, genome = "hg19")

# set up the grouping variables and colours
groups <- pal[1:length(unique(mapt.mut))]
names(groups) <- levels(factor(mapt.mut))
cols <- groups[as.character(factor(mapt.mut))]
samps <- 1:length(mapt.mut)



# draw the plot for the top DMR
par(mfrow=c(1,1))
DMR.plot(ranges=results.ranges, dmr=2, CpGs=betaVals, phen.col=cols, what = "Beta",
         arraytype = "EPIC", pch=16, toscale=TRUE, plotmedians=TRUE, 
         genome="hg19", samps=samps)


#==================================================#


