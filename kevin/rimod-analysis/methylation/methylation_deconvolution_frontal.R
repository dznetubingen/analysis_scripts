####################################
## Frontal Methylation Data Analysis
####################################

## Deconvolution analysis

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

################
# Parsing and QC
################

data.dir <- "~/rimod/Methylation/frontal_methylation_0818"
setwd(data.dir)
list.files(data.dir)
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

# Get annotation
annEpicObj <- getAnnotationObject(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
annEpic <- getAnnotation(annEpicObj)

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
group <- as.character(targets$Group)
qcReport(RGset, sampNames = targets$Sample_Name, sampGroups = group, pdf="QC_report.pdf")

# remove samples with bad detection p-values
keep <- colMeans(detP) < 0.05
RGset <- RGset[,keep]



#####
# Normalization
#####
# the result from the quantro function can be used to choose the normalization procedure

# Apply preprocessQuantile function
mSetSq <- preprocessQuantile(RGset)
#mSetRaw <- preprocessRaw(RGset)
#mSetFn <- preprocessFunnorm(RGset)

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
#

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

#####
# Differential methylation analysis
#####
# Get M-values for analysis and Beta-values for visualization
mVals <- getM(mSetSqFlt)
betaVals <- getBeta(mSetSqFlt)

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

#####
# Deconvolution Analysis 
#####
library(RefFreeEWAS)

cm <- RefFreeCellMix(Y=betaVals, K=8)
