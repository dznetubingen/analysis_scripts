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
library(sva)
library(stringr)

# Get annotation
annEpicObj <- getAnnotationObject(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
annEpic <- getAnnotation(annEpicObj)

###############
# Read in Data
###############
# Data directory
data.dir <- "~/rimod/Methylation/frontal_methylation_0818/"

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
md <- read.csv("~/rimod/files/FTD_Brain_corrected.csv")

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
# Apply preprocessQuantile function as it is more suited for samples that don't have globally different patterns
# expression profiles
mSetFn <- preprocessQuantile(RGset)
mSetRaw <- preprocessRaw(RGset) # for plotting

# Create density plots of the Beta values to compare raw with normalized values
png("raw_density_plot.png", width=800, height=500)
densityPlot(getBeta(mSetRaw), main ="Raw")
dev.off()
png("normalized_quantile_density_plot", width=800, height=500)
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


####
# Apply fitering based on variability
mvals_sd <- apply(mVals, 1, sd)
beta_sd <- apply(betaVals, 1, sd)
# Filter out CpGs with a betaValue SD less than 0.1
keep <- beta_sd > 0.05
mVals <- mVals[keep,]
betaVals <- betaVals[keep,]

# Apply additional filtering using DMRcate
mVals <- rmSNPandCH(mVals, rmcrosshyb = TRUE, rmXY = TRUE)
betaVals <- rmSNPandCH(betaVals, rmcrosshyb = TRUE, rmXY = TRUE) 

# Create matrix for external use
cnames <- paste(targets$Sample_Name2, targets$Group, sep="_")
mvals_ext <- mVals
colnames(mvals_ext) <- cnames
write.table(mvals_ext, "mVals_matrix_frontal_methylation.txt", quote=F, sep="\t")

bvals <- betaVals
colnames(bvals) <- cnames
write.table(bvals, "betaVals_matrix_frontal_metyhlation.txt", quote=F, sep="\t")


# Create Matrix
targets$Group[targets$Group == "FTD_MAPT"] <- "FTD-MAPT"
targets$Group <- make.names(targets$Group)
G <- factor(targets$Group)
A <- targets$age
design.matrix <- model.matrix(~ 0 + G + A)
colnames(design.matrix) <- c("FTD.C9", "FTD.GRN", "FTD.MAPT", "NDC", "AGE")

## matrix for age prediction
gender <- as.character(targets$gender)
gender[gender == "F"] <- 1
gender[gender == "M"] <- 0
horvath <- data.frame(Age = A, Female = as.numeric(gender), Tissue = rep("Brain FCTX", length(A)))
write.table(horvath, "age_prediction_annotation.csv", sep=",", row.names = F, quote=F)

## Run SVA
mod1 <- model.matrix(~ 0 + G)
mod0 <- cbind(mod1[,1])
svs <- sva(mVals, mod1, mod0)$sv
colnames(svs) <- c("SV1", "SV2", "SV3")

design.matrix <- model.matrix(~ 0 + G + svs)
colnames(design.matrix) <- c("FTD.C9", "FTD.GRN", "FTD.MAPT", "NDC", "SV1", "SV2", "SV3")

# Create contrasts matrix
conts <- c("FTD.C9-NDC", "FTD.GRN-NDC", "FTD.MAPT-NDC")
contrast.matrix <- makeContrasts(contrasts = conts, levels = design.matrix)
# Limma fitting
fit <- lmFit(mVals, design = design.matrix)
cont.fit <- contrasts.fit(fit = fit, contrasts = contrast.matrix)
fit2 <- eBayes(cont.fit)
res <- decideTests(fit2)
summary(res)

# Extract and save DMPs
keep_cols = c(1,2,3,4, 36, 37, 38, 39, 40, 41)
annEpicSub <- annEpic[match(rownames(mVals),annEpic$Name), c(1:4, 12:19, 24:ncol(annEpic))]
# MAPT
dmps_mapt <- topTable(fit2, num=Inf, coef="FTD.MAPT-NDC", genelist = annEpicSub)
#dmps_mapt <- dmps_mapt[, keep_cols]
write.table(dmps_mapt, "DMPs_mapt.ndc_quant.txt", sep="\t")
# GRN
dmps_grn <- topTable(fit2, num=Inf, coef="FTD.GRN-NDC", genelist = annEpicSub)
#dmps_grn <- dmps_grn[, keep_cols]
write.table(dmps_grn, "DMPs_grn.ndc_quant.txt", sep="\t")
# C9orf72
dmps_c9 <- topTable(fit2, num=Inf, coef="FTD.C9-NDC", genelist = annEpicSub)
#dmps_c9 <- dmps_c9[, keep_cols]
write.table(dmps_c9, "DMPs_c9orf72.ndc_quant.txt", sep="\t")


# plot the top 4 most significantly differentially methylated CpGs 
# plot C9orf72 DMPs around C9orf72 locus
tmp <- dmps_c9[dmps_c9$chr == "chr9",]
par(mfrow=c(2,2))
sapply(rownames(tmp)[1:4], function(cpg){
  plotCpg(betaVals, cpg=cpg, pheno=targets$Group, ylab = "Beta values")
})
par(mfrow=c(1,1))
plotCpg(betaVals, cpg = "cg07052794", pheno=targets$Group, ylab="Beta values")

# Save mvals

write.table(mVals, "mVals_quant.txt", sep="\t", quote=F)
#============================================#



#==== Differential Region analysis ====#

## DMR MAPT-NDC
mapt.anno <- cpg.annotate(object = mVals, datatype = 'array', what='M', analysis.type = 'differential',
                          design = design.matrix, contrasts = TRUE, cont.matrix = contrast.matrix,
                          coef="FTD.MAPT-NDC", arraytype = "EPIC")
dmrs = dmrcate(mapt.anno, lambda=1000, C=2)
dmr.mapt <- dmrs$results
write.table(dmr.mapt, "DMR_MAPTvsNDC.txt", sep="\t")

## DMR GRN-NDC
grn.anno <- cpg.annotate(object = mVals, datatype = 'array', what='M', analysis.type = 'differential',
                         design = design.matrix, contrasts = TRUE, cont.matrix = contrast.matrix,
                         coef="FTD.GRN-NDC", arraytype = "EPIC")
dmrs = dmrcate(grn.anno, lambda=1000, C=2)
dmr.grn <- dmrs$results
write.table(dmr.grn, "DMR_GRNvsNDC.txt", sep="\t")

## DMR c9orf72-NDC
c9.anno <- cpg.annotate(object = mVals, datatype = 'array', what='M', analysis.type = 'differential',
                        design = design.matrix, contrasts = TRUE, cont.matrix = contrast.matrix,
                        coef="FTD.C9-NDC", arraytype = "EPIC")
dmrs = dmrcate(c9.anno, lambda=1000, C=2)
dmr.c9 <- dmrs$results
write.table(dmr.c9, "DMR_C9vsNDC.txt", sep="\t")

#=======================================#


#==========================================#
# Extract genes near DMPs

# mapt
mapt.dmp <- dmps_mapt[dmps_mapt$adj.P.Val <= 0.05,]
mapt.genes <- mapt.dmp$GencodeBasicV12_NAME
mapt.genes <- str_split(mapt.genes, pattern=";", simplify = T)[,1]

# mapt.up
mapt.dmp.up <- mapt.genes[mapt.dmp$logFC > 0]
mapt.dmp.down <- mapt.genes[mapt.dmp$logFC < 0]
mapt.dmp.up <- mapt.dmp.up[!mapt.dmp.up == ""]
mapt.dmp.down <- mapt.dmp.down[!mapt.dmp.down == ""]

grn.dmp <- dmps_grn[dmps_grn$adj.P.Val <= 0.05,]
grn.genes <- grn.dmp$GencodeBasicV12_NAME
grn.genes <- str_split(grn.genes, pattern=";", simplify = T)[,1]

# mapt.up
grn.dmp.up <- grn.genes[grn.dmp$logFC > 0]
grn.dmp.down <- grn.genes[grn.dmp$logFC < 0]
grn.dmp.up <- grn.dmp.up[!grn.dmp.up == ""]
grn.dmp.down <- grn.dmp.down[!grn.dmp.down == ""]

# remove duplicates
mapt.dmp.up <- mapt.dmp.up[!duplicated(mapt.dmp.up)]
mapt.dmp.down <- mapt.dmp.down[!duplicated(mapt.dmp.down)]
grn.dmp.up <- grn.dmp.up[!duplicated(mapt.dmp.up)]
grn.dmp.down <- grn.dmp.down[!duplicated(grn.dmp.down)]


write.table(mapt.dmp.up, "MAPT_DMPs_UP.txt", quote=F, sep="\t", row.names=F, col.names = F)
write.table(mapt.dmp.down, "MAPT_DMPs_DOWN.txt", quote=F, sep="\t", row.names=F, col.names = F)
write.table(grn.dmp.up, "GRN_DMPs_UP.txt", quote=F, sep="\t", row.names=F, col.names = F)
write.table(grn.dmp.down, "GRN_DMPs_DOWN.txt", quote=F, sep="\t", row.names=F, col.names = F)