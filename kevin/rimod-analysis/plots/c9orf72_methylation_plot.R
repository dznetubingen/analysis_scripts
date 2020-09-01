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
md <- read.csv("/Users/kevin/dzne/rimod_package/files/FTD_Brain.csv")
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


# Create Matrix
targets$Group[targets$Group == "FTD_MAPT"] <- "FTD-MAPT"
targets$Group <- make.names(targets$Group)
G <- factor(targets$Group)
A <- targets$age
design.matrix <- model.matrix(~ 0 + G + A)
colnames(design.matrix) <- c("FTD.C9", "FTD.GRN", "FTD.MAPT", "NDC", "AGE")


## Run SVA
mod1 <- model.matrix(~ 0 + G)
mod0 <- cbind(mod1[,1])
svs <- sva(mVals, mod1, mod0)$sv
svs <- svs[, c(1,2,3)]
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
# C9orf72
dmps_c9 <- topTable(fit2, num=Inf, coef="FTD.C9-NDC", genelist = annEpicSub)


# plot the top 4 most significantly differentially methylated CpGs 
# plot C9orf72 DMPs around C9orf72 locus
tmp <- dmps_c9[dmps_c9$chr == "chr9",]
tmp <- tmp[grepl("C9orf72", tmp$GencodeBasicV12_NAME),]

par(mfrow=c(2,2))
sapply(rownames(tmp)[1:4], function(cpg){
  plotCpg(betaVals, cpg=cpg, pheno=targets$Group, ylab = "Beta values")
  print(cpg)
})
par(mfrow=c(1,1))

plotCpg(betaVals, cpg = "cg07052794", pheno=targets$Group, ylab="Beta values")

# Look into them samples
sites <- rownames(tmp)[1:4]
beta.tmp <- betaVals[sites,]
plot(beta.tmp[1,])

plot(beta.tmp[2,])
text(beta.tmp[2,])

plot(beta.tmp[3,])
text(beta.tmp[3,])

plot(beta.tmp[4,])
text(beta.tmp[4,])

targets$age
targets$Sample_Name[c(24, 48)]
#============================================#

# Make plot for figure
# make df
vals <- betaVals[sites[2],]
group <- targets$Group
df <- data.frame(Beta = betaVals[sites[2],], Group = targets$Group)
mypal <- c("#7570B3", "#db6e1a","#67e08a", "#616665")

ggplot(df, aes(x=Group, fill=Group, y=Beta)) + 
  geom_boxplot() +
  geom_point() +
  scale_fill_manual(values = mypal) +
  theme_minimal() +
  labs(x = "", y = "Beta value") +
  theme(axis.text.x = element_blank()) 

ggsave("c9orf72_methylation.png", width=3.5, height=3.5)