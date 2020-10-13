#############################################
# Create BED file methylation object for the Rimod
# Data resource
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
# remove samples with bad detection p-values
keep <- colMeans(detP) < 0.01
RGset <- RGset[,keep]
targets <- targets[keep,]
#=====================================================#



#=== Normalization ===#
# Apply preprocessQuantile function as it is more suited for samples that don't have globally different patterns
# expression profiles
mSetFn <- preprocessQuantile(RGset)



### NO FILITERING ###
library(stringr)
# Get them BETA vals
betaVals <- getBeta(mSetFn)


# Create matrix for external use
cnames <- paste(targets$Sample_Name2, targets$Group, sep="_")
bvals <- betaVals
colnames(bvals) <- cnames
write.table(bvals, "betaVals_funnorm_unfiltered_matrix_frontal_metyhlation.txt", quote=F, sep="\t")

# Extract and save DMPs
keep_cols = c(1,2,3,4, 36, 37, 38, 39, 40, 41)
annEpicSub <- annEpic[match(rownames(betaVals),annEpic$Name), c(1:4, 12:19, 24:ncol(annEpic))]

# get mapping and map to new IDs
idmap <- read.table("~/rimod/RiMod_ID_mapping.txt",  header=T)
colnames(bvals) <- str_split(colnames(bvals), pattern="_", simplify = T)[,1]
colnames(bvals)[colnames(bvals) == "14412"] <- "0A144"
all(colnames(bvals) %in% idmap$old_id)

idmap <- idmap[idmap$old_id %in% colnames(bvals),]
idmap <- idmap[match(colnames(bvals), idmap$old_id),]
all(idmap$old_id == colnames(bvals))
colnames(bvals) <- idmap$new_id

bvals <- merge(bvals, annEpicSub, by.x="row.names", by.y="row.names")
####
# Create the BED files
####

setwd("methylation_bed_files/")
samples <- idmap$new_id

for (i in 1:length(samples)) {
  sample <- samples[i]
  df <- data.frame(betaValue = as.numeric(bvals[,colnames(bvals) == sample]))
  df$start = bvals$pos
  df$end <- bvals$pos + 1
  df$strand = bvals$strand
  df$name <- rep(sample, nrow(df))
  df$chr <- bvals$chr
  df <- df[, c(6, 2, 3, 5, 1, 4)]
  
  write.table(df, paste("methylation", sample, "frontal", "bed", sep="."), col.names = F, row.names = F, quote=F, sep="\t")
}



