#######################
# C9orf72 Methylation results exploration
# This script's main purpose was to examine if there is any kind of differential methylation around
# the C9orf72 locus
###################
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

# Data directory
data.dir <- "~/rimod/Methylation/frontal_methylation_0818"
setwd(data.dir)

dmp <- read.table("DMPs_c9orf72.ndc_quant.txt", header=T)
dmp <- dmp[dmp$P.Value <= 0.01,]
dmp <- dmp[dmp$chr == "chr9", ]
