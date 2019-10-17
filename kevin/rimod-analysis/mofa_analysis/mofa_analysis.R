# Mofa testing
library(MOFA)
library(MOFAdata)
library(stringr)
setwd("/data/kevin/mofa_files/")
# Load design
design <- read.csv("rimod_design.csv")

# Load datasets

# RNA-seq
rna <- read.table("deseq_vst_values_2019-08-12_07.58.35.txt", sep="\t", header=T, row.names = 1)
rownames(rna) <- str_split(rownames(rna), pattern = "[.]", simplify = T)[,1]

# sRNA-seq
mirna <- read.table("deseq_rLog_values_frontal_smRNA.txt", sep="\t", header=T, row.names = 1)

# Methylation
met <- read.table("mVals_matrix_frontal_methylation.txt", sep="\t", header=T)

# CAGE-seq
cage <- read.table("cageseq_rlog_frontal.txt", sep="\t", header=T, row.names=1)


# Subset design and pad sample names
design <- design[design$REGION == "frontal",]
samples <- as.character(design$SAMPLEID)
samples <- str_pad(samples, 5, side="left", pad="0")
samples <- paste("sample", samples, sep="_")
design$samples <- samples

### Adjust RNA-seq data
rna_samples <- colnames(rna)
rna_samples <- gsub("X", "", rna_samples)
rna_samples <- str_split(rna_samples, pattern="_", simplify = T)[,1]
rna_samples[rna_samples == "A144"] <- "A144_12"
rna_samples <- paste("sample", rna_samples, sep="_")
colnames(rna) <- rna_samples
# Filter by value
rna_sums <- apply(rna, 1, sum)
rna <- rna[rna_sums > 200,]

# Adjust sRNA-seq data
# first get batch list
mirna_batch <- colnames(mirna) # create vector of batches
mirna_batch[grepl("X", mirna_batch)] <- "batch1"
mirna_batch[grepl("sample", mirna_batch)] <- "batch2"

mir_samples <- colnames(mirna)
mir_samples <- gsub("X", "", mir_samples)
mir_samples <- gsub("sample_", "", mir_samples)
mir_samples <- str_sub(mir_samples, 1, 5)
mir_samples <- paste("sample", mir_samples, sep="_")
colnames(mirna) <- mir_samples
mirna <- mirna[colnames(mirna) %in% samples]
mirna_batch <- mirna_batch[colnames(mirna) %in% samples]
batch2_samples <- mir_samples[mirna_batch == 'batch2']

# Adjust methylation data
met_samples <- colnames(met)
met_samples <- gsub("X", "", met_samples)
met_samples <- str_sub(met_samples, 1, 5)
met_samples <- paste("sample", met_samples, sep="_")
colnames(met) <- met_samples
met <- met[colnames(met) %in% samples]

# Adjust cage data
cage_samples <- colnames(cage)
cage_samples <- gsub("sample_", "", cage_samples)
cage_samples <- str_split(cage_samples, pattern="_", simplify = T)[,1]
cage_samples[cage_samples == "A144"] <- "A144_12"
cage_samples <- paste("sample", cage_samples, sep="_")
colnames(cage) <- cage_samples

# Get list of all available samples
avail_samples <- union(colnames(rna), union(colnames(met), union(colnames(cage), colnames(mirna))))
test_samples <- samples[samples %in% avail_samples]

# Fill out datasets with NAs where necessary
for (s in avail_samples){
  
  # Fill out RNA-seq samples
  if (!s %in% colnames(rna)){
    rna$tmp <- rep(NA, nrow(rna))
    colnames(rna)[ncol(rna)] <- s
  }
  
  # Fill out miRNA samples
  if (!s %in% colnames(mirna)){
    mirna$tmp <- rep(NA, nrow(mirna))
    colnames(mirna)[ncol(mirna)] <- s
  }
  
  # Fill out Methylation samples
  
  if (!s %in% colnames(met)){
    met$tmp <- rep(NA, nrow(met))
    colnames(met)[ncol(met)] <- s
  }
  
  # Fill out CAGE samples
  if (!s %in% colnames(cage)){
    cage$tmp <- rep(NA, nrow(cage))
    colnames(cage)[ncol(cage)] <- s
  }
}


# Bring all matrices in the same order
rna <- rna[avail_samples]
mirna <- mirna[avail_samples]
met <- met[avail_samples]
cage <- cage[avail_samples]
data_list = list(rna, mirna, met, cage)


####
# MOFA ANALYSIS
####

# Create the MOFA object
mofa <- createMOFAobject(data_list)
viewNames(mofa) <- c("RNAseq", "smRNAseq", "Methylation", "CAGEseq")

# Define DataOptions
DataOptions <- getDefaultDataOptions()

# Define model options
ModelOptions <- getDefaultModelOptions(mofa)
ModelOptions$numFactors <- 11

# Define training options
TrainOptions <- getDefaultTrainOptions()
TrainOptions$DropFactorThreshold <- 0.01
TrainOptions$seed <- 1578
TrainOptions$maxiter <- 1000


# Prepare MOFA
mofa <- prepareMOFA(
  mofa,
  DataOptions = DataOptions,
  ModelOptions = ModelOptions,
  TrainOptions = TrainOptions
)


# Regress out batch
mirna_batch <- rep("batch1", ncol(mirna))
mirna_batch[colnames(mirna) %in% batch2_samples] <- "batch2"
mofa <- regressCovariates(mofa, views="smRNAseq", covariates = mirna_batch)

# regress out gender
#design <- design[design$samples %in% avail_samples,]
#design <- design[match(avail_samples, design$samples),]
#mofa <- regressCovariates(mofa, views = c("RNAseq", "smRNAseq", "Methylation", "CAGEseq"), covariates = design$GENDER)


# Perform training
mofa <- runMOFA(mofa)

###
# Explore MOFA results
###

# calculation of variance explained
r2 <- calculateVarianceExplained(mofa)
r2$R2Total
plotVarianceExplained(mofa)

# Factor correlation
plotFactorCor(mofa)


# Reactome enrichment analysis
data("reactomeGS")
data("MSigDB_v6.0_C2_human")

gsea.rna <- runEnrichmentAnalysis(mofa,
                              view = "RNAseq",
                              feature.sets = reactomeGS,
                              alpha = 0.01)



gsea.cage <- runEnrichmentAnalysis(mofa,
                                  view = "CAGEseq",
                                  feature.sets = reactomeGS,
                                  alpha = 0.01)


plotEnrichmentBars(gsea.rna)
plotEnrichmentBars(gsea.cage)

plotEnrichment(mofa, gsea.rna, 3, alpha = 0.01, max.pathways = 20)
plotEnrichment(mofa, gsea.cage, 2, alpha = 0.01, max.pathways = 20)

## msigDB enrichment
data("MSigDB_v6.0_C5_human")

gsea.rna.msig <- runEnrichmentAnalysis(mofa,
                                  view = "RNAseq",
                                  feature.sets = MSigDB_v6.0_C5_human,
                                  alpha = 0.01)

gsea.cage.msig <- runEnrichmentAnalysis(mofa,
                                       view = "CAGEseq",
                                       feature.sets = MSigDB_v6.0_C5_human,
                                       alpha = 0.01)

plotEnrichmentBars(gsea.rna.msig)
plotEnrichmentBars(gsea.cage.msig)

plotEnrichment(mofa, gsea.rna.msig, 4, alpha = 0.01, max.pathways = 20)
plotEnrichment(mofa, gsea.cage.msig, 4, alpha = 0.01, max.pathways = 20)


test_samples <- colnames(cage)
samples <- as.character(design$SAMPLEID)
samples <- str_pad(samples, 5, side="left", pad="0")
samples <- paste("sample", samples, sep="_")
design$samples <- samples
design <- design[match(test_samples, samples),]


plotFactorScatters(mofa,
                  factors=1:5,
                  color_by = design$DISEASE.CODE)

plotFactorScatter(mofa,
                   factors=c(2, 3),
                   color_by = design$DISEASE.CODE, dot_size = 5)



### get weights
mofa.weights <- getWeights(mofa,
                           as.data.frame = T,
                           view="smRNAseq")

lf1 <- mofa.weights[mofa.weights$factor == "LF1",]
lf1 <- lf1[order(lf1$value, decreasing = T),]



# testing
plotFactorCor(mofa)



