# Mofa testing
library(MOFA)
library(stringr)
setwd("~/rimod/integrative_analysis/MOFA/")
# Load design
design <- read.csv("rimod_design.csv")

# Load datasets
# Proteomics
prot <- read.table("frontal/proteomics_swath_data_2017-04-27.tsv", sep="\t", header = T, row.names = 2)
prot <- prot[,-1]
# RNA-seq
rna <- read.table("frontal/RiMOD_RNAseq_frontal_CPM.matrix.txt", sep="\t", header=T)
# sRNA-seq
mirna <- read.table("frontal/miRNA_counts_merged_flt5meancounts_250218.txt", sep="\t", header=T, row.names = 1)

# Subset design and pad sample names
design <- design[design$REGION == "frontal",]
samples <- as.character(design$SAMPLEID)
samples <- str_pad(samples, 5, side="left", pad="0")
samples <- paste("sample", samples, sep="_")

### Adjust RNA-seq data
rna_samples <- colnames(rna)
rna_samples <- gsub("X", "", rna_samples)
rna_samples <- as.character(sapply(rna_samples, FUN = function(x){strsplit(x, split="_")[[1]][[1]]}))
rna_samples[rna_samples == "A144"] <- "A144_12"
rna_samples <- paste("sample", rna_samples, sep="_")
colnames(rna) <- rna_samples
# Filter by value
rna_sums <- apply(rna, 1, sum)
rna <- rna[rna_sums > 100,]



# Adjust sRNA-seq data
mir_samples <- colnames(mirna)
mir_samples <- gsub("RNAomeTb", "", mir_samples)
mir_samples <- substr(mir_samples, 1, 5)
mir_samples <- paste("sample", mir_samples, sep="_")
colnames(mirna) <- mir_samples
mirna <- mirna[colnames(mirna) %in% samples]

# Get list of all available samples
avail_samples <- union(colnames(rna), colnames(mirna))
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
}

# Bring all matrices in the same order
rna <- rna[avail_samples]
mirna <- mirna[avail_samples]
data_list = list(rna, mirna)

# Create the MOFA object
mofa <- createMOFAobject(data_list)
viewNames(mofa) <- c("RNAseq", "smRNAseq")
#plotTilesData(mofa)

# Define DataOptions
DataOptions <- getDefaultDataOptions()

# Define model options
ModelOptions <- getDefaultModelOptions(mofa)

# Define training options
TrainOptions <- getDefaultTrainOptions()
TrainOptions$DropFactorThreshold <- 0.01
TrainOptions$seed <- 1578

# Prepare MOFA
mofa <- prepareMOFA(
  mofa,
  DataOptions = DataOptions,
  ModelOptions = ModelOptions,
  TrainOptions = TrainOptions
)

# Perform training
mofa <- runMOFA(mofa)


# calculation of variance explained
r2 <- calculateVarianceExplained(mofa)
r2$R2Total
plotVarianceExplained(mofa)