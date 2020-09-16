#####
# Make mapping from Brainbank IDs to new IDs
# RiMod-FTD
#####
library(stringr)
setwd("~/rimod/")

cage <- read.csv("files/FTD_Brain_corrected.csv")
rna <- read.table("RNAseq/rimod_frontal_rnaseq_metadata_v2.txt", sep="\t", header=T)
met <- read.table("Methylation/frontal_methylation_0818/FTD_methylation_July2018.csv", sep=",", header=T)
mir <- read.table("smallRNA/frontal/rimod_human_frontal_smRNAseq_metadata.txt", sep="\t", header=T)

cage <- as.character(cage$SAMPLEID)
rna <- as.character(rna$SampleID)
met <- as.character(met$X)
mir <- as.character(mir$id)

# padding: add 0s to the left 
cage <- str_pad(cage, width=5, side="left", pad="0")
rna <- str_pad(rna, width=5, side="left", pad="0")
met <- str_pad(met, width=5, side="left", pad="0")
mir <- str_pad(mir, width=5, side="left", pad="0")

# adjust the A144/12 sample, set it to '0A144' (has length of 5)
cage[grepl("A144", cage)] <- "0A144"
rna[grepl("A144", rna)] <- "0A144"
met[grepl("A144", met)] <- "0A144"
mir[grepl("A144", mir)] <- "0A144"

# Remove duplicates
cage <- cage[!duplicated(cage)]
rna <- rna[!duplicated(rna)]
met <- met[!duplicated(met)]
mir <- mir[!duplicated(mir)]

# create the union 
samples <- union(cage, union(rna, union(met, mir)))

# make sample matching 
new_names <- paste("rimod", c(1:length(samples)), sep="")

df <- data.frame(old_id = samples, new_id = new_names)

write.table(df, "RiMod_ID_mapping.txt", sep="\t", quote=F, row.names = F)
