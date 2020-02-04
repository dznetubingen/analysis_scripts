##############################################################################
# Collect Pipeline Results from iPSC-neural smRNA-seq datasets processed with OASIS
##############################################################################
library(stringr)
setwd("~/rimod/smallRNA/iPSC/")

out <- "~/rimod/smallRNA/frontal/processed_data/Rimod_smRNAseq_2019/data/counts/"

# Read extra data (2019, frmo Tübingen)
files2 <- list.files(out, pattern='*allspeciesCounts*')
# create initial data frame
df <- read.table(paste0(out, files2[1]), sep="\t")
colnames(df) <- c("gene", files2[1])

# iterate over remaining files
for (i in 2:length(files2)) {
  tmp <- read.table(paste0(out, files2[i]), sep="\t")
  colnames(tmp) <- c("gene", files2[i])
  
  # merge
  df <- merge(df, tmp, by='gene')
  
}
rownames(df) <- df$gene
df <- df[, -1]
tub_df <- df


# Now Split in iPSC data and Human data
human_samples <- c("final_5bp_trimmed_sample_11076_F_allspeciesCounts.txt",
                   "final_5bp_trimmed_sample_98169F_allspeciesCounts.txt", "final_5bp_trimmed_sample_10316_F_allspeciesCounts.txt", "final_5bp_trimmed_sample_11082_F_allspeciesCounts.txt",
                   "final_5bp_trimmed_sample_95231F_allspeciesCounts.txt", "final_5bp_trimmed_sample_11021_F_allspeciesCounts.txt", "final_5bp_trimmed_sample_12042_F_allspeciesCounts.txt",
                   "final_5bp_trimmed_sample_98061_F_allspeciesCounts.txt", "final_5bp_trimmed_sample_00116F_allspeciesCounts.txt")
tub_ips <- tub_df[,!colnames(tub_df) %in% human_samples]
df <- tub_ips

# subset to human miRNAs
df <- df[grepl("hsa", rownames(df)),]
# adjust colnames
colnames(df) <- gsub("final_5bp_trimmed_sample_", "", colnames(df))
colnames(df) <- gsub("_allspeciesCounts.txt", "", colnames(df))


###
# Get metadata
###
md <- read.csv("ipsc_neuron_metadata.csv", sep=",", stringsAsFactors = F)
md$sample <- gsub("-", "_", md$sample)
md$sample <- gsub(" ", "", md$sample)

# manually adjust all the names to match ...
colnames(df)[1] <- "8355_CL6"
colnames(df)[2] <- "BRECL_CL26"
colnames(df)[3] <- "BRECL_CL70"
colnames(df)[4] <- "DN19_Cl1"
colnames(df)[5] <- "DN19_Cl2"
colnames(df)[6] <- "F_FD_07Cl6"
colnames(df)[7] <- "F_FD_07Cl23"
colnames(df)[8] <- "LZ_FD_13Cl10"
colnames(df)[9] <- "GM23279"
colnames(df)[10] <- "GM23280"
colnames(df)[11] <- "K7.1"
colnames(df)[12] <- "LZ_FD_13Cl3" # not sure here whether that is actually correct
colnames(df)[13] <- "NAS6"
colnames(df)[14] <- "NAS7"
colnames(df)[15] <- "ND32955_Cl32"
colnames(df)[16] <- "ND41870"
colnames(df)[17] <- "ND55"
colnames(df)[18] <- "ND66"
colnames(df)[19] <- "ST127615Cl49"
colnames(df)[20] <- "ST127615_Cl4"


write.table(df, "iPSCNeurons_smRNAseq_counts.txt", sep="\t", quote=F, col.names = NA)

write.table(md, "iPSCNeurons_smRNAseq_metadata_formatted.txt", sep="\t", quote=F, row.names = F)






