##############################################################################
# Collect Pipeline Results from frontal smRNA-seq datasets processed with OASIS
##############################################################################
library(stringr)
setwd("~/rimod/smallRNA/frontal/")

out1 <- "~/rimod/smallRNA/frontal/processed_data/Rimod_smRNAseq_2018/data/counts/"
out2 <- "~/rimod/smallRNA/frontal/processed_data/Rimod_smRNAseq_2019/data/counts/"



## Read old data (2018, form Göttingen)
files1 <- list.files(out1, pattern='*allspeciesCounts*')
# create initial data frame
df <- read.table(paste0(out1, files1[1]), sep="\t")
colnames(df) <- c("gene", files1[1])

# iterate over remaining files
for (i in 2:length(files1)) {
  tmp <- read.table(paste0(out1, files1[i]), sep="\t")
  colnames(tmp) <- c("gene", files1[i])
  # merge
  df <- merge(df, tmp, by='gene')
}

rownames(df) <- df$gene
df <- df[, -1]
got_df <- df

# Read extra data (2019, frmo Tübingen)
files2 <- list.files(out2, pattern='*allspeciesCounts*')
# create initial data frame
df <- read.table(paste0(out2, files2[1]), sep="\t")
colnames(df) <- c("gene", files2[1])

# iterate over remaining files
for (i in 2:length(files2)) {
  tmp <- read.table(paste0(out2, files2[i]), sep="\t")
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
                   "final_5bp_trimmed_sample_98061_F_allspeciesCounts.txt")
tub_hs <- tub_df[,colnames(tub_df) %in% human_samples]
tub_ips <- tub_df[,!colnames(tub_df) %in% human_samples]

hs_df <- merge(got_df, tub_hs, by="row.names")
rownames(hs_df) <- hs_df$Row.names
hs_df <- hs_df[, -1]



# adjust samplenames
samples <- colnames(hs_df)
samples <- gsub("RNAomeTb", "", samples)
samples <- gsub("_allspeciesCounts.txt", "", samples)
samples <- gsub("_mm_smallrna_sr_Farah", "", samples)
samples <- gsub("final_5bp_trimmed_", "", samples)
colnames(hs_df) <- samples



# now create design file with batch information included
md <- read.csv("~/rimod/smallRNA/FTD_Brain.csv")
mds <- as.character(md$SAMPLEID)
mds <- str_pad(mds, width=5, side='left', pad='0')
md$id <- mds



# Assign Disease Codes
dc <- rep("NDC", length(samples))
dc[grepl("MAPT", samples)] <- "FTD.MAPT"
dc[grepl("GRN", samples)] <- "FTD.GRN"
dc[grepl("C9", samples)] <- "FTD.C9"
# manually for new samples
for (i in 39:length(dc)){
  s <- samples[i]
  s <- gsub("sample", "", gsub("_", "", gsub("F", "", s)))
  if (s %in% md$id){
    tmp <- md[md$id == s,]
    tmp <- as.character(tmp$DISEASE.CODE)[1]
    dc[i] <- tmp
  }
}
dc[dc == 'control'] = 'NDC'
dc <- gsub("-", ".", dc)

# create batch variable
batch1 <- rep("gottingen", 38)
batch2 <- rep("tubingen", length(human_samples))
batch <- c(batch1, batch2)

sample_names <- samples

# get IDs
samples <- gsub("sample_", "", samples)
samples <- str_sub(samples, start=1, end = 5)


df <- data.frame(dc = dc, sample = sample_names, id = samples, batch = batch)

write.table(df, "rimod_human_frontal_smRNAseq_metadata.txt", sep="\t", quote=F, col.names=NA)
write.table(hs_df, "rimod_human_frontal_smRNAseq_counts.txt", sep="\t", quote=F, col.names = NA)



# iPSC samples
samples <- gsub("_allspeciesCounts.txt", "", gsub("final_5bp_trimmed_", "", colnames(tub_ips)))
colnames(tub_ips) <- samples
write.table(tub_ips, "rimod_iPSC_frontal_smRNAseq_counts.txt", sep="\t", quote=F, col.names=NA)
