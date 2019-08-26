##############################################################################
# Collect Pipeline Results from temporal smRNA-seq datasets processed with OASIS
##############################################################################
library(stringr)
setwd("~/rimod/smallRNA/temporal/")

out <- "~/rimod/smallRNA/temporal/processed/rimod_temporal_sRNA/data/counts/"

## Read data
files <- list.files(out, pattern='*allspeciesCounts*')

# create initial data frame
df <- read.table(paste0(out, files[1]), sep="\t")
colnames(df) <- c("gene", files[1])

# iterate over remaining files
for (i in 2:length(files)) {
  tmp <- read.table(paste0(out, files[i]), sep="\t")
  colnames(tmp) <- c("gene", files[i])
  # merge
  df <- merge(df, tmp, by='gene')
}

rownames(df) <- df$gene
df <- df[, -1]

# Remove non-human mirnas
df <- df[grepl("hsa", rownames(df)),]




# adjust samplenames
samples <- colnames(df)
samples <- gsub("hs_smallreb_sr_sakibs_", "", samples)
samples <- gsub("_allspeciesCounts.txt", "", samples)
colnames(df) <- samples


# now create design file with batch information included
md <- read.csv("~/rimod/smallRNA/FTD_Brain.csv")
mds <- as.character(md$SAMPLEID)
mds <- str_pad(mds, width=5, side='left', pad='0')
md$id <- mds

# grep samples IDs
samples <- colnames(df)
samples <- str_split(samples, pattern="temporal", simplify = T)[,1]

# subset dataframes
keep <- samples %in% md$id
samples <- samples[keep]
df <- df[,keep]
md <- md[md$id %in% samples,]
md <- md[!duplicated(md$id),]

design <- data.frame(sample_id = md$id, age = md$AGE, gender=md$GENDER, dc = md$DISEASE.CODE)

# save data
write.table(design, "rimod_human_temporal_smRNAseq_metadata.txt", sep="\t", quote=F, col.names = NA)
write.table(df, "rimod_human_temporal_smRNAseq_counts.txt", sep="\t", quote=F, col.names = NA)

