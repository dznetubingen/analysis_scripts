##############
# Create Master Sample File for David Craig, RiMod-FTD database
##############
library(stringr)
setwd("~/rimod/rimod_resource_data/")

# format md file
md <- read.csv("../files/FTD_Brain_corrected.csv", stringsAsFactors = F)
md$sample <- str_pad(md$SAMPLEID, width=5, side="left", pad="0")
md$sample[md$sample == "A144_12"] <- "0A144"

# load mapping file
rimod <- read.table("RiMod_ID_mapping.txt", header=T, stringsAsFactors = F)

all(md$sample %in% rimod$old_id)
md <- merge(md, rimod, by.x="sample", by.y="old_id")

##
# make RNA-seq data file
# load RNA-seq metadata and merge with old_id
rna <- read.table("~/rimod/RNAseq/rimod_frontal_rnaseq_metadata_v2.txt", header=T, sep="\t")
rna <- merge(rna, rimod, by.x="SampleID", by.y="old_id")

rna_data <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2020-05-04_15.45.57/deseq_normalized_counts_2020-05-04_15.45.57.txt")
rownames(rna_data) <- str_split(rownames(rna_data), pattern="[.]", simplify = T)[,1]
colnames(rna_data) <- gsub("X", "", colnames(rna_data))
colnames(rna_data) <- str_split(colnames(rna_data), pattern="_", simplify = T)[,1]
colnames(rna_data) <- str_pad(colnames(rna_data), width=5, side="left", pad="0")
rna <- rna[match(colnames(rna_data), rna$SampleID),]
rna$SampleID == colnames(rna_data)
# TRUE
colnames(rna_data) <- rna$new_id

write.table(rna_data, "RNAseq.frontal.version0.tsv", sep="\t", quote=F, col.names = NA)

#====== end RNA-seq ==========#

####
# make miRNA data file
# load smallRNA-seq data
mirna <- read.table("~/rimod/smallRNA/frontal/rimod_human_frontal_smRNAseq_metadata.txt", stringsAsFactors = F)
mirna$sample_id <- str_pad(mirna$id, width=5, side="left", pad="0")
all(mirna$sample_id %in% rimod$old_id)

mirna <- merge(mirna, rimod, by.x="sample_id", by.y="old_id")

# load data
srna <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_normalized_counts_temporal_smRNA.txt")
colnames(srna) <- gsub("X", "", colnames(srna))
all(colnames(srna) %in% mirna$sample)
mirna <- mirna[match(colnames(srna), mirna$sample),]
all(mirna$sample == colnames(srna))
# TRUE
colnames(srna) <- mirna$new_id

write.table(srna, "smRNAseq.frontal.version0.tsv", sep="\t", quote=F, col.names = NA)


#=========== end smRNA-seq ==============#

####
# Make methylation data file
# load methylation data
met <- read.table("~/rimod/Methylation/frontal_methylation_0818/betaVals_matrix_frontal_metyhlation.txt")
metcols <- gsub("X", "", colnames(met))
metcols <- str_split(metcols, pattern="_", simplify = T)[,1]
metcols[metcols == "14412"] <- "0A144"
all(metcols %in% rimod$old_id) 
colnames(met) <- metcols

tmp <- rimod[rimod$old_id %in% metcols,]
tmp <- tmp[match(colnames(met), tmp$old_id),]
all(tmp$old_id == colnames(met))
# TRUE
colnames(met) <- tmp$new_id

write.table(met, "Methylation.frontal.version0.tsv", sep="\t", quote=F, col.names = NA)

#=========== end methylation ==========#

####
# CAGE-seq bed file creation
####
fro.files <- list.files("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/CAGEseq/bed_files/", full.names = T)
fro.files <- fro.files[grepl("_fro", fro.files)]

# CPM normalization
cpm_norm <- function(x){
  (x / sum(x)) * 1000000
  return(x)
}

for (i in 1:length(fro.files)) {
  print(i)
  f = fro.files[i]
  # get the correct name
  fname <- str_split(f, pattern="//", simplify = T)[,2]
  fname <- str_split(fname, pattern="_", simplify = T)[,2]
  
  if (!fname %in% md$sample){
    print(fname)
  }
  
  tmp <- md[md$sample == fname,]
  tmp <- tmp[!duplicated(tmp$sample),]
  new_fname <- tmp$new_id
  
  f = read.table(f)
  # perform CPM normalization
  f$V5 <- cpm_norm(f$V5)
  
  write.table(f, paste0("cageseq/CAGEseq.", new_fname, ".bed"), row.names = F, col.names = F, quote=F, sep="\t")
}



#========== end CAGE-seq ==========#


##
# Create the super file
df <- md[!duplicated(md$new_id),]

# Check if that covers all them samples
all(rna$new_id %in% df$new_id)
all(mirna$new_id %in% df$new_id)
all(metcols %in% df$sample)
# YES

# cleaning
mutation <- df$GENE
mutation[is.na(mutation)] <- "N/A"

pmd <- df$PMD.MIN.
pmd[is.na(pmd)] <- "N/A"
pmd[pmd == "#VALUE!"] <- "N/A"

ph <- df$PH
ph[is.na(ph)] <- "N/A"

# create dataframe with only necessary IDs
df <- data.frame(Sample_UID = df$new_id,
                 age = df$AGE,
                 gender = df$GENDER,
                 disease_code = make.names(df$DISEASE.CODE),
                 mutation = mutation,
                 pmd_min = pmd,
                 ph = ph)

write.table(df, "RiMod_master_sample_file.txt", sep="\t", row.names = F, quote=F)
