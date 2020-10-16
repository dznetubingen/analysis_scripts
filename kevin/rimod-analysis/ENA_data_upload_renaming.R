
library(stringr)

setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/data_upload/")

# Load the master sample table
md <- read.table("RiMod_master_sample_file.txt", sep="\t", header=T)
number <- as.numeric(gsub("rimod", "", md$Sample_UID))
md <- md[order(number),]

write.table(md, "ordered_master_table.txt", sep="\t", quote=F, row.names = F)

# load id-mapping table
idmap <- read.table("RiMod_ID_mapping.txt", sep="\t", header=T, stringsAsFactors = F)


####
# smRNA Tübingen renaming
#####
# Rename smRNA-seq Tübingen data 
setwd("smrnaseq_frontal_tübingen/final_trimmed/")
files <- list.files()
nf <- files
nf <- gsub("final_5bp_trimmed_sample_", "", nf)
nf <- gsub("F.fastq.gz", "", nf)
nf <- gsub("_", "", nf)
tmp <- idmap[idmap$old_id %in% nf,]
tmp <- tmp[match(nf, tmp$old_id),]
nf <- tmp$new_id
nf <- paste("smRNAseq_tuebingen_frontal_human_", nf, ".fastq.gz", sep="")

# do the renaming
file.rename(files, nf)

#=== end renaming smRNA Tübingen ===#

####
# smRNA Göttingen renaming
#####
setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/data_upload/smrnaseq_frontal_göttingen/")
files <- list.files()
samples <- files
samples <- gsub("RNAomeTb", "", samples)
samples <- gsub("NDC", "FTD", samples)
samples <- str_split(samples, pattern="FTD", simplify = T)[,1]
# rename 103277 and 110140
samples[samples == "103277"] <- "10327"
samples[samples == "110140"] <- "11014"
all(samples %in% idmap$old_id)
tmp <- idmap[idmap$old_id %in% samples,]
tmp <- tmp[match(samples, tmp$old_id),]
nf <- tmp$new_id
nf <- paste("smRNAseq_goettingen_frontal_human_", nf, ".fastq.gz", sep="")

# do the renaming
file.rename(files, nf)

#======== end renaming smRNA Göttingen ===#




####
# Rename frontal CAGE-seq data from human post-mortem brain
####
setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/data_upload/frontal_cageseq/cageseq_fastq/")
files <- list.files()
samples <- files
samples <- gsub("_fro.fastq.gz", "", samples)
# rename A144_12
samples[samples == "A144_12"] <- "0A144"
all(samples %in% idmap$old_id)
tmp <- idmap[idmap$old_id %in% samples,]
tmp <- tmp[match(samples, tmp$old_id),]
nf <- tmp$new_id
nf <- paste("CAGEseq_frontal_", nf, ".fastq.gz", sep="")

# do the renaming
file.rename(files, nf)

#========== end renaming CAGE-seq frontal data ===========#


####
# Rename frontal 
