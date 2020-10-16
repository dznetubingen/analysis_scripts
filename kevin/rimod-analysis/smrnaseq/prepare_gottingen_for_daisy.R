library(stringr)
setwd("~/rimod/smallRNA/daisy/")

# format frontal data
frontal <- read.table("rimod_human_frontal_smRNAseq_metadata.txt")
frontal$id <- str_pad(frontal$id, pad="0", side="left", width=5)
colnames(frontal) <- c("Disease Code", "Name", "SampleID", "Batch", "Age", "Gender", "PMD", "PH")
frontal <- frontal[, -4]

# change the names to fit with the fastq files
fq <- list.files("frontal/raw_data/")
fnames <- frontal$Name
fnames <- paste("RNAomeTb", fnames, sep="")
f <- str_split(fnames, pattern="GFM", simplify = T)
fnames <- f[,1]
fnames <- paste(fnames, "GFM_mm_smallrna_sr_Farah", f[,2], sep="")

test <- paste(fnames, ".fastq.gz", sep="")
test %in% fq
fq %in% test

frontal$Name <- fnames

# change gender
gender <- as.character(frontal$Gender)
gender[gender == "F"] <- "Female"
gender[gender == "M"] <- "Male"
frontal$Gender <- gender
write.table(frontal, "smrnaseq_frontal_md_daisy.txt", sep="\t", quote=F, row.names = F)


# temporal data
fq <- list.files("temporal/smrnaseq_gottingen_temporal/")
fq.df <- as.data.frame(str_split(fq, pattern="temporal", simplify = T))

md <- read.table("rimod_human_temporal_smRNAseq_metadata.txt", sep="\t", stringsAsFactors = F, header=T, row.names = 1)
md$sample_id <- str_pad(md$sample_id, pad="0", width=5, side="left")
md <- rbind(md, data.frame(sample_id = "11014", age = 58, gender = "M", dc = "FTD-GRN"))

fq.df <- fq.df[match(md$sample_id, fq.df$V1),]

md$Name <- paste(fq.df$V1, "temporal", fq.df$V2, sep="")
#md$Name <- gsub(".fastq.gz", "", md$Name)
colnames(md) <- c("SampleID", "Age", "Gender", "Disease Code", "Name")

gender <- md$Gender
gender[gender == "M"] <- "Male"
gender[gender == "F"] <- "Female"
md$Gender <- gender

write.table(md, "smrnaseq_temporal_md_daisy.txt", sep="\t", quote=F, row.names = F)
