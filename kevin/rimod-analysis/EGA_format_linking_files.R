# Link files and samples for EGA
library(stringr)
setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/data_upload/EGA_data_upload/")


# smrnaseq göttingen
files = list.files("smrnaseq_goettingen_encrypted/")
files <- files[!grepl("md5", files)]
samples <- gsub("smRNAseq_goettingen_frontal_human_", "", gsub(".fastq.gz.gpg", "", files))
files <- paste0("smrnaseq_goettingen_encrypted/", files)
df <- data.frame("Sample alias" = samples,
                 "Fastq File" = files,
                 "Checksum" = rep("", length(files)),
                 "Unencrypted checksum" = rep("", length(files)), 
                 check.names = F)

write.table(df, "goettingen_linking_file.csv", sep=",", row.names = F, quote=F)
dim(df)


# cageseq
files <- list.files("cageseq_encrypted/")
files <- files[!grepl("md5", files)]
samples <- gsub("CAGEseq_frontal_", "", gsub(".fastq.gz.gpg", "", files))
files <- paste0("smrnaseq_goettingen_encrypted/", files)
df <- data.frame("Sample alias" = samples,
                 "Fastq File" = files,
                 "Checksum" = rep("", length(files)),
                 "Unencrypted checksum" = rep("", length(files)), 
                 check.names = F)
dim(df)
write.table(df, "cageseq_linking_file.csv", sep=",", row.names = F, quote=F)


# rnaseq
files <- list.files("rnaseq_encrypted//")
files <- files[!grepl("md5", files)]
files1 <- files[grepl("_1.f", files)]
files2 <- files[grepl("_2.f", files)]
samples1 <- files1
samples2 <- files2
samples1 <- str_split(samples1, pattern="_", simplify = T)[,1]
samples2 <- str_split(samples2, pattern="_", simplify = T)[,1]
files1 <- paste0("rnaseq_encrypted/", files1)
files2 <- paste0("rnaseq_encrypted/", files2)

all(samples1 == samples2)

df <- data.frame("Sample alias" = samples1,
                 "First Fastq File" = files1,
                 "First Checksum" = rep("", length(files1)),
                 "First Unencrypted checksum" = rep("", length(files1)),
                 "Second Fastq File" = files2,
                 "Second Checksum" = rep("", length(files1)),
                 "Second Unencrypted checksum" = rep("", length(files1)),
                 check.names = F)
dim(df)
write.table(df, "rnaseq_linking_file.csv", sep=",", row.names = F, quote=F)


# microglia
files <- list.files("microglia_encrypted/")
files <- files[!grepl("md5", files)]
samples <- gsub("_all.fastq.gz.gpg", "", files)
files <- paste0("microglia_encrypted/", files)

df <- data.frame("Sample alias" = samples,
                 "Fastq File" = files,
                 "Checksum" = rep("", length(files)),
                 "Unencrypted checksum" = rep("", length(files)), 
                 check.names = F)
dim(df)
write.table(df, "microglia_linking_file.csv", sep=",", row.names = F, quote=F)

# neurons
files <- list.files("neurons_encrypted/")
files <- files[!grepl("md5", files)]
samples <- gsub("_all.fastq.gz.gpg", "", files)
files <- paste0("neurons_encrypted/", files)
df <- data.frame("Sample alias" = samples,
                 "Fastq File" = files,
                 "Checksum" = rep("", length(files)),
                 "Unencrypted checksum" = rep("", length(files)), 
                 check.names = F)
dim(df)
write.table(df, "neurons_linking_file.csv", sep=",", row.names = F, quote=F)
