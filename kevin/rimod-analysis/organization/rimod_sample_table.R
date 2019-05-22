################################################################
# Create a sample table for rimod with informaiton about all 
# samples and the data available for them
################################################################
library(stringr)
setwd("~/rimod/overview/")

# Read Brain FTD file
md <- read.csv("~/rimod/files/FTD_Brain.csv")
cage_fro <- md[md$REGION == 'frontal',]
cage_tem <- md[md$REGION == 'temporal',]
# cage frontal samples
cage_fro <- as.character(cage_fro$SAMPLEID)
samples <- as.character(sapply(cage_fro, function(x){strsplit(x, split="_")[[1]][[1]]}))
cage_fro = str_pad(cage_fro, 5, side='left', pad="0")

# cage temporal samples
cage_tem <- as.character(cage_tem$SAMPLEID)
samples <- as.character(sapply(cage_tem, function(x){strsplit(x, split="_")[[1]][[1]]}))
cage_tem = str_pad(cage_tem, 5, side='left', pad="0")

# collect all the sample IDs
samples <- as.character(md$SAMPLEID)
samples <- as.character(sapply(samples, function(x){strsplit(x, split="_")[[1]][[1]]}))
samples = str_pad(samples, 5, side='left', pad="0")
samples <- unique(samples)

# NOTE: all samples are split by "_", this removes a part of some IDs, however doesn't affect their uniqueness
# e.g. 'A144_12' becomes '0A144' because of padding. This makes it easier to keep IDs consistent throughout samples

#####
# Frontal RNA-seq samples
#####
rna <- read.table("~/rimod/RNAseq/Count_matrix.txt", sep="\t", header=T, nrows = 5)
rna_samples <- colnames(rna)
rna_samples <- as.character(sapply(rna_samples, function(x){strsplit(x, split="_")[[1]][[1]]}))
rna_samples <- str_pad(gsub("X", "", rna_samples), 5, side='left', pad="0")


####
# Frontal sRNA-seq samples
####
# first batch of samples from Göttingen
fro_srna <- list.files("~/rimod/smallRNA/frontal/fastq/")
fro_srna <- as.character(sapply(fro_srna, function(x){strsplit(x, split="_")[[1]][[1]]}))
fro_srna <- gsub("RNAomeTb", "", fro_srna) # remove weird string at begonning
fro_srna <- as.character(sapply(fro_srna, function(x){strsplit(x, split="FTD")[[1]][[1]]})) # fix FTD sample names
fro_srna <- as.character(sapply(fro_srna, function(x){strsplit(x, split="NDC")[[1]][[1]]})) # fix NDC sample names
fro_srna <- str_pad(fro_srna, width=5, side='left', pad='0')

# Tübingen batch
tue_srna <- list.files("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/sRNA_frontal/processed_fastq/")
tue_srna <- gsub(".fastq.gz", "", as.character(sapply(tue_srna, function(x){strsplit(x, split="_sample_")[[1]][[2]]})))
# remove the cell line samples
for (char in c("N", "S", "BRE", "bre", "ST", "LZ", "K", "GM", "gm", "FD")){
  tue_srna <- tue_srna[!grepl(char, tue_srna)]
}
tue_srna <- str_pad(gsub("_", "", gsub("F", "", tue_srna)), width=5, side='left', pad='0')

# merge göttingen and tübingen frontal sRNA
fro_srna <- c(fro_srna, tue_srna)

####
# Temporal sRNA-seq samples
####
tem_srna <- list.files("~/rimod/smallRNA/temporal/fastq")
tem_srna <- as.character(sapply(tem_srna, function(x){strsplit(x, split="_")[[1]][[1]]}))
tem_srna <- str_pad(gsub("temporal", "", tem_srna), width = 5, side='left', pad='0')

####
# Frontal Methylation samples
####
fro_met <- read.table("~/rimod/Methylation/frontal_methylation_0818/mVals_matrix_frontal_methylation.txt", sep="\t", nrows = 5, header=T)
fro_met <- gsub("X", "", colnames(fro_met))
fro_met <- str_pad(as.character(sapply(fro_met, function(x){strsplit(x, split="_")[[1]][[1]]})), width=5, side='left', pad='0')


####
# (Frontal?) Proteomics samples
####
prot <- read.table("~/rimod/Proteomics/RIMOD_data_2017-04-27/metadata_2017-04-27.tsv", sep="\t", header=T)
prot_fro <- prot[prot$Brain.area == "frontal",]
prot_tem <- prot[prot$Brain.area == 'temporal',]
prot_fro <- str_pad(prot_fro$Sample.ID, width=5, side='left', pad='0')
prot_tem <- str_pad(prot_tem$Sample.ID, width=5, side='left', pad='0')



################################
# Merge Everything to a Table  #
################################
rimod <- data.frame(samples)

# merge cage frontal samples
tmp <- samples
tmp[tmp %in% cage_fro] <- "yes"
tmp[!tmp=='yes'] <- 'no'
rimod$CAGE_fro <- tmp

# merge cage temporal samples
tmp <- samples
tmp[tmp %in% cage_tem] <- "yes"
tmp[!tmp=='yes'] <- 'no'
rimod$CAGE_tem <- tmp

# merge RNA-seq
tmp <- samples
tmp[tmp %in% rna_samples] <- "yes"
tmp[!tmp=='yes'] <- "no"
rimod$RNAseq_fro <- tmp

# merge frontal sRNA
tmp <- samples
tmp[tmp %in% fro_srna] <- "yes"
tmp[!tmp=='yes'] <- 'no'
rimod$sRNAseq_fro <- tmp

# merge temporal sRNA
tmp <- samples
tmp[tmp %in% tem_srna] <- "yes"
tmp[!tmp=='yes'] <- 'no'
rimod$sRNAseq_tem <- tmp

# merge methylation
tmp <- samples
tmp[tmp %in% fro_met] <- "yes"
tmp[!tmp=='yes'] <- 'no'
rimod$Methyl_fro <- tmp

# merge proteomics fro
tmp <- samples
tmp[tmp %in% prot_fro] <- "yes"
tmp[!tmp=='yes'] <- 'no'
rimod$Prot_fro <- tmp

# merge proteoimcs tem
tmp <- samples
tmp[tmp %in% prot_tem] <- "yes"
tmp[!tmp=='yes'] <- 'no'
rimod$Prot_tem <- tmp

write.table(rimod, "~/rimod/overview/rimod_sample_table.txt", sep="\t", quote=F, row.names = F)




