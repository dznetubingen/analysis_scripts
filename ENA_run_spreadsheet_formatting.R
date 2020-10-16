setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/data_upload/")

md <- read.table("data-2268886247575672299.tsv", sep="\t",  header=T)


df <- read.table("smrnaseq_frontal_tübingen/experiment_single_fastq_spreadsheet_template.tsv", sep="\t", header=T)


#####
# frontal smRNA-seq TÜBINGEN upload
######
# load files
files <- list.files("smrnaseq_frontal_tübingen/final_trimmed/")
fq <- files[!grepl("md5", files)]
samples <- gsub("smRNAseq_tuebingen_frontal_human_", "", gsub(".fastq.gz", "", fq))

# fill out common forms
library_strategy <- rep("miRNA-Seq", length(fq))
library_selection <- rep("other", length(fq))
instrument <- rep("Illumina NextSeq 550", length(fq))
library_name <- rep("NEXTFLEX Small RNA-seq Kit v3", length(fq))
source <- rep("TRANSCRIPTOMIC", length(fq))
empty <- rep("", length(fq))

df <- data.frame(sample_alias = samples,
                 instrument_model = instrument,
                 library_name = library_name,
                 library_source = source,
                 library_selection = library_selection,
                 library_strategy = library_strategy,
                 design_description = empty,
                 library_construction_protocol = empty,
                 file_name = fq)

write.table(df, "smrnaseq_tuebingen_spreadsheet.tsv", sep="\t", quote=F, row.names = F)

#===== end smRNA-seq Tübingen ========#

######
# frontal smRNA-seq GÖTTINGEN upload
######
# load files
files <- list.files("smrnaseq_frontal_göttingen/")
fq <- files[!grepl("md5", files)]
samples <- gsub("smRNAseq_goettingen_frontal_human_", "", gsub(".fastq.gz", "", fq))


library_strategy <- rep("miRNA-Seq", length(fq))
library_selection <- rep("other", length(fq))
instrument <- rep("Illumina NextSeq 550", length(fq))
library_name <- rep("NEBNext Small RNA Library Prep Set for Illumina ", length(fq))
source <- rep("TRANSCRIPTOMIC", length(fq))
empty <- rep("", length(fq))

df <- data.frame(sample_alias = samples,
                 instrument_model = instrument,
                 library_name = library_name,
                 library_source = source,
                 library_selection = library_selection,
                 library_strategy = library_strategy,
                 design_description = empty,
                 library_construction_protocol = empty,
                 file_name = fq)

write.table(df, "smrnaseq_goettingen_spreadsheet.tsv", sep="\t", quote=F, row.names = F)
#======== end smRNA-seq göttingen =======#

#########
# frontal CAGE-seq 
#########
files <- list.files("frontal_cageseq/cageseq_fastq/") 
fq <- files[!grepl("md5", files)]
samples <- gsub("CAGEseq_frontal_", "", gsub(".fastq.gz", "", fq))


library_strategy <- rep("OTHER", length(fq))
library_selection <- rep("CAGE", length(fq))
instrument <- rep("Illumina HiSeq 2500", length(fq))
library_name <- rep("Deep-CAGEseq (Takahashi et al., Nature Protocol, 2012)", length(fq))
source <- rep("TRANSCRIPTOMIC", length(fq))
empty <- rep("", length(fq))

df <- data.frame(sample_alias = samples,
                 instrument_model = instrument,
                 library_name = library_name,
                 library_source = source,
                 library_selection = library_selection,
                 library_strategy = library_strategy,
                 design_description = empty,
                 library_construction_protocol = empty,
                 file_name = fq)
write.table(df, "cageseq_frontal_spreadsheet.tsv", sep="\t", quote=F, row.names = F)

#=== end CAGE-seq spreadsheet formatting ===#


#######
# miRNA mimic experiment iPSC-neurons
#######
files <- list.files("neurons/fastq/") 
fq <- files[!grepl("md5", files)]
samples <- gsub(".fastq.gz", "", fq)

library_strategy <- rep("RNA-Seq", length(fq))
library_selection <- rep("other", length(fq))
instrument <- rep("Illumina NextSeq 550", length(fq))
library_name <- rep("TrueSeq Standard mRNA Library Prep", length(fq))
source <- rep("TRANSCRIPTOMIC", length(fq))
empty <- rep("", length(fq))

df <- data.frame(sample_alias = samples,
                 instrument_model = instrument,
                 library_name = library_name,
                 library_source = source,
                 library_selection = library_selection,
                 library_strategy = library_strategy,
                 design_description = empty,
                 library_construction_protocol = empty,
                 file_name = fq)


write.table(df, "iPSC_neurons_mirna_experiment_spreadsheet.tsv", sep="\t", quote=F, row.names = F)

#=== end iPSC neurons ===#


#######
# miRNA mimic experiment iPSC-microglia
#######
files <- list.files("microglia/Fastq/") 
fq <- files[!grepl("md5", files)]
samples <- gsub(".fastq.gz", "", fq)

library_strategy <- rep("RNA-Seq", length(fq))
library_selection <- rep("other", length(fq))
instrument <- rep("Illumina NextSeq 550", length(fq))
library_name <- rep("TrueSeq Standard mRNA Library Prep", length(fq))
source <- rep("TRANSCRIPTOMIC", length(fq))
empty <- rep("", length(fq))

df <- data.frame(sample_alias = samples,
                 instrument_model = instrument,
                 library_name = library_name,
                 library_source = source,
                 library_selection = library_selection,
                 library_strategy = library_strategy,
                 design_description = empty,
                 library_construction_protocol = empty,
                 file_name = fq)


write.table(df, "iPSC_microglia_mirna_experiment_spreadsheet.tsv", sep="\t", quote=F, row.names = F)

#=== end iPSC microglia ===#