#################
# RiMod Frontal Neuron ChIP-seq analysis
# Using DiffBind
#################
library(DiffBind)
library(stringr)
setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/chipseq/frontal_sakib/")


##
# Create Metadata file
##
chip.md <- read.csv("../chipseq_frontal_neuron_md.txt")
chip.md$Human_ID <- str_pad(chip.md$Human_ID, width = 5, side = "left", pad = "0")
# rimod overall sample sheet
md <- read.csv("~/rimod/files/FTD_Brain.csv")
md$SAMPLEID <- str_pad(md$SAMPLEID, width = 5, side = "left", pad = "0")
md <- md[md$SAMPLEID %in% chip.md$Human_ID,]
md <- md[md$REGION == "frontal",]
md <- data.frame(sample=md$SAMPLEID, age=md$AGE, sex=md$GENDER, mutation=md$GENE)

# fill outo missing sample 11014
#tmp <- data.frame(sample="11014", age=58, sex="M", mutation="Ser82Val")
#md <- rbind(md, tmp)
md <- merge(chip.md, md, by.x="Human_ID", by.y="sample")

# adjust sample name
md$Sample_name <- str_split(md$Sample_name, pattern="sr_", simplify = T)[,2]



# get design file to match with chip-seq file
design <- read.table("../analysis_neuron_290120/rimod_chipseq/design_rimod_chipseq_frontal_neuron.csv", sep=",", header=T)
design$sample <- paste(design$group, paste0("R",design$replicate), sep="_")
design <- design[!grepl("Input", design$sample),]
design$fastq_1 <- str_split(design$fastq_1, pattern="sakibm_", simplify = T)[,2]
design$fastq_1 <- str_split(design$fastq_1, pattern="_S", simplify = T)[,1]

design$fastq_1 <- gsub("/home/kevin/Raw_FASTQ_files_H3K4me3_Frontal_FTLD/", "", design$fastq_1)

design$fastq_1 <- gsub(".fastq.gz", "", design$fastq_1)
design <- design[, c(3, 7)]

# final merge of design and md
md <- merge(md, design, by.x="Sample_name", by.y="fastq_1")
rownames(md) <- md$sample

# Fit metadata for DiffBind
data_dir <- "/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/chipseq/frontal_sakib/results/bwa/mergedLibrary/macs/narrowPeak/"
bam_dir <- "/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/chipseq/frontal_sakib/results/bwa/mergedLibrary/"
md$Peaks <- paste(data_dir, md$sample, "_peaks.xls", sep="")
md$bamReads <- paste(bam_dir, md$sample, ".mLb.clN.sorted.bam", sep="")
md$PeakCaller <- rep("macs", nrow(md))
md$PeakFormat <- rep("macs", nrow(md))
md$SampleID <- md$sample
md$Condition <- md$group
write.csv(md, "chipseq_samplesheet.csv")

####
# MAPT analysis
####
md.mapt <- md[md$group %in% c("GRN", "NDC"),]
write.csv(md.mapt, "MAPT_chipseq_samplesheet.csv")
# load
mapt <- dba(sampleSheet = "MAPT_chipseq_samplesheet.csv")
# count
mapt <- dba.count(mapt)
# contrast
mapt <- dba.contrast(mapt)
# analyze
mapt <- dba.analyze(mapt)
mapt.res <- dba.report(mapt)

# make object

mapt.mask <- dba.mask(frontal, DBA_CONDITION, "MAPT")

frontal <- dba.count(frontal)
