#################
# RiMod Frontal Neuron ChIP-seq analysis
#################
library(DESeq2)
library(stringr)
library(IHW)
setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/chipseq/frontal_sakib/")


# Load data
cts <- read.table("results/bwa/mergedLibrary/macs/narrowPeak/consensus/H3K4me3/deseq2/H3K4me3.consensus_peaks.featureCounts.txt", header=T, sep="\t")
peak_id = paste(cts$Geneid, cts$Chr, cts$Start, cts$End, cts$Strand, sep="_")
rownames(cts) <- peak_id
cts <- cts[, c(-1, -2, -3, -4, -5, -6)]
colnames(cts) <- str_split(colnames(cts), pattern="[.]", simplify = T)[,1]
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
design <- design[design$sample %in% colnames(cts),]
design$fastq_1 <- str_split(design$fastq_1, pattern="sakibm_", simplify = T)[,2]
design$fastq_1 <- str_split(design$fastq_1, pattern="_S", simplify = T)[,1]

design$fastq_1 <- gsub("/home/kevin/Raw_FASTQ_files_H3K4me3_Frontal_FTLD/", "", design$fastq_1)

design$fastq_1 <- gsub(".fastq.gz", "", design$fastq_1)
design <- design[, c(3, 7)]

# final merge of design and md
md <- merge(md, design, by.x="Sample_name", by.y="fastq_1")
rownames(md) <- md$sample
md <- md[match(colnames(cts), md$sample),]

# testing
md$group[14] <- "NDC"

#==========================================================================#

# Make DESeq2 object
dds <- DESeqDataSetFromMatrix(cts,
                              colData = md,
                              design = ~ sex + group)
# specifcy control group
dds$group <- relevel(dds$group, ref = "NDC")

# pre-filtering
count_cutoff <- 10
sample_cutoff <- 4
dds <- estimateSizeFactors(dds)
keep <- rowSums((counts(dds, normalized=TRUE) >= count_cutoff)) >= sample_cutoff
table(keep)

dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Exctract results
pval_cut <- 0.05
res.mapt <- results(dds, c("group", "MAPT", "NDC"), filterFun = )
res.mapt <- na.omit(res.mapt)
deg.mapt <- res.mapt[res.mapt$padj <= pval_cut,]
print(dim(deg.mapt))

### GRN - control
res.grn <- results(dds, c("group", "GRN", "NDC"), filterFun = ihw)
res.grn <- na.omit(res.grn)
deg.grn <- res.grn[res.grn$padj <= pval_cut,]
print(dim(deg.grn))

### C9orf72 - control
res.c9 <- results(dds, c("group", "C9orf72", "NDC"), filterFun = ihw)
res.c9 <- na.omit(res.c9)
deg.c9 <- res.c9[res.c9$padj <= pval_cut,]
print(dim(deg.c9))


####
# Plotting
####
vst <- varianceStabilizingTransformation(dds)
plotPCA(vst, intgroup="group")

library(limma)
des <- model.matrix( ~ md$group)
res <- removeBatchEffect(assay(vst), batch=md$sex, design=des)
assay(vst) <- res
plotPCA(vst, intgroup="group", ntop=10000)


# Perform annotation
anno <- read.table("results/bwa/mergedLibrary/macs/narrowPeak/consensus/H3K4me3/H3K4me3.consensus_peaks.annotatePeaks.txt", sep="\t", header=T)
colnames(anno)[1] <- "interval"
#test <- deg.grn[deg.c9$log2FoldChange < 0,]
test <- rownames(deg.grn)
test <- str_split(test, pattern="_chr", simplify = T)[,1]
test <- anno[anno$interval %in% test,]
write.table(as.character(test$Nearest.PromoterID), "~/tmp_stuff/grn_chipseq.txt", quote=F, row.names = F, col.names = F)
