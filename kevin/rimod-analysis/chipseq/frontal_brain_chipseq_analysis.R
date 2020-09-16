#################
# RiMod Frontal Neuron ChIP-seq analysis
#################
library(DESeq2)
library(stringr)
setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/rimod/chipseq/brain_cemil/")


# Load data
cts <- read.table("results/bwa/mergedLibrary/macs/narrowPeak/consensus/H3K4me3/deseq2/H3K4me3.consensus_peaks.featureCounts.txt", header=T, sep="\t")
peak_id = paste(cts$Geneid, cts$Chr, cts$Start, cts$End, cts$Strand, sep="_")
rownames(cts) <- peak_id
cts <- cts[, c(-1, -2, -3, -4, -5, -6)]
colnames(cts) <- str_split(colnames(cts), pattern="[.]", simplify = T)[,1]
##
# Create Metadata file
##
chip.md <- read.csv("../design_chipseq_frontal_neuron_cemil.csv")
chip.md$sample <- str_split(chip.md$fastq_1, pattern="/", simplify = T)[,5]
chip.md <- chip.md[!grepl("Input", chip.md$group),]
group <- str_split(chip.md$sample, pattern="_", simplify = T)[,1]
chip.md$sample_id <- str_sub(group, start=1, end=5)
group[grepl("C9", group)] <- "C9"
group[grepl("GRN", group)] <- "GRN"
group[grepl("MAPT", group)] <- "MAPT"
group[grepl("NDC", group)] <- "NDC"
id <- str_split(chip.md$sample, pattern="_", simplify = T)[,7]
sample <- paste0("Brain_", group, "_R", id)
chip.md$sample <- sample


# rimod overall sample sheet
md <- read.csv("~/rimod/files/FTD_Brain_corrected.csv")
md$SAMPLEID <- str_pad(md$SAMPLEID, width = 5, side = "left", pad = "0")
md <- md[md$SAMPLEID %in% chip.md$sample_id,]
md <- md[md$REGION == "frontal",]
md <- data.frame(sample=md$SAMPLEID, age=md$AGE, sex=md$GENDER, mutation=md$GENE)

md <- merge(chip.md, md, by.x="sample_id", by.y="sample")

# match files
cts <- cts[,colnames(cts) %in% md$sample]
cts <- cts[, match(md$sample, colnames(cts))]

#==========================================================================#

# Make DESeq2 object
dds <- DESeqDataSetFromMatrix(cts,
                              colData = md,
                              design = ~ sex + group)
# specifcy control group
dds$group <- relevel(dds$group, ref = "Brain_NDC")

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
res.mapt <- results(dds, c("group", "Brain_MAPT", "Brain_NDC"))
res.mapt <- na.omit(res.mapt)
deg.mapt <- res.mapt[res.mapt$padj <= pval_cut,]
print(dim(deg.mapt))

### GRN - control
res.grn <- results(dds, c("group", "Brain_GRN", "Brain_NDC"))
res.grn <- na.omit(res.grn)
deg.grn <- res.grn[res.grn$padj <= pval_cut,]
print(dim(deg.grn))

### C9orf72 - control
res.c9 <- results(dds, c("group", "Brain_C9", "Brain_NDC"))
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
<<<<<<< HEAD
write.table(as.character(test$Nearest.PromoterID), "~/tmp_stuff/grn_chipseq.txt", quote=F, row.names = F, col.names = F)
=======
write.table(as.character(test$Nearest.PromoterID), "~/tmp_stuff/grn_chipseq.txt", quote=F, row.names = F, col.names = F)
>>>>>>> a3109a993c6d5817cd92ef53a98e72c0ea2d4fba
