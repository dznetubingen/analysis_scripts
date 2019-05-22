#################################################
## Integration of miRNA data with results from
## regulon analysis
#################################################
library(biomaRt)

# Testing

setwd("~/rimod/integrative_analysis/regulon_miRNA/mapt/")

# Load regulons
regs <- readRDS("~/rimod/CAGE/cage_analysis/regulon_analysis_250418/mapt/mapt_sig_regulons.rds")
regs.down <- readRDS("~/rimod/CAGE/cage_analysis/regulon_down_analysis_300418/mapt/mapt_sig_regulons.rds")
# Load DEGs
deg <- read.table("deseq_result_mapt.ndc_2018-02-28_09.22.50.txt", sep="\t", header=T, row.names=1)
deg <- deg[deg$padj <= 0.05,]
deg.up <- deg[deg$log2FoldChange >= 1,]
deg.down <- deg[deg$log2FoldChange <= -1,]
# load expression values
mirna <- read.table("../deseq_rLog_values2018-02-28_09.22.50.txt", sep="\t", row.names=1, header=T)
cage <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2018-04-26_genewise/deseq_rLog_values_2018-04-26_14.22.04.txt", sep="\t", row.names=1, header=T)
# Read targets
print("Loading targets ...")
targets <- read.csv("mapt_deg_miRWalk_miRNA_Targets.csv")
targets <- targets[targets$TargetScan == 1,]
# Load mart
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get ENSEMBL IDs for TFs
up.tfs <- names(regs)
down.tfs <- names(regs.down)
up.tfs <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "hgnc_symbol", values = up.tfs, mart = ensembl)
down.tfs <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "hgnc_symbol", values = down.tfs, mart = ensembl)


### Check overlap of down miRNAs-targets with up-regulated regulons
down.mirs <- rownames(deg.down)
dm.targets <- targets[targets$mirnaid %in% down.mirs,]
dmts.transcript <- as.character(levels(factor(dm.targets$ensemblid)))
bm <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id"), filters = "ensembl_transcript_id", values = dmts.transcript, mart = ensembl)
dmts.genes <- bm$ensembl_gene_id
dm.overlap <- dmts.genes[dmts.genes %in% unlist(regs)]
dm.tf.overlap <- up.tfs[up.tfs$ensembl_gene_id %in% dmts.genes,]

### Check overlap of up miRNAs-targest with down-regulated regulons
up.mirs <- rownames(deg.up)
um.targets <- targets[targets$mirnaid %in% up.mirs,]
upts.transcripts <- as.character(levels(factor(um.targets$ensemblid)))
bm <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id"), filters = "ensembl_transcript_id", values = upts.transcripts, mart = ensembl)
umts.genes <- bm$ensembl_gene_id
um.overlap <- umts.genes[umts.genes %in% unlist(regs)]
up.tf.overlap <- down.tfs[down.tfs$ensembl_gene_id %in% umts.genes, ]

