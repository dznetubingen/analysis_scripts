setwd("~/rimod/smallRNA/")

# iPSC neuron count data
ips <- read.table("iPSC/iPSCNeurons_smRNAseq_counts.txt", sep="\t", header=T, row.names=1)
rs <- rowSums(ips,)
rs <- rs[rs > 1000]


# Check highly expressed MAPT DEGs
deg <- read.table("frontal/analysis/analysis_0719/deseq_result_mapt.ndc_frontal_smRNAseq.txt", sep="\t", header=T, row.names = 1)
deg <- deg[deg$padj <= 0.05,]
deg <- deg[rownames(deg) %in% names(rs),]
mapt <- deg

# GRN
deg <- read.table("frontal/analysis/analysis_0719/deseq_result_grn.ndc_frontal_smRNAseq.txt", sep="\t", header=T, row.names = 1)
deg <- deg[deg$padj <= 0.05,]
deg <- deg[rownames(deg) %in% names(rs),]
grn <- deg

# C9orf72
deg <- read.table("frontal/analysis/analysis_0719/deseq_result_c9.ndc_frontal_smRNAseq.txt", sep="\t", header=T, row.names = 1)
deg <- deg[deg$padj <= 0.05,]
deg <- deg[rownames(deg) %in% names(rs),]
c9 <- deg


# look into common DEGs
cmn <- intersect(rownames(mapt), intersect(rownames(grn), rownames(c9)))
cmn <- cmn[cmn %in% names(rs)]


# Get specific target lists
extractTargetGenes <- function(deg, te_table, name="MAPT", out_dir = "frontal/analysis/mir_targets/"){
  mirs <- rownames(deg)
  for (m in mirs){
    tmp <- te_table[te_table$mirna == m,]
    targets <- as.character(tmp$targets)
    write.table(targets, paste0(out_dir, name, m, "_targets.txt", sep="_"), sep="\t", quote=F, row.names=F, col.names=F)
  }
}

# mapt
te_table <- read.table("frontal/analysis/target_mrna_correlation_analysis_0819/MAPT_miRNA_target_edge_table.txt", sep="\t", header=T)
extractTargetGenes(mapt, te_table = te_table, name="MAPT_")


# grn
te_table <- read.table("frontal/analysis/target_mrna_correlation_analysis_0819/GRN_miRNA_target_edge_table.txt", sep="\t", header=T)
extractTargetGenes(mapt, te_table = te_table, name="GRN_")


# c9
te_table <- read.table("frontal/analysis/target_mrna_correlation_analysis_0819/C9_miRNA_target_edge_table.txt", sep="\t", header=T)
extractTargetGenes(mapt, te_table = te_table, name="C9_")
