#########
# RiMod TF enrichment results analysis
#########

setwd("~/rimod/CAGE/cage_analysis/tf_enrichment_analysis_050420/")

# Load results
pval <- 0.01
mapt.up <- read.table("frontal_tf_activity/mapt/results_up/enrichment/homer2_enrichment_result.txt", header=T, sep="\t", comment.char = "")
mapt.down <- read.table("frontal_tf_activity/mapt/results_down/enrichment/homer2_enrichment_result.txt", header=T, sep="\t", comment.char = "")
grn.up <- read.table("frontal_tf_activity/grn/results_up/enrichment/homer2_enrichment_result.txt", header=T, sep="\t", comment.char = "")
grn.down <- read.table("frontal_tf_activity/grn/results_down/enrichment/homer2_enrichment_result.txt", header=T, sep="\t", comment.char = "")
# keep only significant TFs
mapt.up <- mapt.up[mapt.up$p.value < pval,]
mapt.down <- mapt.down[mapt.down$p.value < pval,]
grn.up <- grn.up[grn.up$p.value < pval,]
grn.down <- grn.down[grn.down$p.value < pval,]

####
# Make Target lists
####
library(biomaRt)
library(stringr)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Get table of target genes
getTargetGenes <- function(targets){
  targets <- targets[, c(4, 6)]
  targets <- targets[grepl("promoter", targets$Annotation),]
  
  # format the annotation
  tmp <- as.character(targets$Annotation)
  tmp <- str_split(tmp, pattern="[(]", simplify = T)[,2]
  tmp <- gsub(")", "", tmp)
  tmp <- str_split(tmp, pattern="[.]", simplify = T)[,1]
  targets$Annotation <- tmp
  
  # load genes from biomaRt
  bm <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id"), filters="ensembl_transcript_id", values=targets$Annotation, mart=ensembl)
  targets <- merge(targets, bm, by.x="Annotation", by.y="ensembl_transcript_id")
  targets <- targets[!duplicated(targets),]
  colnames(targets) <- c("TargetTranscript", "TF", "TargetGene")
  return(targets)
}

# MAPT up
targets <- read.table("frontal_tf_activity/mapt/results_up/tf_targets/tf_target_mapping_annot.txt", sep="\t", header=T)
targets.mapt.up <- getTargetGenes(targets)
write.table(targets.mapt.up, "MAPT_up_targets.txt", sep="\t", quote=F)
# MAPT down
targets <- read.table("frontal_tf_activity/mapt/results_down/tf_targets/tf_target_mapping_annot.txt", sep="\t", header=T)
targets.mapt.down <- getTargetGenes(targets)
write.table(targets.mapt.down, "MAPT_down_targets.txt", sep="\t", quote=F)
# GRN up
targets <- read.table("frontal_tf_activity/grn/results_up/tf_targets/tf_target_mapping_annot.txt", sep="\t", header=T)
targets.grn.up <- getTargetGenes(targets)
write.table(targets.grn.up, "GRN_up_targets.txt", sep="\t", quote=F)
# GRN down
targets <- read.table("frontal_tf_activity/grn/results_down/tf_targets/tf_target_mapping_annot.txt", sep="\t", header=T)
targets.grn.down <- getTargetGenes(targets)
write.table(targets.grn.down, "GRN_down_targets.txt", sep="\t", quote=F)

#======== end target list creation ===============#



#####
# Make TF-module overlap heatmaps
####
library(pheatmap)
library(viridis)
# Subset TFs for significant ons
targets.mapt.up <- targets.mapt.up[targets.mapt.up$TF %in% mapt.up$Motif.Name,]
targets.mapt.down <- targets.mapt.down[targets.mapt.down$TF %in% mapt.down$Motif.Name,]
targets.grn.up <- targets.grn.up[targets.grn.up$TF %in% grn.up$Motif.Name,]
targets.grn.down <- targets.grn.down[targets.grn.down$TF %in% grn.down$Motif.Name,]

# Load modules
mod.mapt.down <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_mapt_filtered_down_modules.txt", header=T)
mod.mapt.up <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_mapt_filtered_up_modules.txt", header=T)
mod.grn.down <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_grn_filtered_down_modules.txt", header=T)
mod.grn.up <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_grn_filtered_up_modules.txt", header=T)


makeModuleOverlapDF <- function(mod, tgt, prefix="MAPT_UP"){
  
  # get hgnc symbols
  bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=tgt$TargetGene, mart=ensembl)
  tgt <- merge(tgt, bm, by.x = "TargetGene", by.y="ensembl_gene_id")
  
  tfs <- as.character(levels(factor(tgt$TF)))
  mods <- as.character(mod$CLUSTER_NAME)
  
  df <- data.frame(Module = mods)
  rownames(df) <- mods
  
  for (tf in tfs) {
    tf.targets <- tgt[tgt$TF == tf,]$hgnc_symbol
    # calculate module overlap
    mod_ovl <- c()
    for (m in mods) {
      tmp <- as.character(mod[mod$CLUSTER_NAME == m,]$CLUSTER_GENES)
      tmp <- str_split(tmp, pattern=",")[[1]]
      # calc IOU
      u <- length(union(tf.targets, tmp))
      i <- length(intersect(tf.targets, tmp))
      iou <- i / u
      mod_ovl <- c(mod_ovl, iou)
    }
    # save in dataframe
    df$tmp <- mod_ovl
    colnames(df)[ncol(df)] <- tf
  }
  df <- df[, -1]
  df <- t(df)
  rs <- rowSums(df)
  df <- df[order(rs, decreasing = T),]
  colnames(df) <- paste(prefix, colnames(df), sep="_")
  return(df)
}

# MAPT
mapt.up <- makeModuleOverlapDF(mod = mod.mapt.up, tgt = targets.mapt.up, prefix="MAPT_UP")
mapt.down <- makeModuleOverlapDF(mod = mod.mapt.down, tgt = targets.mapt.down, prefix="MAPT_DOWN")

# GRN
grn.up <- makeModuleOverlapDF(mod = mod.grn.up, tgt=targets.grn.up, prefix="GRN_UP")
grn.down <- makeModuleOverlapDF(mod = mod.grn.down, tgt=targets.grn.down, prefix="GRN_DOWN")

# subset for genes for which we have actual expression values
mat <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2020-05-04_15.45.57/deseq_vst_values_2020-05-04_15.45.57.txt", sep="\t", header=T, row.names=1)
genes <- str_split(rownames(mat), pattern="[.]", simplify = T)[,1]
genes <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=genes, mart=ensembl)
genes <- genes$hgnc_symbol

# subsetting
mapt.up <- mapt.up[rownames(mapt.up) %in% genes,]
mapt.down <- mapt.down[rownames(mapt.down) %in% genes,]
grn.up <- grn.up[rownames(grn.up) %in% genes,]
grn.down <- grn.down[rownames(grn.down) %in% genes,]

# subset by differential expression
deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2020-05-04_15.45.57/deseq_result_grn.ndc_fro_2020-05-04_15.45.57.txt", sep="\t", header=T)
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filter = "ensembl_gene_id", values=deg$X, mart=ensembl)
deg <- merge(deg, bm, by.x="X", by.y="ensembl_gene_id")
deg <- deg[deg$padj <= 0.05,]
grn.deg.up <- grn.up[rownames(grn.up) %in% deg$hgnc_symbol,]
grn.deg.down <- grn.down[rownames(grn.down) %in% deg$hgnc_symbol,]

# mapt
deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2020-05-04_15.45.57/deseq_result_mapt.ndc_fro_2020-05-04_15.45.57.txt", sep="\t", header=T)
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filter = "ensembl_gene_id", values=deg$X, mart=ensembl)
deg <- merge(deg, bm, by.x="X", by.y="ensembl_gene_id")
deg <- deg[deg$padj <= 0.05,]

mapt.deg.up <- mapt.up[rownames(mapt.up) %in% deg$hgnc_symbol,]
mapt.deg.down <- mapt.down[rownames(mapt.down) %in% deg$hgnc_symbol,]

# subset a number of TFs
#n_subset = 30
#mapt.up <- mapt.up[1:n_subset,]
#mapt.down <- mapt.down[1:n_subset,]
#grn.up <- grn.up[1:n_subset,]
#grn.down <- grn.down[1:n_subset,]

# make the actual heatmaps
pheatmap(mapt.deg.up, color=viridis(200), cluster_rows = F, cluster_cols = F, filename = "MAPT_up_TFModuleOvl.png", width=3, height=5)
pheatmap(mapt.deg.down, color=viridis(200), cluster_rows = F, cluster_cols = F, filename = "MAPT_down_TFModuleOvl.png", width=3, height=5)

pheatmap(grn.deg.up, color=viridis(200), cluster_rows = F, cluster_cols = F, filename = "GRN_up_TFModuleOvl.png", width=3, height=5)
pheatmap(grn.deg.down, color=viridis(200), cluster_rows = F, cluster_cols = F, filename = "GRN_down_TFModuleOvl.png", width=3, height=5)


