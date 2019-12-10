#########
# Calculate overlap of predicted regulator transcription factors from CAGE-seq (tf-activity) and RNA-seq (ChEA3)
#########

setwd("~/rimod/integrative_analysis/tf_cage_rnaseq/")


rna.tf <- read.table("~/rimod/RNAseq/analysis/tf_enrichment_chea3/grn_up_Integrated_meanRank.tsv", sep="\t", header=T)
cage.tf <- read.table("~/rimod/CAGE/cage_analysis/tf_enrichment_analysis/frontal_tf_activity/grn/results_up/enrichment/homer2_enrichment_result.txt", header=T, sep="\t",
                      comment.char = "")


# Calculate common TFs as predicted from RNA-seq and CAGE-seq analyses
integrate_tfs <- function(cage, rna, pvalue=0.01, tf_cutoff=50){
  cage.all <- cage
  # subset enriched TFs by P-value
  cage <- cage[cage$p.value < 0.01,]
  
  # consider only first 100 TFs from RNA-seq
  rna <- rna[1:tf_cutoff,]
  
  # Only RNA-seq
  rna.only <- rna[!rna$TF %in% cage.all$Motif.Name,]
  
  # Create final sets of TFs
  common.tf <- intersect(cage$Motif.Name, rna$TF)
  final.tf <- union(common.tf, rna.only$TF)
  
  final.df <- rna[rna$TF %in% final.tf,]
  
  return(final.df)
}

###
# MAPT
###

# up
cage.tf.up <- read.table("~/rimod/CAGE/cage_analysis/tf_enrichment_analysis/frontal_tf_activity/mapt/results_up/enrichment/homer2_enrichment_result.txt", header=T, sep="\t",
                         comment.char = "")
rna.tf.up <- read.table("~/rimod/RNAseq/analysis/tf_enrichment_chea3/mapt_up_Integrated_meanRank.tsv", sep="\t", header=T)
up.tfs <- integrate_tfs(cage.tf.up, rna.tf.up)

# down
cage.tf.down <- read.table("~/rimod/CAGE/cage_analysis/tf_enrichment_analysis/frontal_tf_activity/mapt/results_down/enrichment/homer2_enrichment_result.txt", header=T, sep="\t",
                         comment.char = "")
rna.tf.down <- read.table("~/rimod/RNAseq/analysis/tf_enrichment_chea3/mapt_down_Integrated_meanRank.tsv", sep="\t", header=T)
down.tfs <- integrate_tfs(cage.tf.down, rna.tf.down)

write.table(up.tfs, "MAPT_common_TFs_up.txt", sep="\t", quote=F, row.names = F)
write.table(down.tfs, "MAPT_common_TFs_down.txt", sep="\t", quote=F, row.names = F)

###
# GRN
###

# up
cage.tf.up <- read.table("~/rimod/CAGE/cage_analysis/tf_enrichment_analysis/frontal_tf_activity/grn/results_up/enrichment/homer2_enrichment_result.txt", header=T, sep="\t",
                         comment.char = "")
rna.tf.up <- read.table("~/rimod/RNAseq/analysis/tf_enrichment_chea3/grn_up_Integrated_meanRank.tsv", sep="\t", header=T)
up.tfs <- integrate_tfs(cage.tf.up, rna.tf.up)

# down
cage.tf.down <- read.table("~/rimod/CAGE/cage_analysis/tf_enrichment_analysis/frontal_tf_activity/grn/results_down/enrichment/homer2_enrichment_result.txt", header=T, sep="\t",
                           comment.char = "")
rna.tf.down <- read.table("~/rimod/RNAseq/analysis/tf_enrichment_chea3/grn_down_Integrated_meanRank.tsv", sep="\t", header=T)
down.tfs <- integrate_tfs(cage.tf.down, rna.tf.down)

write.table(up.tfs, "GRN_common_TFs_up.txt", sep="\t", quote=F, row.names = F)
write.table(down.tfs, "GRN_common_TFs_down.txt", sep="\t", quote=F, row.names = F)

###
# C9orf72
###

# up
cage.tf.up <- read.table("~/rimod/CAGE/cage_analysis/tf_enrichment_analysis/frontal_tf_activity/c9orf72/results_up/enrichment/homer2_enrichment_result.txt", header=T, sep="\t",
                         comment.char = "")
rna.tf.up <- read.table("~/rimod/RNAseq/analysis/tf_enrichment_chea3/c9orf72_up_Integrated_meanRank.tsv", sep="\t", header=T)
up.tfs <- integrate_tfs(cage.tf.up, rna.tf.up)

# down
cage.tf.down <- read.table("~/rimod/CAGE/cage_analysis/tf_enrichment_analysis/frontal_tf_activity/c9orf72/results_down/enrichment/homer2_enrichment_result.txt", header=T, sep="\t",
                           comment.char = "")
rna.tf.down <- read.table("~/rimod/RNAseq/analysis/tf_enrichment_chea3/c9orf72_down_Integrated_meanRank.tsv", sep="\t", header=T)
down.tfs <- integrate_tfs(cage.tf.down, rna.tf.down)

write.table(up.tfs, "C9orf72_common_TFs_up.txt", sep="\t", quote=F, row.names = F)
write.table(down.tfs, "C9orf72_common_TFs_down.txt", sep="\t", quote=F, row.names = F)




