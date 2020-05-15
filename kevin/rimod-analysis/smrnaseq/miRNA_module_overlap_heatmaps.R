###########
# Create miRNA-module overlap heatmaps
###########

setwd("~/rimod/smallRNA/frontal/analysis/mirna_module_overlap_0520/")


###
# Load necessary data
###
# Load target mappings
mir.mapt <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/MAPT_miRNA_target_edge_table.txt", sep="\t", header=T)
mir.grn <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/GRN_miRNA_target_edge_table.txt", sep="\t", header=T)

mapt.deg <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_result_mapt.ndc_frontal_smRNAseq.txt", sep="\t", header=T)
grn.deg <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_result_grn.ndc_frontal_smRNAseq.txt", sep="\t", header=T)
mapt.deg <- mapt.deg[mapt.deg$padj <= 0.05,]
grn.deg <- grn.deg[grn.deg$padj <= 0.05,]


# Load modules
mod.mapt.down <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_mapt_filtered_down_modules.txt", header=T)
mod.mapt.up <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_mapt_filtered_up_modules.txt", header=T)
mod.grn.down <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_grn_filtered_down_modules.txt", header=T)
mod.grn.up <- read.table("~/rimod/RNAseq/analysis/human_base/rnaseq_grn_filtered_up_modules.txt", header=T)

#===== end of data loading =====#

####
# Function to create IOU heatmaps for miRNAs and modules
####
makeModuleOverlapDF <- function(mod, mirs, tgt_mapping, prefix="MAPT_UP"){
  df <- data.frame(Module = mods)
  rownames(df) <- mods
  
  for (mir in mirs){
    mir.targets <- as.character(tgt_mapping[tgt_mapping$mirna == mir,]$targets)
    if (length(mir.targets) > 0){
      
      mod_ovl <- c()
      for (m in mods) {
        tmp <- as.character(mod[mod$CLUSTER_NAME == m,]$CLUSTER_GENES)
        tmp <- str_split(tmp, pattern=",")[[1]]
        
        u <- length(union(mir.targets, tmp))
        i <- length(intersect(mir.targets, tmp))
        iou <- i / u
        mod_ovl <- c(mod_ovl, i)
      }
      
      df$tmp <- mod_ovl
      colnames(df)[ncol(df)] <- mir
    }
  }
  
  df <- df[,-1]
  df <- t(df)
  
  # remove miRNAs without targets
  rs <- rowSums(df)
  df <- df[rs > 0,]
  df <- df[order(rowSums(df), decreasing = T),]
  colnames(df) <- paste(prefix, colnames(df), sep="_")
  
  return(df)
}
#==== end funciton ====#

####
# Create the heatmaps
mapt.down.mods <- makeModuleOverlapDF(mod.mapt.down, mapt.deg$X, mir.mapt, prefix="MAPT_DOWN")
mapt.up.mods <- makeModuleOverlapDF(mod.mapt.up, mapt.deg$X, mir.mapt, prefix="MAPT_UP")
grn.down.mods <- makeModuleOverlapDF(mod.grn.down, grn.deg$X, mir.grn, prefix="GRN_DOWN")
grn.up.mods <- makeModuleOverlapDF(mod.grn.up, grn.deg$X, mir.grn, prefix="GRN_UP")

# make the actual heatmaps
pheatmap(mapt.up.mods, color=viridis(200), cluster_rows = F, cluster_cols = F, filename = "MAPT_upModules_miRNAovl.png", width=3, height=5)
pheatmap(mapt.down.mods, color=viridis(200), cluster_rows = F, cluster_cols = F, filename = "MAPT_downModules_miRNAovl.png", width=3, height=5)

pheatmap(grn.up.mods, color=viridis(200), cluster_rows = F, cluster_cols = F, filename = "GRN_upModules_miRNAovl.png", width=3, height=5)
pheatmap(grn.down.mods, color=viridis(200), cluster_rows = F, cluster_cols = F, filename = "GRN_downModules_miRNAovl.png", width=3, height=5)