#########################
# DE analysis of deconvolve expression
########################
library(DESeq2)

setwd("~/rimod/RNAseq/analysis/deconvolution/")

out_dir = "cell_specific_degs/"

# Load deconvolved 'expression' data
mapt <- read.table("mapt_celltype_expression.txt", sep="\t", header=T, row.names = 1)
grn <- read.table("grn_celltype_expression.txt", sep="\t", header=T, row.names = 1)
c9 <- read.table("c9_celltype_expression.txt", sep="\t", header=T, row.names = 1)
cont <- read.table("control_celltype_expression.txt", sep="\t", header=T, row.names = 1)


celltypes <- colnames(mapt)[1:5]

###
# MAPT
###
ftd.df <- mapt
group_name = "MAPT"
for (ct in celltypes){
  df = data.frame(control = cont[ct], ftd = ftd.df[ct])
  colnames(df) <- c("control", "FTD")
  df <- round(df)
  
  md = data.frame(group = c("control", "FTD"))
  
  dds <- DESeqDataSetFromMatrix(df,
                                colData = md,
                                design = ~ group)
  
  dds$group <- relevel(dds$group, ref = "control")
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  dds <- DESeq(dds)
  res <- na.omit(results(dds, c("group", "FTD", "control")))
  
  deg <- res[res$pvalue <= 0.05,]
  deg.up <- deg[deg$log2FoldChange > 0,]
  deg.down <- deg[deg$log2FoldChange < 0,]
  
  
  write.table(rownames(deg.up), paste0(out_dir, "DEGs_Up_", group_name, "_", ct, ".txt"), quote=F, row.names = F, col.names = F)
  write.table(rownames(deg.down), paste0(out_dir, "DEGs_Down_", group_name, "_", ct, ".txt"), quote=F, row.names = F, col.names = F)
}


###
# GRN
###
ftd.df <- grn
group_name = "GRN"
for (ct in celltypes){
  df = data.frame(control = cont[ct], ftd = ftd.df[ct])
  colnames(df) <- c("control", "FTD")
  df <- round(df)
  
  md = data.frame(group = c("control", "FTD"))
  
  dds <- DESeqDataSetFromMatrix(df,
                                colData = md,
                                design = ~ group)
  
  dds$group <- relevel(dds$group, ref = "control")
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  dds <- DESeq(dds)
  res <- na.omit(results(dds, c("group", "FTD", "control")))
  
  deg <- res[res$pvalue <= 0.05,]
  deg.up <- deg[deg$log2FoldChange > 0,]
  deg.down <- deg[deg$log2FoldChange < 0,]
  
  
  write.table(rownames(deg.up), paste0(out_dir, "DEGs_Up_", group_name, "_", ct, ".txt"), quote=F, row.names = F, col.names = F)
  write.table(rownames(deg.down), paste0(out_dir, "DEGs_Down_", group_name, "_", ct, ".txt"), quote=F, row.names = F, col.names = F)
}


###
# C9orf72
###
ftd.df <- c9
group_name = "C9orf72"
for (ct in celltypes){
  df = data.frame(control = cont[ct], ftd = ftd.df[ct])
  colnames(df) <- c("control", "FTD")
  df <- round(df)
  
  md = data.frame(group = c("control", "FTD"))
  
  dds <- DESeqDataSetFromMatrix(df,
                                colData = md,
                                design = ~ group)
  
  dds$group <- relevel(dds$group, ref = "control")
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  dds <- DESeq(dds)
  res <- na.omit(results(dds, c("group", "FTD", "control")))
  
  deg <- res[res$pvalue <= 0.05,]
  deg.up <- deg[deg$log2FoldChange > 0,]
  deg.down <- deg[deg$log2FoldChange < 0,]
  
  
  write.table(rownames(deg.up), paste0(out_dir, "DEGs_Up_", group_name, "_", ct, ".txt"), quote=F, row.names = F, col.names = F)
  write.table(rownames(deg.down), paste0(out_dir, "DEGs_Down_", group_name, "_", ct, ".txt"), quote=F, row.names = F, col.names = F)
}
