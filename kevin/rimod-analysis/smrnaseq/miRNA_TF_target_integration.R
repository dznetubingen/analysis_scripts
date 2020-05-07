#############
# Create TF-miRNA target mappings using motif instance finding results from CAGE-seq data
#############
library(stringr)
setwd("~/rimod/CAGE/cage_analysis/tf_enrichment_analysis_050420/")

####
# Function to calculate miRNAs potentially targeted by TFs
#######
calculate_tf_mirna_targets <- function(gtf, mir, tfbs, bp_cutoff=1000){
  
  # format GTF
  colnames(gtf) <- paste("V", c(1:9), sep="")
  ids <- str_split(gtf$V9, pattern=";", simplify = T)
  genes <- c()
  for (i in 1:nrow(ids)){
    for (j in 1:ncol(ids)){
      if (grepl("gene_name", ids[i,j])){
        genes <- c(genes, str_split(ids[i,j], pattern=" ", simplify=T)[,3])
      }
    }
  }
  gtf$gene <- genes
  
  # format miRNAs
  mir <- mir[mir$padj <= 0.05,]
  mirs <- rownames(mir)
  mirs <- gsub("-5p", "", gsub("-3p", "", mirs))
  mirs <- gsub("hsa", "", gsub("-", "", mirs))
  mirs <- gsub("p", "", mirs)
  mirs <- toupper(mirs)
  mir$names <- mirs
  gtf <- gtf[gtf$gene %in% mir$names,]
  gtf <- gtf[gtf$V3 == "gene",]
  
  # format TFBS
  tfbs$chr <- str_split(tfbs$V1, pattern="_", simplify = T)[,1]
  tfbs$start <- as.numeric(str_split(tfbs$V1, pattern="_", simplify = T)[,2])
  
  # Overlap TFBS data with miRNAs
  df <- data.frame(miRNA="dummy", TF="dummy", offset = 0, miRNA_start = 0, peak_start = 0, count = 0, score = 0)
  for (i in 1:nrow(gtf)) {
    m <- gtf[i,]$gene
    m.chr <- gtf[i,]$V1
    m.strand <- gtf[i,]$V7
    m.start <- gtf[i,]$V4
    print(m)
    
    # subset chromosome
    tmp <- tfbs[tfbs$chr == m.chr,]
    tmp <- tmp[tmp$V5 == m.strand,]
    
    for (j in 1:nrow(tmp)){
      t.start <- tmp[j,]$start
      
      # check for proximity
      offset = bp_cutoff*2
      if (m.strand == "+"){
        offset <- m.start - t.start
        peak_before_tss <- m.start > t.start
      }
      else if (m.strand == "-"){
        offset <- t.start - m.start
        peak_before_tss <- m.start < t.start
      }
      if ((offset < bp_cutoff) && peak_before_tss){
        df.tmp <- data.frame(miRNA=m, TF=tmp$V4, offset = abs(m.start - t.start), miRNA_start = m.start, peak_start = t.start, count = 1, score = tmp$V6)
        df <- rbind(df, df.tmp)
      }
    }
  }
  
  df <- df[-1,]
  
  return(df)
}

#==== end function ====#

# mapt up
gtf <- read.table("gencode.v34.miRNA.annotation.gtf", sep="\t")
mir <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_result_mapt.ndc_frontal_smRNAseq.txt")
mir <- mir[mir$log2FoldChange > 0,]
tfbs <- read.table("frontal_tf_activity/mapt/results_up/tf_targets/motif_instances_homer2.txt", sep="\t")
mapt.up <- calculate_tf_mirna_targets(gtf, mir, tfbs, bp_cutoff=2000)
mapt.up <- mapt.up[-1,]
# mapt down
gtf <- read.table("gencode.v34.miRNA.annotation.gtf", sep="\t")
mir <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_result_mapt.ndc_frontal_smRNAseq.txt")
mir <- mir[mir$log2FoldChange < 0,]
tfbs <- read.table("frontal_tf_activity/mapt/results_down/tf_targets/motif_instances_homer2.txt", sep="\t")
mapt.down <- calculate_tf_mirna_targets(gtf, mir, tfbs, bp_cutoff=2000)
mapt.down <- mapt.down[-1,]

# grn up
gtf <- read.table("gencode.v34.miRNA.annotation.gtf", sep="\t")
mir <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_result_grn.ndc_frontal_smRNAseq.txt")
mir <- mir[mir$log2FoldChange > 0,]
tfbs <- read.table("frontal_tf_activity/grn/results_up/tf_targets/motif_instances_homer2.txt", sep="\t")
grn.up <- calculate_tf_mirna_targets(gtf, mir, tfbs, bp_cutoff=2000)
grn.up <- grn.up[-1,]

# grn down
gtf <- read.table("gencode.v34.miRNA.annotation.gtf", sep="\t")
mir <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/deseq_result_grn.ndc_frontal_smRNAseq.txt")
mir <- mir[mir$log2FoldChange < 0,]
tfbs <- read.table("frontal_tf_activity/grn/results_down/tf_targets/motif_instances_homer2.txt", sep="\t")
grn.down <- calculate_tf_mirna_targets(gtf, mir, tfbs, bp_cutoff=2000)
grn.down <- grn.down[-1,]

# Subset the results by score and save
mapt.up <- mapt.up[mapt.up$score > 15,]
mapt.down <- mapt.down[mapt.down$score > 15,]
grn.up <- grn.up[grn.up$score > 15,]
grn.down <- grn.down[grn.down$score > 15,]

write.table(mapt.up, "MAPT_upMiRNAs_TFBS_score15.txt", sep="\t", quote=F)
write.table(mapt.down, "MAPT_downMiRNAs_TFBS_score15.txt", sep="\t", quote=F)
write.table(grn.up, "GRN_upMiRNAs_TFBS_score15.txt", sep="\t", quote=F)
write.table(grn.down, "GRN_downMiRNAs_TFBS_score15.txt", sep="\t", quote=F)

