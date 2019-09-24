##############################
# Combined Pathway analysis
# Compare results from different omics to look for conistently deregulated pathways
##############################
library(stringr)
library(biomaRt)
setwd("~/rimod/integrative_analysis/pathway_analysis_integration/")


saveGenelist <- function(df, term_name, path="~/rimod/RNAseq/analysis/deconvolution/investigative_genelists/"){
  
  tmp <- df[df$term_name == term_name,]
  genes <- c()
  if (nrow(tmp) > 1){
    for (i in 1:nrow(tmp)){
      genes <- c(genes, str_split(tmp$intersections[i], pattern=",", simplify = T)[1,])
    }
  }
  else {
    genes <- c(genes, str_split(tmp$intersections, pattern=",", simplify = T)[1,])
  }
  genes <- genes[!duplicated(genes)]
  
  write.table(genes, paste0(path, term_name, "_intersection_genes.txt"), quote=F, row.names = F)
  
  return(genes)
}




# Load RNA-seq
rna.mapt.up <- read.csv("gProfiler_RNAseq_MAPT_up.csv")
rna.grn.up <- read.csv("gProfiler_RNAseq_GRN_up.csv")
rna.c.up <- read.csv("gProfiler_RNAseq_C9orf72_NS_up.csv")
rna.mapt.down <- read.csv("gProfiler_RNAseq_MAPT_down.csv")
rna.grn.down <- read.csv("gProfiler_RNAseq_GRN_down.csv")
rna.c.down <- read.csv("gProfiler_RNAseq_C9orf72_NS_down.csv")

# Load CAGEseq
cage.mapt.up <- read.csv("gProfiler_CAGEseq_MAPT_up.csv")
cage.mapt.down <- read.csv("gProfiler_CAGEseq_MAPT_down.csv")
cage.grn.up <- read.csv("gProfiler_CAGEseq_GRN_up.csv")
cage.grn.down <- read.csv("gProfiler_CAGEseq_GRN_down.csv")
cage.c.up <- read.csv("gProfiler_CAGEseq_C9orf72_up.csv")

# Load miRNA
mir.mapt.tup <- read.csv("gProfiler_miRNA_MAPT_upTargets.csv")
mir.mapt.tdown <- read.csv("gProfiler_miRNA_MAPT_downTargets.csv")
mir.grn.tup <- read.csv("gProfiler_miRNA_GRN_upTargets.csv")
mir.grn.tdown <- read.csv("gProfiler_miRNA_GRN_downTargets.csv")
mir.c.tup <- read.csv("gProfiler_miRNA_C9orf72_upTargets.csv")
mir.c.tdown <- read.csv("gProfiler_miRNA_C9orf72_downTargets.csv")

# Load alternative splicing
as.mapt <- read.csv("gProfiler_AS_MAPT_dPSI0.2.csv")
as.grn <- read.csv("gProfiler_AS_GRN_dPSI0.2.csv")
as.c <- read.csv("gProfiler_AS_C9orf72_dPSI0.2.csv")

# Promotor shifting
ps.mapt <- read.csv("gProfiler_PromShift_MAPT.csv")
ps.grn <- read.csv("gProfiler_PromShift_GRN.csv")
ps.c <- read.csv("gProfiler_PromShift_C9orf72.csv")


####
# MAPT comparison
####
mapt.up <- intersect(rna.mapt.up$term_name, cage.mapt.up$term_name)
mapt.down <- intersect(rna.mapt.down$term_name, cage.mapt.down$term_name)
mapt.down.mir <- intersect(mapt.down, mir.mapt.tdown$term_name)

####
# GRN comparison
####
grn.up <- intersect(rna.grn.up$term_name, cage.grn.up$term_name)
grn.down <- intersect(rna.grn.down$term_name, cage.grn.down$term_name)

###
# C9orf72 comparison
###
c.up <- intersect(rna.c.up$term_name, cage.c.up$term_name)
c.down <- rna.c.down$term_name



# Common up-regulated pythways
cmn.up <- intersect(grn.up, mapt.up)
cmn.up.df <- rna.grn.up[rna.grn.up$term_name %in% cmn.up,]

# Save Up-genelists
saveGenelist(cmn.up.df, "extracellular matrix organization")

# Common down-regulated pythways
cmn.down <- intersect(grn.down, mapt.down)
cmn.down.df <- rna.grn.down[rna.grn.down$term_name %in% cmn.down,]

saveGenelist(cmn.down.df, "Oxidative phosphorylation")
saveGenelist(cmn.down.df, "vesicle")
saveGenelist(cmn.down.df, "envelope")



###
# miRNA target comparison
# up-miRNAs --> down-regulated genes
tdown <- intersect(mir.c.tdown$term_name, intersect(mir.grn.tdown$term_name, mir.mapt.tdown$term_name))
tdown.df <- mir.c.tdown[mir.c.tdown$term_name %in% tdown,]

# save genelists
saveGenelist(tdown.df, "Membrane Trafficking")
saveGenelist(tdown.df, "hsa-miR-19b-3p")
saveGenelist(tdown.df, "protein binding")

# get de-regulated miRNAs
grn <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/DEGs_P0.05_LFC0.8_grn.ndc_frontal_smRNAseq_miRNAs.txt", header=T)
mapt <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/DEGs_P0.05_LFC0.8_mapt.ndc_frontal_smRNAseq_miRNAs.txt", header=T)
c9 <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/DEGs_P0.05_LFC0.8_c9.ndc_frontal_smRNAseq_miRNAs.txt", header=T)
cmn.mirs <- intersect(grn$x, intersect(mapt$x, c9$x))
mir.enr <- tdown.df[tdown.df$term_name %in% cmn.mirs,]

# down-miRNAs
tup <- intersect(mir.c.tup$term_name, intersect(mir.grn.tup$term_name, mir.mapt.tup$term_name))
tup.df <- mir.c.tup[mir.c.tup$term_name %in% tup,]

####
# MAPT specific analysis
####
mapt.up <- intersect(rna.mapt.up$term_name, cage.mapt.up$term_name)
mapt.down <- intersect(rna.mapt.down$term_name, cage.mapt.down$term_name)

####
# GRN specific analysis
####
grn.up <- intersect(rna.grn.up$term_name, cage.grn.up$term_name)
grn.down <- intersect(rna.grn.down$term_name, cage.grn.down$term_name)

####
# C9orf72 specific analysis
####
c9.up <- intersect(rna.c.up$term_name, cage.c.up$term_name)
c9.down <- rna.c.down$term_name


##########
# Active TF investigation
##########
# Have a closer look at active transcription factors (in frontal lobe)
tf_dir = "~/rimod/CAGE/cage_analysis/chea_tf_regulators/frontal/"

# laod data
tf.mapt.up <- read.table("~/rimod/CAGE/cage_analysis/chea_tf_regulators/frontal/Integrated_meanRank_MAPT_up.tsv", sep="\t", header=T)
tf.mapt.down <- read.table("~/rimod/CAGE/cage_analysis/chea_tf_regulators/frontal/Integrated_meanRank_MAPT_down.tsv", sep="\t", header=T)
tf.grn.up <- read.table("~/rimod/CAGE/cage_analysis/chea_tf_regulators/frontal/Integrated_meanRank_GRN_up.tsv", sep="\t", header=T)
tf.grn.down <- read.table("~/rimod/CAGE/cage_analysis/chea_tf_regulators/frontal/Integrated_meanRank_GRN_down.tsv", sep="\t", header=T)
tf.c.up <- read.table("~/rimod/CAGE/cage_analysis/chea_tf_regulators/frontal/Integrated_meanRank_c9orf72_up.tsv", sep="\t", header=T)
tf.c.down <- read.table("~/rimod/CAGE/cage_analysis/chea_tf_regulators/frontal/Integrated_meanRank_c9orf72_down.tsv", sep="\t", header=T)

# Select the first 50 TFs
cutoff = 50
tf.mapt.up <- tf.mapt.up[1:cutoff,]
tf.mapt.down <- tf.mapt.down[1:cutoff,]
tf.grn.up <- tf.grn.up[1:cutoff,]
tf.grn.down <- tf.grn.down[1:cutoff,]
tf.c.up <- tf.c.up[1:cutoff,]
tf.c.down <- tf.c.down[1:cutoff,]


# subset with differntially regulated TFs
mapt.deg <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2019-08-15_14.49.08_frontal/deseq_result_mapt.ndc_2019-08-15_14.49.08.txt", sep="\t", header=T)
grn.deg <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2019-08-15_14.49.08_frontal/deseq_result_grn.ndc_2019-08-15_14.49.08.txt", sep="\t", header=T)
c9.deg <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2019-08-15_14.49.08_frontal/deseq_result_c9.ndc_2019-08-15_14.49.08.txt", sep="\t", header=T)
mapt.deg <- mapt.deg[mapt.deg$padj <= 0.05,]
grn.deg <- grn.deg[grn.deg$padj <= 0.05,]
c9.deg <- c9.deg[c9.deg$padj <= 0.05,]

# ensembl to HGNC
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
genesToHgnc <- function(genes, mart){
  bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=genes, mart = mart)
  return(bm$hgnc_symbol[!duplicated(bm$hgnc_symbol)])
}

mapt.degs <- genesToHgnc(mapt.deg$X, ensembl) 
grn.degs <- genesToHgnc(grn.deg$X, ensembl)
c9.degs <- genesToHgnc(c9.deg$X, ensembl)

# Get intersection of TFs with DEGs
tf.mapt.up <- tf.mapt.up[tf.mapt.up$TF %in% mapt.degs,]
tf.mapt.down <- tf.mapt.down[tf.mapt.down$TF %in% mapt.degs,]
tf.grn.up <- tf.grn.up[tf.grn.up$TF %in% grn.degs,]
tf.grn.down <- tf.grn.down[tf.grn.down$TF %in% grn.degs,]
tf.c.up <- tf.c.up[tf.c.up$TF %in% c9.degs,]
tf.c.down <- tf.c.down[tf.c.down$TF %in% c9.degs,]





cmn.tf.up = intersect(tf.mapt.up$TF, intersect(tf.grn.up$TF, tf.c.up$TF))
cmn.tf.down <- intersect(tf.mapt.down$TF, intersect(tf.grn.down$TF, tf.c.down$TF))
