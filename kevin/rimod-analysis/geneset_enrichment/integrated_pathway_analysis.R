##############################
# Combined Pathway analysis
# Compare results from different omics to look for conistently deregulated pathways
##############################
setwd("~/rimod/integrative_analysis/pathway_analysis_integration/")

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
cmn.up.df <- rna.grn.up[rna.grn.up$term_name %in% cmn.up,][, c(1:3)]

# Common down-regulated pythways
cmn.down <- intersect(grn.down, mapt.down)
cmn.down.df <- rna.grn.down[rna.grn.down$term_name %in% cmn.down,][, c(1:3)]


###
# miRNA target comparison
# up-miRNAs --> down-regulated genes
tdown <- intersect(mir.c.tdown$term_name, intersect(mir.grn.tdown$term_name, mir.mapt.tdown$term_name))
tdown.df <- mir.c.tdown[mir.c.tdown$term_name %in% tdown,][,c(1:3)]
# get de-regulated miRNAs
grn <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/DEGs_P0.05_LFC0.8_grn.ndc_frontal_smRNAseq_miRNAs.txt", header=T)
mapt <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/DEGs_P0.05_LFC0.8_mapt.ndc_frontal_smRNAseq_miRNAs.txt", header=T)
c9 <- read.table("~/rimod/smallRNA/frontal/analysis/analysis_0719/DEGs_P0.05_LFC0.8_c9.ndc_frontal_smRNAseq_miRNAs.txt", header=T)
cmn.mirs <- intersect(grn$x, intersect(mapt$x, c9$x))
mir.enr <- tdown.df[tdown.df$term_name %in% cmn.mirs,]

# down-miRNAs
tup <- intersect(mir.c.tup$term_name, intersect(mir.grn.tup$term_name, mir.mapt.tup$term_name))
tup.df <- mir.c.tup[mir.c.tup$term_name %in% tup,][,c(1:3)]

