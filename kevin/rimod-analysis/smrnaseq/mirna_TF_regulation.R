################
# Look for TFBS of active TFs around miRNA promoters
###############
library(stringr)
base_dir <- "/Users/kevin/dzne/rimod_analysis/mirna_tf_regulation/"
setwd(base_dir)

# load GTF file
gtf <- read.table("/Users/kevin/dzne/gencode.mirna.gtf", sep="\t", stringsAsFactors = F)
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
#gtf <- gtf[gtf$V3 == "gene",]

# load miRNA DE result
mir <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/analysis_0719/deseq_result_grn.ndc_frontal_smRNAseq.txt", sep="\t", header=T, row.names=1, stringsAsFactors = F)
mir <- mir[mir$padj <= 0.05,]

# manually modify miRNA names
mirs <- rownames(mir)
mirs <- gsub("-5p", "", gsub("-3p", "", mirs))
mirs <- gsub("hsa", "", gsub("-", "", mirs))
mirs <- gsub("p", "", mirs)
mirs <- toupper(mirs)
mir$names <- mirs

# Overlap miRNAs with GTF to get miRNA positions
gtf <- gtf[gtf$gene %in% mir$names,]


# Load TF information
motifs <- read.table("/Users/kevin/dzne/tf_enrichment_analysis_050420/frontal_tf_activity/grn/results_up/tf_targets/motif_instances_homer2.txt", sep="\t")


tf <- tf[grepl("promoter", tf$Annotation),]






# first look for any kind of overlap
tag <- as.character(tf$Annotation)
tag <- str_split(tag, pattern="[(]", simplify = T)[,2]
tag <- gsub(")", "", tag)
tf$tag <- tag

gtf <- read.table("/Users/kevin/dzne/gencode.mirna.gtf", sep="\t", stringsAsFactors = F)
gtf <- gtf[gtf$V3 == "transcript",]
ids <- str_split(gtf$V9, pattern=";", simplify = T)
genes <- c()
for (i in 1:nrow(ids)){
  for (j in 1:ncol(ids)){
    if (grepl("transcript_id", ids[i,j])){
      genes <- c(genes, str_split(ids[i,j], pattern=" ", simplify=T)[,3])
    }
  }
}
gtf$transcript <- genes

ovl <- intersect(tag, genes)

test <- gtf[gtf$transcript %in% ovl,]
tf <- tf[tf$tag %in% ovl,]

tf2 <- tf[tf$tag == "ENST00000582661.1",]
