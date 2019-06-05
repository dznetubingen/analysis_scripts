########################################################
# Create Design files for MAPT Mutation comparisons
########################################################
library(stringr)
setwd("~/rimod/RNAseq/as_analysis/leafcutter/group_files/")

md <- read.table("~/rimod/RNAseq/rnaseq_frontal_md.txt", sep="\t", header=T)
rownames(md) <- md$sampleid

# MAPT
mapt <- read.table("mapt_group_file.txt", sep="\t")
samples <- str_split(as.character(mapt$V1), pattern="_", simplify = TRUE)[,1]
md <- md[samples,]

# get mutations
muts <- md[md$group == 'case',]
muts <- levels(factor(muts$mutation))

for (m in muts) {
  case <- md[md$group == 'case',]
  case <- case[case$mutation == m,]
  case$mutated_gene = paste(as.character(case$mutated_gene), m, sep="_")
  control <- md[md$group == 'control',]
  tmp <- rbind(case, control)
  
  # generate group file
  gf <- mapt[samples %in% rownames(tmp),]
  gf$rin <- tmp$RIN
  gf$age <- tmp$AGE
  gf$gender <- tmp$GENDER
  
  filename <- paste("mapt", m, "group_file_confounders.txt", sep="_")
  write.table(gf, filename, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
}

