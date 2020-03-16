#############################
# Integration of methylation and promoter shifting results
#############################



setwd("~/rimod/integrative_analysis/methylation_promoter_analysis/")


# load promshift
ps <- read.table("~/rimod/CAGE/cage_analysis/promotor_shifting/frontal/grn_shifting_promotors.txt", sep="\t", header=T, stringsAsFactors = F)
ps <- ps[ps$fdr.KS <= 0.05,]

# load CpGs
met <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_grn.ndc_quant.txt", sep="\t", header=T, stringsAsFactors = F)
met <- met[met$adj.P.Val <= 0.05,]

threshold <- 500
# For each promoeter that is shifting, check for DMPs in the vicinity
df <- data.frame("PromIdx" = 1, "CpG" = "asdf")

i <- 1

for (i in 1:nrow(ps)) {
  print(i)
  # extract promoter information
  prom <- ps[i,]
  chr <- prom$chr
  start <- prom$start
  end <- prom$end
  strand <- prom$strand
  
  # subset relevant CpGs
  tmp <- met[met$chr == chr,]
  tmp <- tmp[tmp$strand == strand,]
  
  # Iterate over all CpGs
  for (j in 1:nrow(tmp)) {
    cpg <- tmp[j,]
    if (abs(start - cpg$pos) <= threshold){
      print("case")
      df <- rbind(df, data.frame("PromIdx" = i, "CpG" = cpg$Name))
    }
    else if (abs(end - cpg$pos) <= threshold){
      df <- rbind(df, data.frame("PromIdx" = i, "CpG" = cpg$Name))
      print("case")
    }
  }
  
}
df <- df[-1,]

idxs <- as.numeric(df$PromIdx)
ps[idxs,]

idx = df$PromIdx[1]
prom <- ps[idx,]
cpgs <- as.character(df[df$PromIdx == idx,]$CpG)
cpgs <- met[met$Name %in% cpgs,]
cpgs
