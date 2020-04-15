library(biomaRt)
library(stringr)

# Load GTF and make gene table
gtf <- read.table("~/dzne/rimod_package/gencode.v33.annotation.gtf", sep="\t")
gtf <- gtf[gtf$V3 == "gene",]
genes <- gtf$V9
genes <- str_split(genes, pattern=";", simplify=T)

as <- read.table("~/dzne/rimod_package/as_analysis/psi/mapt_control.deltapsi.tsv", sep="\t", header=T, stringsAsFactors = F)
#as <- read.table("~/rimod/RNAseq/as_analysis/majiq/majiq_psi/grn_control.deltapsi.tsv", sep="\t", header=T, stringsAsFactors = F)


cpg <- read.table("~/dzne/rimod_package/frontal_methylation_0818/DMPs_mapt.ndc_quant.txt", sep="\t", header = T, stringsAsFactors = F)
#cpg <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_grn.ndc_quant.txt", sep="\t", header = T, stringsAsFactors = F)

cpg <- cpg[cpg$adj.P.Val <= 0.05,]

distance <- 500

# Iterate over all junctions
df <- data.frame(Junction = "dummy", CpG = "dummy", Gene = "dummy")

for (i in 1:nrow(as)){
  print(paste("Juncion:", i))
  # Get the significant junctions at index i
  tmp <- as[i,]
  pvals <- as.numeric(str_split(tmp$P..dPSI...0.20..per.LSV.junction, pattern=";")[[1]])
  juncs <- str_split(tmp$Junctions.coords, pattern=";")[[1]]
  gene = str_split(tmp$Gene.ID, pattern="[.]", simplify = T)[,1]
  
  if (any(pvals) <= 0.05){
    
  }
  keep <- pvals <= 0.05
  juncs <- juncs[keep]
  
  for (j in 1:length(juncs)) {
    
  }

  if (!chr %in% c("chrX", "chrY")){
    # for each junction, check if it overlaps a CpG site
    if (length(juncs) > 0){
      for (j in 1:length(juncs)) {
        j <- as.numeric(str_split(juncs[j], pattern="-")[[1]])
        start <- j[1]
        end <- j[2]
        
        # Check all CpG sites on that chromosome
        cpg.chr <- cpg[cpg$chr == chr,]
        for (c in 1:nrow(cpg.chr)) {
          pos <- cpg.chr[c,]$pos
          diff.start <- pos - start
          diff.end <- pos - end
          if (abs(diff.start) < distance) {
            tmp.df <- data.frame(Junction = i, CpG = cpg.chr[c,]$Name, Gene = tmp$Gene.Name)
            df <- rbind(df, tmp.df)
          } else if (abs(diff.end) < distance){
            tmp.df <- data.frame(Junction = i, CpG = cpg.chr[c,]$Name, Gene = tmp$Gene.Name)
            df <- rbind(df, tmp.df)
          }
        }
      }
    }
  }

}


# VAC14
vac <- as[as$Gene.Name == "VAC14",]
vac.cpg <- cpg[cpg$Name == "cg12377220",]

# GIT2
gits <- df[df$Gene == "GIT2",]
git <- as[as$Gene.Name == "GIT2",]
git.cpg <- cpg[cpg$Name %in% gits$CpG,]
