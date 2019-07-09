##################################
# Analysis of MAJIQ results
##################################
library(stringr)
setwd("~/rimod/RNAseq/as_analysis/majiq/")

##
# Function to generate a keep vector
##
generateKeepVector <- function(dpsi, dpsi_cutoff = 0.05){
  keep_vec = rep(FALSE, length(dpsi))
  count = 0
  for (elem in dpsi) {
    count = count + 1
    tmp = c()
    for (e in elem){
      tmp = c(tmp, as.numeric(e))
    }
    if (any(abs(tmp) > dpsi_cutoff)){
      keep_vec[count] <- TRUE
    }
  }  
  return(keep_vec)
}

# helper function to get DPSI values in nice format
getDPSI <- function(mat){
  dpsi <- as.character(mat$`E(dPSI) per LSV junction`)
  dpsi <- sapply(dpsi, function(x){strsplit(x, split=";")})
  return(dpsi)
}

# helper function to split genes
splitGenenames <- function(genes){
  res <- as.character(str_split(genes, pattern="[.]", simplify = TRUE)[,1])
  return(res)
}


# remove dot from gene

# Define dPSI cutoff
dpsi_cutoff = 0.1

# MAPT results
mapt <- read.table("majiq_psi/mapt_control.deltapsi.tsv", sep="\t", header=T, check.names=F)
mapt <- mapt[generateKeepVector(getDPSI(mapt), dpsi_cutoff = dpsi_cutoff),]

# GRN results
grn <- read.table("majiq_psi/grn_control.deltapsi.tsv", sep="\t", header=T, check.names=F)
grn <- grn[generateKeepVector(getDPSI(grn), dpsi_cutoff = dpsi_cutoff),]

# MAPT results
c9orf72 <- read.table("majiq_psi/c9orf72_control.deltapsi.tsv", sep="\t", header=T, check.names=F)
c9orf72 <- c9orf72[generateKeepVector(getDPSI(c9orf72), dpsi_cutoff = dpsi_cutoff),]

# Save differentially splice genes
mapt_genes <- splitGenenames(mapt$`Gene ID`)
grn_genes <- splitGenenames(grn$`Gene ID`)
c9orf72_genes <- splitGenenames(c9orf72$`Gene ID`)

write.table(mapt_genes, paste0("mapt_AS_genes_dPSI_", dpsi_cutoff, ".txt"), sep="\t", quote=F, row.names = FALSE)
write.table(grn_genes, paste0("grn_AS_genes_dPSI_", dpsi_cutoff, ".txt"), sep="\t", quote=F, row.names = FALSE)
write.table(c9orf72_genes, paste0("c9orf72_AS_genes_dPSI_", dpsi_cutoff, ".txt"), sep="\t", quote=F, row.names = FALSE)




