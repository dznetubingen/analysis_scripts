#############
# Flipper Tursiops Variant Analysis
#############
library(vcfR)

setwd("~/flipper/snp_analysis/tursiops/")

# Load VCF file


# Load dna
dna <- ape::read.dna("tursiops_203_ragoo.fasta", format="fasta")

# Subset data to largest scaffold (should be Chromosome 1)
scf_name = "MRVK01001658.1_RaGOO"
scf1 <- dna[grepl(scf_name, names(dna))]


###
# S1
###

# Creat chromR object
vcf1 <- read.vcfR("scaffold1/scf1.genomeS1.vcf.gz")
chrom1 <- create.chromR(name="Scaff1", vcf=vcf1, seq=scf1)

# Plot and mask
chrom1 <- masker(chrom1, min_QUAL = 20, min_MQ=55, max_MQ = 61)

# chromosome QC
chrom1 <- proc.chromR(chrom1, verbose = TRUE, win.size=10000)
chromoqc(chrom1)

#==============================#


###
# S2
###

# Creat chromR object
vcf2 <- read.vcfR("scaffold1/scf1.genomeS2.vcf.gz")
chrom2 <- create.chromR(name="Scaff1", vcf=vcf2, seq=scf1)

# Plot and mask
chrom2 <- masker(chrom2, min_QUAL = 20, min_MQ=55, max_MQ = 61)

# chromosome QC
chrom2 <- proc.chromR(chrom2, verbose = TRUE, win.size=10000)
chromoqc(chrom2)

#==============================#


###
# S1
###

# Creat chromR object
vcf3 <- read.vcfR("scaffold1/scf1.genomeS3.vcf.gz")
chrom <- create.chromR(name="Scaff1", vcf=vcf3, seq=scf1)

# Plot and mask
chrom3 <- masker(chrom3, min_QUAL = 20, min_MQ=55, max_MQ = 61)

# chromosome QC
chrom3 <- proc.chromR(chrom3, verbose = TRUE, win.size=10000)
chromoqc(chrom3)

#==============================#
