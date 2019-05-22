###################################
## Regulon Comparison
##
## Compare regulons from different 
##
###################################

# load libs
library(VennDiagram)
library(topGO)

# Change this later ...
setwd("~/rimod/CAGE/cage_analysis/tf_enrichment_pipeline/")

# Load regulons
c9regs <- readRDS("c9orf72/c9orf72_sig_regulons.rds")
maptregs <- readRDS("mapt/mapt_sig_regulons.rds")
grnregs <- readRDS("grn/grn_sig_regulons.rds")

# Get sizes
c9regs.sizes <- sapply(c9regs, length)
maptregs.sizes <- sapply(maptregs, length)
grnregs.sizes <- sapply(grnregs, length)

cmn <- intersect(names(c9regs), intersect(names(maptregs), names(grnregs)))

# GRN and MAPT
gm.cmn <- intersect(names(maptregs), names(grnregs))
mrs.sub <- maptregs.sizes[names(maptregs.sizes) %in% gm.cmn]
grs.sub <- grnregs.sizes[names(grnregs.sizes) %in% gm.cmn]

# exclusive regulons
c9.only <- c9regs[!names(c9regs) %in% c(names(maptregs), names(grnregs))]
mapt.only <- maptregs[!names(maptregs) %in% c(names(c9regs), names(grnregs))]
grn.only <- grnregs[!names(grnregs) %in% c(names(maptregs), names(c9regs))]

# gene based analysis
mgenes <- as.character(unlist(maptregs))
cgenes <- as.character(unlist(c9regs))
ggenes <- as.character(unlist(grnregs))

mg.only <- mgenes[!mgenes %in% c(cgenes, ggenes)]
cg.only <- cgenes[!cgenes %in% c(mgenes, ggenes)]
gg.only <- ggenes[!ggenes %in% c(mgenes, cgenes)]

cmn.genes <- intersect(mgenes, intersect(cgenes, ggenes))

write(mg.only, "mapt_only_genes.txt")
write(cg.only, "c9orf72_only_genes.txt")
write(gg.only, "grn_only_genes.txt")#
write(cmn.genes, "common_genes.txt")
