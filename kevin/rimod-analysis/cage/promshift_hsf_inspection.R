library(stringr)

hsf <- read.table("/Users/kevin/dzne/rimod_analysis/hsf1_dependent_transactivation.tsv", sep="\t", header=T)
hsf <- as.character(hsf$MoleculeName)
hsf <- str_split(hsf, pattern=" ", simplify = T)[,2]

mapt <- read.table("/Users/kevin/dzne/rimod_package/analysis/deconvolution/cell_type_specificity/MAPT_cell_composition_filterd_DEGs.txt", sep="\t", header=T)
grn <- read.table("/Users/kevin/dzne/rimod_package/analysis/deconvolution/cell_type_specificity/GRN_cell_composition_filterd_DEGs.txt", sep="\t", header=T)
c9 <- read.table("/Users/kevin/dzne/rimod_package/analysis/deconvolution/cell_type_specificity/C9orf72_cell_composition_filterd_DEGs.txt", sep="\t", header=T)

mapt <- mapt[mapt$hgnc_symbol %in% hsf,]
grn <- grn[grn$hgnc_symbol %in% hsf,]
c9 <- c9[c9$hgnc_symbol %in% hsf,]