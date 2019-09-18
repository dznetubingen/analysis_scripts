#######################################
# PSEA-like cell type specific differential expression
# on rimod frontal RNA-seq data
#######################################
library(PSEA)
library(stringr)
setwd("~/rimod/RNAseq/analysis/deconvolution/")

# Load expresion values
mat <- read.table("frontal_lengthScaledTPM_counts.txt", sep="\t", row.names=1, header=T)
rownames(mat) <- str_split(rownames(mat), pattern="[.]", simplify = T)[,1]
colnames(mat) <- gsub("X", "", colnames(mat))

# Load deconvolution results
fracs <- read.table("cdn_predictions.txt", sep="\t", row.names = 1, header=T)
rownames(fracs) <- gsub("X", "", rownames(fracs))
fracs$Neurons <- fracs$InNeurons + fracs$ExNeurons
fracs <- fracs[, c(-2, -8)]

# Load MD 
md <- read.table("~/rimod/RNAseq/rnaseq_frontal_md.txt", sep="\t", header=T, stringsAsFactors = F)
md <- md[match(rownames(fracs), md$ids),]
md <- md[-19,]

# keep mapt and control
keep <- as.character(unlist(md['mutated_gene'])) %in% c('MAPT', 'control')
# testing for differential expression
md <- md[keep,]
fracs = fracs[as.character(unlist(md['ids'])),]
mat <- mat[, as.character(unlist(md['ids']))]

# create neural difference
group <- as.character(unlist(md['group']))
group[group == 'case'] <- 1
group[group == 'control'] <- 0
group <- as.numeric(group)


# # Create references
neuron_reference <- fracs$Neurons
astro_reference <- fracs$Astrocytes
oligo_reference <- fracs$Oligodendrocytes
micro_reference <- fracs$Microglia
unknown_reference <- fracs$Unknown
endo_reference <- fracs$Endothelial
opc_reference <- fracs$OPC

neural_diff <- neuron_reference * group

# expression
exp <- as.numeric(mat["ENSG00000186868",])
model2 <- lm(exp ~ neuron_reference + astro_reference + oligo_reference + micro_reference)

crplot(model2, "neuron_reference",  newplot = F)
summary(model2)
