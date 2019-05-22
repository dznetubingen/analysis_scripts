###=============================================================================###
### miRNA-target regulon integration 
###
# Regulon based approach to integrate miRNA expression with gene expression
#
# First regulons are defined as miRNA target clusters with significant negative correlation
# Then overlap of these regulons is compared with the state of the regulators and their regulons
#
#
#

# load libs
library(pheatmap)
library(viridis)
library(stringr)
library(igraph)
source("~/scripts/utility_funs.R")

setwd("~/rimod/smallRNA/analysis/network_analysis/regulon_analysis/")

######################################
### Load and modify necessary data ###
######################################
# Load rLog transformed expression values
rna <- read.table("~/rimod/RNAseq/deseq_analysis_220118/rLog_expression_values_RNAseq.txt", sep="\t", header=T, check.names = F, row.names=1)
mirna <- read.table("~/rimod/smallRNA/analysis/deseq_analysis_170118/rLog_expression_values_sRNA.txt", sep="\t", header=T, row.names=1)
# load phenodata
design <- read.table("~/rimod/smallRNA/smallRNA_design_file.txt", sep="\t", header=T)
design$sample <- str_pad(design$sample, 5, side = "left", pad = "0")
# Exclude bad samples
design <- design[!is.na(design$gender),]
design <- design[!design$sample == "10166",] # bad RIN
mirna <- mirna[,substr(colnames(mirna), 8, 12) %in% design$sample]
# Load target prediction file (currently TargetScan)
targets <- read.csv("~/resources/miRNA_target_predictions/mirwalk/miRWalk_miRNA_targets_FTD.NDC.DEGs.csv")
## Data formatting 
# Bring data in suitable format
samples <- as.character(design$ids)
# Extract sample IDs from RNA count table headers
rna.samples <- substr(colnames(rna), 1, 5)
# Only consider mRNA samples that are available in miRNA dataset
rna.keep <- rna.samples %in% design$sample
rna <- rna[,rna.keep]
rna.samples <- rna.samples[rna.keep]
# Bringt data in compatible format and order everything
mir.samples <- substr(colnames(mirna), 8, 12)
mirna.order <- match(rna.samples, mir.samples)
mir.samples <- mir.samples[mirna.order]
mirna <- mirna[mirna.order]
design <- design[match(rna.samples, design$sample),]
## Load DEG files
mirna.deg <- read.table("~/rimod/smallRNA/analysis/deseq_analysis_170118/DEGs_smallRNA_padj0.05.txt", sep="\t", header=T, row.names = 1)
rna.deg <- read.table("~/rimod/RNAseq/deseq_analysis_220118/DEGs_RNAseq_padj0.05.txt", sep="\t", header=T, row.names=1)
## Load PPI information
ints <- read.table("~/resources/ppi/adjusted_int_table_900_symbols.txt", sep="\t", header=T)
############################################
### data loading and formatting complete ###
############################################

################################################
### miRNA-target Correlation Analysis ##########
################################################
# Refine regulons based on correlation of miRNAs and their targets
cor_method = "pearson"
lfc_cutoff = 0.5
cor_cutoff = -0.4
pval_cutoff = 0.05
## Data frame to store network in
network <- data.frame(source = "dec", 
                      target = "dec", 
                      interaction = "dec", 
                      intValue = 1, 
                      lfcSource = 1, 
                      interactionType = "dec",
                      sourceType = "dec")
# Only consider genes with sufficient LFC
mirna.sig <- mirna.deg[abs(mirna.deg$log2FoldChange) >= lfc_cutoff,]
rna.sig <- rna.deg[abs(rna.deg$log2FoldChange) >= lfc_cutoff,]
## For every significant miRNA, calculate the correlation with its
## predicted targets
for (m in rownames(mirna.sig)){
  # Get miRNA counts
  mir.counts <- mirna[rownames(mirna) == m,]
  # Get DE targets
  mir.targets <- as.character(targets[targets$mirnaid == m,]$genesymbol)
  mir.targets <- mir.targets[mir.targets %in% rownames(rna.sig)]
  mir.target.counts <- na.omit(rna[rownames(rna) %in% mir.targets,])
  
  if (length(mir.targets) > 0){
    
    # Calculate the correlation
    cand.genes <- c()
    cand.cors <- c()
    for (r in 1:nrow(mir.target.counts)){
      tmp.cor <- cor.test(as.numeric(mir.counts), as.numeric(mir.target.counts[r,]), method=cor_method)
      if (tmp.cor$p.value <= pval_cutoff && tmp.cor$estimate <= cor_cutoff){
        cand.genes <- c(cand.genes, rownames(mir.target.counts)[r])
        cand.cors <- c(cand.cors, tmp.cor$estimate)
      }
    }
    
    
    # Add candidate miRNA-target connections to network
    if (length(cand.genes) > 0){
      for (i in 1:length(cand.genes)) {
        source.lfc <- mirna.deg[rownames(mirna.deg) == m,]$log2FoldChange
        net.entry <- data.frame(source = m, 
                                target = cand.genes[i], 
                                interaction = "inhibition", 
                                intValue = cand.cors[i], 
                                lfcSource = source.lfc, 
                                interactionType = "miRNA.target",
                                sourceType = "miRNA")
        network <- rbind(network, net.entry)
      }
    }
  }
  
}
# Remove decoy row
network <- network[-1,]

##############################
### Define regulons
##############################
regs <- levels(as.factor(as.character(network$source)))
regulons <- list()
for (i in 1:length(regs)) {
  reg <- regs[i]
  reg.targets <- as.character(network[network$source == reg,]$target)
  regulons$tmp <- reg.targets
  names(regulons)[i] <- reg
}

#############
## Create Network with regulons as vertices and edges as Jaccard distance between edges

# Calculate Jaccard coefficient
jaccard <- function(a, b){
  int <- length(intersect(a,b))
  jac <- int / (length(a) + length(b) -int)
  return(jac)
}

# Create matrix to store values in 
no.regs <- length(regulons)
jac.mat <- matrix(0, no.regs, no.regs)
colnames(jac.mat) <- names(regulons)
rownames(jac.mat) <- names(regulons)
# Calculate jaccard for every pair
for (i in 1:no.regs){
  a <- regulons[[i]]
  for (j in 1:no.regs){
    if (j != i){
      if (jac.mat[i,j] == 0){
        b = regulons[[j]]
        jac <- jaccard(a,b)
        jac.mat[i,j] <- jac
        jac.mat[j,i] <- jac
      }
    }
  }
}

# Create edges for nodes with a positive jaccard score
jac.cutoff <- 0.05
edges <- c()
jacs <- c()
for (i in 1:no.regs){
  for (j in (i+1):no.regs){
    if (j <= no.regs){
      jac <- jac.mat[i,j]
      if (jac > jac.cutoff){
        e <- c(names(regulons)[i], names(regulons[j]))
        edges <- c(edges, e)
        jacs <- c(jacs, jac)
      }
    }
  }
}


# Create graph 
g <- graph(edges = edges)
E(g)$weight <- jacs
# Sizes of regulons and LFC of regulators
size <- c()
lfcs <- c()
new.regs <- regulons[names(regulons) %in% V(g)$name]
new.regs <- new.regs[match(V(g)$name, names(new.regs))]
for (i in 1:length(new.regs)){
  reg.name <- names(new.regs)[i]
    size <- c(size, length(new.regs[[i]]))
    lfc <- as.numeric(mirna.deg[rownames(mirna.deg) == reg.name,]$log2FoldChange)
    lfcs <- c(lfcs, lfc)
}
V(g)$size <- size
V(g)$lfc <- lfcs
# LFC of regulators

#write_graph(g, file = "regulon_graph.gml", format = "gml")

########## Add TF-miRNA(regulon) information to network (create new network)
## Load TF-miRNA mapping (currently from RegNetwork)
tf.mir <- read.table("~/resources/TF_interactions/RegNetwork/human/human.source", sep="\t", stringsAsFactors = F)
colnames(tf.mir) <- c("source", "n1", "target", "n2")
# Now load proteomics data and only consider TFs that are measured in the brain (frontal and temporal for now)
prot <- read.table("~/rimod/Proteomics/RIMOD_data_2017-04-27/swath_data_2017-04-27.tsv", sep="\t", header=T)
prot.ids <- prot[,1:2]
prot <- prot[,c(-1,-2)]
rownames(prot) <- prot.ids$GENE_SYMBOL
reg.names <- names(new.regs)
tf.mir <- tf.mir[tf.mir$source %in% rownames(prot),]
tf.mir <- tf.mir[tf.mir$target %in% reg.names,]
for (i in 1:nrow(tf.mir)){
  a <- as.character(tf.mir[i,1])
  b <- as.character(tf.mir[i,3])
  e <- c(a,b)
  edges <- c(edges, e)
}
tf_regulon_graph <- graph(edges)
write_graph(tf_regulon_graph, file = "~/rimod/smallRNA/analysis/network_analysis/regulon_analysis/TF_regulon_graph.gml", format = "gml")




####################################################
### Collect common targets of connected regulons ###
####################################################
reg.names <- names(new.regs)
mirna.deg.reg <- mirna.deg[rownames(mirna.deg) %in% reg.names,]
mirna.deg.reg.up <- mirna.deg.reg[mirna.deg.reg$log2FoldChange > 0,]
mirna.deg.reg.down <- mirna.deg.reg[mirna.deg.reg$log2FoldChange < 0,]
reg.up <- new.regs[reg.names %in% rownames(mirna.deg.reg.up)]
reg.down <- new.regs[reg.names %in% rownames(mirna.deg.reg.down)]

reg.up.genes <- table(as.character(unlist(reg.up)))
reg.down.genes <- table(as.character(unlist(reg.down)))

reg.up.genes <- reg.up.genes[reg.up.genes > 1]
reg.down.genes <- reg.down.genes[reg.down.genes > 1]

# write.table(as.character(names(reg.up.genes)), "upMir_common_targets.txt", quote=F, row.names = F)
# write.table(as.character(names(reg.down.genes)), "downMir_common_targets.txt", quote=F, row.names = F)


###====================================================================================================================###
# Based on the results from the Regulon network construction                                                           ###
# create a network from the Regulon target genes (i.e. genes that are possibly regulated by at least 2 miRNAs)         ### 
# This network can be extended by some degree to investigate downstream effects                                        ###
###====================================================================================================================###

# Only consider co-regulated genes in the network object
coreg_genes <- c(names(reg.up.genes), names(reg.down.genes))
nw.coreg <- network[network$target %in% coreg_genes,]

coreg_edges <- c()
for (i in 1:nrow(nw.coreg)){
  s = as.character(nw.coreg[i,1])
  t = as.character(nw.coreg[i,2])
  e = c(s, t)
  coreg_edges <- c(coreg_edges, e)
}
# Create graph
coreg.graph <- graph(coreg_edges)

# Collect fold changes
lfcs <- c()
verts <- V(coreg.graph)$name
for (i in 1:length(verts)){
  v <- verts[i]
  if (grepl("hsa", v)){
    lfc <- as.numeric(mirna.deg[rownames(mirna.deg) == v,]$log2FoldChange)
  }
  else {
    lfc <- as.numeric(rna.deg[rownames(rna.deg) == v,]$log2FoldChange)
  }

  lfcs <- c(lfcs, lfc)
}
V(coreg.graph)$lfc <- lfcs

write_graph(coreg.graph, file = "~/rimod/smallRNA/analysis/network_analysis/regulon_analysis/regulon_target_graph.gml", format = "gml")
