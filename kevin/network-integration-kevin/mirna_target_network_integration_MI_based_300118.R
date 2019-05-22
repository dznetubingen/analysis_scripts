############################################
## miRNA-mRNA target correlation analysis ##
############################################
# first created on 30-01-2018
# Integration of DE miRNAs with predicted targets using
# - mutual information as measure of relatedness
# - target gene state (DE or not)
# - PPI information
# - planned: proteomics data
#
# Currently, for sRNA-seq an RNA-seq data, the rLog transformed and normalized values are used for correlation calculation
# Both sRNA and mRNA were analyzed using DESeq2 with Age and Gender as additional covariates included in the model and in the
# rLog transformation
# note: for the RNAseq data, gender bias can be found despite this in the MDS plot - but possibly no need for correction necessary
#       if only used for correlation
#

# load libs
library(pheatmap)
library(viridis)
library(stringr)
library(infotheo)
library(igraph)
source("~/scripts/utility_funs.R")


######################################
### Load and modify necessary data ###
######################################

# Load rLog transformed expression values
rna <- read.table("~/rimod/RNAseq/deseq_analysis_220118/rLog_expression_values_RNAseq.txt", sep="\t", header=T, check.names = F, row.names=1)
mirna <- read.table("~/rimod/smallRNA/analysis/deseq_analysis_170118/rLog_expression_values_sRNA.txt", sep="\t", header=T, row.names=1)

# Load proteomics data
prot <- read.table("~/rimod/Proteomics/RIMOD_data_2017-04-27/swath_data_2017-04-27.tsv", sep="\t", header=T, stringsAsFactors = F)
rownames(prot) <- prot$GENE_SYMBOL
prot <- prot[,c(-1,-2)]

# load phenodata
design <- read.table("~/rimod/smallRNA/smallRNA_design_file.txt", sep="\t", header=T)
design$sample <- str_pad(design$sample, 5, side = "left", pad = "0")
# Exclude bad samples
design <- design[!is.na(design$gender),]
design <- design[!design$sample == "10166",] # bad RIN
mirna <- mirna[,substr(colnames(mirna), 8, 12) %in% design$sample]

# Load protein phenodata
prot.md <- read.table("~/rimod/Proteomics/RIMOD_data_2017-04-27/metadata_2017-04-27.tsv", sep="\t", header=T)
prot.md <- prot.md[match(prot.md$Sample.Label, colnames(prot)),]
prot.md$Sample.ID <- str_pad(prot.md$Sample.ID, 5, side ="left", pad ="0")

# Load target prediction file (currently TargetScan)
targets <- read.csv("~/resources/miRNA_target_predictions/mirwalk/miRWalk_miRNA_targets_FTD.NDC.DEGs.csv")
# Target filtering
# targets <- targets[targets$TargetScan == 1,] # only consider targets also found in TargetScan


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
# Only consider proteomics information from samples in the other datasets
prot.md <- prot.md[prot.md$Sample.ID %in% design$sample,]
prot.md <- prot.md[prot.md$Brain.area == "frontal",]
prot <- prot[,colnames(prot) %in% prot.md$Sample.Label]
prot.log <- log2(prot + 0.5)

## Load DEG files
mirna.deg <- read.table("~/rimod/smallRNA/analysis/deseq_analysis_170118/DEGs_smallRNA_padj0.05.txt", sep="\t", header=T, row.names = 1)
rna.deg <- read.table("~/rimod/RNAseq/deseq_analysis_220118/DEGs_RNAseq_padj0.05.txt", sep="\t", header=T, row.names=1)

## Load PPI information
ints <- read.table("~/resources/ppi/adjusted_int_table_900_symbols.txt", sep="\t", header=T)

## Load TF-Gene mapping (currently from TRRUST)
tf.gene <- read.table("~/resources/TF_interactions/TRRUST/trrust_rawdata.human.tsv", sep="\t", stringsAsFactors = F)
colnames(tf.gene) <- c("source", "target", "interaction", "literature")
tf.gene$interaction[tf.gene$interaction == "Unknown"] <- "unkown"
tf.gene$interaction[tf.gene$interaction == "Repression"] <- "inhibition"
tf.gene$interaction[tf.gene$interaction == "Activation"] <- "activation"

## Load TF-miRNA mapping (currently from RegNetwork)
tf.mir <- read.table("~/resources/TF_interactions/RegNetwork/human/human.source", sep="\t", stringsAsFactors = F)
colnames(tf.mir) <- c("source", "n1", "target", "n2")

############################################
### data loading and formatting complete ###
############################################

################################################
### miRNA-target integration based on MI #######
################################################
cor_method = "pearson"
lfc_cutoff = 0.5
pval_cutoff = 0.05
mi_cutoff = 0.5

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
  print(m)
  # Get miRNA counts
  mir.counts <- mirna[rownames(mirna) == m,]
  # Get DE targets
  mir.targets <- as.character(targets[targets$mirnaid == m,]$genesymbol)
  mir.targets <- mir.targets[mir.targets %in% rownames(rna.sig)]
  mir.target.counts <- na.omit(rna[rownames(rna) %in% mir.targets,])
  
  # Only if there are target genes available
  if (length(mir.targets) > 0){
    # Calculate the mutual information
    cand.genes <- c()
    cand.mis <- c()
    mir.counts <- discretize(as.numeric(mir.counts))
    for (r in 1:nrow(mir.target.counts)){
      target.counts <- discretize(as.numeric(mir.target.counts[r,]))
      tmp.mi <- mutinformation(mir.counts, target.counts)
      if (tmp.mi >= mi_cutoff){
        cand.genes <- c(cand.genes, rownames(mir.target.counts)[r])
        cand.mis <- c(cand.mis, tmp.mi)
      }
    }
    
    
    
    # Add candidate miRNA-target connections to network
    if (length(cand.genes) > 0){
      for (i in 1:length(cand.genes)) {
        source.lfc <- mirna.deg[rownames(mirna.deg) == m,]$log2FoldChange
        net.entry <- data.frame(source = m, 
                                target = cand.genes[i], 
                                interaction = "inhibition", 
                                intValue = cand.mis[i], 
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

### miRNA-target analys complete

###########################################
### TF-target correlation analysis ########
###########################################
tf.source <- as.character(tf.gene$source)
tf.target <- as.character(tf.gene$target)
keep <- tf.source %in% network$target | tf.target %in% network$target
tf.gene.flt <- tf.gene[keep,]

# Calculate correlation for all TF-gene pairs
for (i in 1:nrow(tf.gene.flt)) {
  tf <- tf.gene.flt[i,1]
  gene <- tf.gene.flt[i,2]
  # Only proceed if both TF and gene are expressed
  if (tf %in% rownames(rna) && gene %in% rownames(rna)){
    tf.cts <- discretize(rna[rownames(rna) == tf,])
    gene.cts <- discretize(rna[rownames(rna) == gene,])
    tf.mi <- mutinformation(tf.cts, gene.cts)
    # If correlation is significant, add to network
    if (tf.mi >= mi_cutoff){
      source.lfc = 0
      if (tf %in% rownames(rna.deg)){
        source.lfc <- rna.deg[rownames(rna.deg) == tf,]$log2FoldChange
      } 
      net.entry <- data.frame(source = tf, 
                              target = gene, 
                              interaction = as.character(tf.gene.flt[i,3]), 
                              intValue = tf.mi, 
                              lfcSource = source.lfc, 
                              interactionType = "tf.gene",
                              sourceType = "gene")
      network <- rbind(network, net.entry)
    }
  }
}

### TF-target correlation complete


##########################################
### Add PPI information to the network ###
##########################################

# Extract relevant PPIs
net.targets <- as.character(network$target)
ints.net <- ints[ints$hgnc_symbol_a %in% net.targets,]
ints.net <- ints.net[ints.net$hgnc_symbol %in% net.targets,]
# Find LFCs for source nodes
source.prots <- ints.net$hgnc_symbol_a
lfcs <- c()
for (i in 1:length(source.prots)){
  prot <- source.prots[i]
  lfc = 0
  if (prot %in% rownames(rna.deg)){
    lfc <- rna.deg[rownames(rna.deg) == prot,]$log2FoldChange
  }
  lfcs <- c(lfcs, lfc)
}
# Create PPI data frame
ppi.df <- data.frame(source = ints.net$hgnc_symbol_a, 
                     target = ints.net$hgnc_symbol, 
                     interaction = ints.net$action, 
                     intValue = ints.net$score/1000,
                     lfcSource = lfcs,
                     interactionType = "ppi",
                     sourceType = "Protein")
# Remove double-edges to avoid bias while clustering
keep <- !duplicated(paste(ppi.df$source, ppi.df$target, sep=""))
ppi.df <- ppi.df[keep,]
# Add to network
network <- rbind(network, ppi.df)

#########################################
### Add TF-miRNA maapings (RegNetwork) ##
#########################################
# Extract TFs that target miRNAs of the network
net.mirs <- as.character(network$source)
net.mirs <- net.mirs[grepl("hsa-", net.mirs)]
net.mirs <- net.mirs[!duplicated(net.mirs)]
tf.mir.sub <- tf.mir[tf.mir$target %in% net.mirs,]
tf.mir.sub <- tf.mir.sub[tf.mir.sub$source %in% rownames(rna.deg),]
# Get LFCs
lfcs <- c()
tfs <- as.character(tf.mir.sub$source)
for (i in 1:length(tfs)){
  tf <- tfs[i]
  lfc = 0
  if (tf %in% rownames(rna.deg)){
    lfc <- rna.deg[rownames(rna.deg) == tf,]$log2FoldChange
  }
  lfcs <- c(lfcs, lfc)
}

tf.mir.df <- data.frame(source = tf.mir.sub$source,
                        target = tf.mir.sub$target,
                        interaction = "unkown",
                        intValue = rep(0.5, nrow(tf.mir.sub)),
                        lfcSource = lfcs,
                        interactionType = "tf.mir",
                        sourceType = "TF")
# Add to network
network <- rbind(network, tf.mir.df)

# Adjust rownames
rownames(network) <- c(1:nrow(network))

# Save files
write.table(network, "~/rimod/smallRNA/analysis/network_analysis/miRNA_gene_network_pval0.05_300118.txt", sep="\t", row.names=F, quote=F)


#########################################
### Include Proteomics Data #############
#########################################

# # Limma analysis of proteomic data
# G <- as.character(prot.md$Disease.Code)
# G[grepl("FTD", G)] <- "FTD"
# G <- as.factor(G)
# age <- as.numeric(prot.md$Age)
# gender <- as.factor(prot.md$Gender)
# pdesign <- model.matrix(~ G + age + gender)
# colnames(pdesign) <- c("FTD", "NDC", "age", "F", "M")
# cont.mat <- makeContrasts(contrasts = "FTD-NDC", levels=pdesign)
# fit <- lmFit(prot.log, design = pdesign)
# fit.cont <- contrasts.fit(fit, contrasts = cont.mat)
# fit2 <- eBayes(fit.cont)
# res <- topTable(fit2, coef = 1)


#############################
### iGraph analysis #########
#############################

# Convert network to igraph format
edges <- c()
for (i in 1:nrow(network)) {
  e <- c(as.character(network[i,1]), as.character(network[i,2]))
  edges <- c(edges, e)
}
g <- graph(edges=edges, directed = T)

l <- layout_with_fr(g)
plot(g, vertex.color = "gold", vertex.size = 5, edge.arrow.size = .5, layout = l, 
     vertex.label.cex = 0.6, vertex.label.color = "black")


### Collect information on vertices and add to the graph
# Collect fold changes
vnames <- V(g)$name
net.lfcs <- c()
for (i in 1:length(vnames)){
  lfc <- 0
  v <- vnames[i]
  if (grepl("hsa-", v)) { # check if miRNA
    lfc <- as.numeric(mirna.deg[rownames(mirna.deg) == v,]$log2FoldChange)
  }
  else if (v %in% rownames(rna.deg)){
    lfc <- as.numeric(rna.deg[rownames(rna.deg) == v,]$log2FoldChange)
  }
  
  net.lfcs <- c(net.lfcs, lfc)
}

# Assign type
types = c()
for (i in 1:length(vnames)){
  v <- vnames[i]
  type = ""
  if (grepl("hsa-",v)){
    type = "miRNA"
  }
  else if (v %in% tf.gene$source){
    type = "TF"
  }
  else {
    type = "Gene"
  }
  types <- c(types, type)
}

# Assign to graph object
V(g)$lfc <- net.lfcs
V(g)$type <- types

### Collect information on edges and add to the graph
E(g)$intType <- as.character(network$interaction)
E(g)$intValue <- network$intValue
E(g)$intType2 <- as.character(network$interactionType)

# Test GML format
write_graph(g, "~/rimod/smallRNA/analysis/network_analysis/miRNA_target_TF_network_pval0.05_300118.gml", format = "gml")

# Do some igraph stuff

g2 <- delete.vertices(g, degree(g) < 3)
l <- layout.fruchterman.reingold(g2)
ebc <- walktrap.community(g2)
V(g2)$color = membership(ebc)

plot(g2, vertex.size = 5, edge.arrow.size = .5, layout = l, vertex.label.cex = 0.6, vertex.label.color = "black")
