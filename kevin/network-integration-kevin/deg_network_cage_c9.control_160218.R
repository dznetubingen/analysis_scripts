#################################################################
## Integrative network generation based on CAGE and miRNA data ##
#################################################################
# first created on 15-02-2018
#
# Network creation based on differentially expressed protein coding and miRNA genes
# Integration of miRNA-target, TF-target and PPI information
#
#
# Currently, for sRNA-seq an CAGE-seq data, the rLog transformed and normalized values are used for correlation calculation
# Both sRNA and CAGE data were analyzed using DESeq2 with Age and Gender as additional covariates included in the model and in the
# rLog transformation
#

# load libs
library(pheatmap)
library(viridis)
library(stringr)
library(igraph)
source("~/scripts/utility_funs.R")


######################################
### Load and modify necessary data ###
######################################

# Load rLog transformed expression values
cage <- read.table("~/rimod/CAGE/analysis/CAGE_rLog_expression_values_dc.contrast_160218.txt", sep="\t", header=T, check.names = F, row.names = 1)
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
# Target filtering
# targets <- targets[targets$TargetScan == 1,] # only consider targets also found in TargetScan


## Data formatting 
# Bring data in suitable format
samples <- as.character(design$sample)
# Extract sample IDs from RNA count table headers
cage.samples <- substr(colnames(cage), 8, 12)
# Only consider mRNA samples that are available in miRNA dataset
cage.keep <- cage.samples %in% design$sample
cage <- cage[,cage.keep]
cage.samples <- cage.samples[cage.keep]
# Bringt data in compatible format and order everything
mir.samples <- substr(colnames(mirna), 8, 12)
mirna.order <- match(cage.samples, mir.samples)
mir.samples <- mir.samples[mirna.order]
mirna <- mirna[mirna.order]
design <- design[match(cage.samples, design$sample),]

## Load DEG files
mirna.deg <- read.table("~/rimod/smallRNA/analysis/deseq_analysi_FTD.C9-NDC_070218/DEGs_sRNA_DESeq_FTD.C9.NDC_070218.txt", sep="\t", header=T, row.names = 1)
cage.deg <- read.table("~/rimod/CAGE/analysis/DEGs_cage_deseq_rimod_C9.control_160218.txt", sep="\t", header=T, row.names=1)

## Load PPI information
ints <- read.table("~/resources/ppi/adjusted_int_table_900_symbols.txt", sep="\t", header=T)

## Load TF-Gene mapping (currently from TRRUST)
org <- read.table("~/resources/TF_interactions/oreganno/ORegAnno_Combined_2016.01.19_hs.tsv", sep="\t", fill = T, stringsAsFactors = F)
org.header <- read.table("~/resources/TF_interactions/oreganno/oreganno_header.tsv", sep="\t", header=T)
colnames(org) <- colnames(org.header)
keep <- c("NEGATIVE OUTCOME", "POSITIVE OUTCOME")
org <- org[org$Outcome %in% keep,] # only keep interactions with known outcome (neutral is probably unknown)
# Don't use miRNA-target interactions
org <- org[!grepl("hsa-miR", org$Regulatory_Element_Symbol),]

## Load TF-miRNA mapping (currently from RegNetwork)
tf.mir <- read.table("~/resources/TF_interactions/RegNetwork/human/human.source", sep="\t", stringsAsFactors = F)
colnames(tf.mir) <- c("symbol_source", "gene_id_source", "target_symbol", "gene_id_target")
tf.mir$regulation <- rep("positive", nrow(tf.mir))

# Just use control and C9orf72 samples
keep <- as.character(design$group %in% c("NDC", "FTD-C9"))
cage <- cage[,keep]
mirna <- mirna[,keep]

############################################
### data loading and formatting complete ###
############################################

###########################################################
### Regulatory correlation analysis                     ###
### Filter regulatory interactions based on correlation ###
###########################################################
cor_method = "spearman"
lfc_cutoff = log2(1.5)
cor_cutoff = 0.5
pval_cutoff = 0.05

## Filter DEGs
cage.deg <- cage.deg[cage.deg$padj <= 0.01,]
cage.deg <- cage.deg[abs(cage.deg$log2FoldChange) >= lfc_cutoff,]
mirna.deg <- mirna.deg[mirna.deg$padj <= 0.01,]
# Only consider genes with sufficient LFC
mirna.sig <- mirna.deg[abs(mirna.deg$log2FoldChange) >= lfc_cutoff,]
cage.sig <- cage.deg[abs(cage.deg$log2FoldChange) >= lfc_cutoff,]

## Data frame to store network in
network <- data.frame(source = "dec", 
                      target = "dec", 
                      interaction = "dec", 
                      intValue = 1)


degenes <- rownames(cage.deg)

#========================================================================#
### TF-gene interactions
for (deg in degenes){
  deg.exp <- as.numeric(cage[rownames(cage) == deg,])
  # acts as TF
  if (deg %in% org$Regulatory_Element_Symbol){
    deg.regs <- org[org$Regulatory_Element_Symbol == deg,]
    deg.regs <- deg.regs[deg.regs$Gene_Symbol %in% degenes,]
    if (nrow(deg.regs) > 0){
      for (i in 1:nrow(deg.regs)){
        b <- deg.regs[i,]$Gene_Symbol
        b.exp <- as.numeric(cage[rownames(cage) == b,])
        ct <- cor.test(deg.exp, b.exp, method = cor_method)
        if (ct$p.value <= pval_cutoff){
          entry <- data.frame(source = deg, target = b, interaction = deg.regs[i,]$Outcome, intValue = ct$estimate)
          network <- rbind(network, entry)
        }
      }
    }
  }
  # regulated by TF
  if(deg %in% org$Gene_Symbol){
    deg.regs <- org[org$Gene_Symbol == deg,]
    deg.regs <- deg.regs[deg.regs$Regulatory_Element_Symbol %in% degenes,]
    if (nrow(deg.regs) > 0){
      for (i in 1:nrow(deg.regs)){
        b <- deg.regs[i,]$Regulatory_Element_Symbol
        b.exp <- as.numeric(cage[rownames(cage) == b,])
        ct <- cor.test(deg.exp, b.exp, method = cor_method)
        if (ct$p.value <= pval_cutoff){
          entry <- data.frame(source = b, target = deg, interaction = deg.regs[i,]$Outcome, intValue = ct$estimate)
          network <- rbind(network, entry)
        }
      }
    }
  }
}

network <- network[-1,]
netdup <- paste(network$source, network$target, sep="")
keep <- !duplicated(netdup)
network <- network[keep,]
# Remove interactions where the action does not match the correlation
keep <- c()
for (i in 1:nrow(network)){
  if (network[i,]$interaction == "POSITIVE OUTCOME"){
    if (network[i,]$intValue > 0){
      keep <- c(keep, i)
    }
  }
  if (network[i,]$interaction == "NEGATIVE OUTCOME"){
    if (network[i,]$intValue < 0 ){
      keep <- c(keep, i)
    }
  }
}
network <- network[keep,]
#========================================================================#



#========================================================================#
#### miRNA gene interactions

for (deg in degenes){
  deg.exp <- as.numeric(cage[rownames(cage) == deg,])
  # acts as TF
  if (deg %in% targets$genesymbol){
    deg.regs <- targets[targets$genesymbol == deg,]
    if (nrow(deg.regs) > 0){
      for (i in 1:nrow(deg.regs)){
        b <- deg.regs[i,]$mirnaid
        b.exp <- as.numeric(mirna[rownames(mirna) == b,])
        ct <- cor.test(deg.exp, b.exp, method = cor_method)
        if (ct$p.value <= pval_cutoff && ct$estimate <= -cor_cutoff){
          entry <- data.frame(source = b, target = deg, interaction = "inhibition", intValue = ct$estimate)
          network <- rbind(network, entry)
        }
      }
    }
  }
}
# remove duplicated edges
netdup <- paste(network$source, network$target, sep="")
keep <- !duplicated(netdup)
network <- network[keep,]

#========================================================================#

#========================================================================#
### Add PPIs
net.ppis <- ints[ints$hgnc_symbol_a %in% degenes,]
net.ppis <- net.ppis[net.ppis$hgnc_symbol %in% degenes,]
net.ppis <- net.ppis[net.ppis$is_directional == "t",]
# remove duplicated edges
netdup <- paste(net.ppis$hgnc_symbol_a, net.ppis$hgnc_symbol)
keep <- !duplicated(netdup)
net.ppis <- net.ppis[keep,]
net.ppis <- data.frame(source = net.ppis$hgnc_symbol_a, target = net.ppis$hgnc_symbol, interaction = net.ppis$action, intValue = net.ppis$score/1000)

network <- rbind(network, net.ppis)

#========================================================================#

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
  else if (v %in% rownames(cage.deg)){
    lfc <- as.numeric(cage.deg[rownames(cage.deg) == v,]$log2FoldChange)
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
  else if (v %in% org$Regulatory_Element_Symbol){
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

write_graph(g, file = "~/rimod/CAGE/analysis/cage_network_analysis/deg_cage_c9.control_network.gml", format = "gml")

