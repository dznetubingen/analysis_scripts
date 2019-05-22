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
cage <- read.table("~/rimod/CAGE/analysis/CAGE_rLog_expression_values_080218.txt", sep="\t", header=T, check.names = F, row.names = 1)
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
mirna.deg <- read.table("~/rimod/smallRNA/analysis/deseq_analysis_170118/DEGs_smallRNA_padj0.05.txt", sep="\t", header=T, row.names = 1)
cage.deg <- read.table("~/rimod/CAGE/analysis/DEGs_cage_deseq_rimod_frontal_shrLFC_080218.txt", sep="\t", header=T, row.names=1)

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

### Melt everything down to three data frames
# 1: expression of genes (mRNA + miRNA)
# 2: DE analysis of genes (mRNA + miRNA)
# 3: regulatory information (TF-miR, TF-gene, miRNA-gene)

# # rLog value dataframe for all genes
# rlog.df <- mirna
# colnames(rlog.df) <- colnames(cage)
# rlog.df <- rbind(rlog.df, cage)
# 
# # Regulatory interaction data frame
# reg.df <- tf.mir
# org.regs <- as.character(org$Outcome)
# org.regs[org.regs == "POSITIVE OUTCOME"] <- "positive"
# org.regs[org.regs == "NEGATIVE OUTCOME"] <- "negative"
# org.df <- data.frame(symbol_source = org$Regulatory_Element_Symbol, 
#                      gene_id_source = org$Regulatory_Element_ID, 
#                      target_symbol = org$Gene_Symbol,
#                      gene_id_target = org$Gene_ID,
#                      regulation = org.regs)



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
        if (ct$p.value <= pval_cutoff && ct$estimat <= -cor_cutoff){
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

write_graph(g, file = "~/rimod/CAGE/analysis/cage_network_analysis/test_graph_deg_network.gml", format = "gml")


#################################
### Network Analysis ############
#################################
#
# Remove nodes with degree < 3, as they lack additional information to confirm their validity in the network context
# Detect communities (clusters) using the Multilevel algorithm
# Extract communites and save them for futher analysis
#

g2 <- delete.vertices(g, degree(g) < 3)
g2 <- as.undirected(g2)
l <- layout_with_fr(g2)
ebc <- cluster_louvain(g2)
V(g2)$color = membership(ebc)

plot(g2, vertex.size = 5, edge.arrow.size = .5, layout = l, vertex.label.cex = 0.6, vertex.label.color = "black")

# plot every cluster separately
for (i in 1:length(table(membership(ebc)))){
  gsub <- delete.vertices(g, !membership(ebc) == i)
  plot(gsub, vertex.size = 4,  edge.arrow.size = .5, layout = l, vertex.label.cex = 0.6, vertex.label.color = "black")
}


# Extract genes belonging to the communities
no.com <- length(table(membership(ebc)))  # number of detected communities
mem <- membership(ebc)
for (i in 1:no.com) {
  com <- mem[mem == i]
  com.genes <- names(com)
  file.name <- paste("~/rimod/smallRNA/analysis/network_analysis/community_analysis/community_", i, "_multilevel.txt", sep="")
  write.table(com.genes, file.name, quote=F, row.names = F)
  
  # Additional generate igraph objects
  g.com <- delete.vertices(g2, !V(g2)$name %in% com.genes)
  file.name.gml <- paste("~/rimod/smallRNA/analysis/network_analysis/community_analysis/community_", i, "_multilevel.gml", sep="")
  write_graph(g.com, file.name.gml, format = "gml")
}












### ======================================================================================================================================= ###
####
## Create rLog expression table from miRNA and RNA genes that are in the network
## for non-integrated network inference (e.g. ARACNE)
################
# x_mirna <- mirna
# x_rna <- rna
# colnames(x_rna) <- colnames(x_mirna)
# x_all <- rbind(x_mirna, x_rna)
# net.genes <- c(as.character(network$source), as.character(network$target))
# x_all <- x_all[rownames(x_all) %in% net.genes,]
# write.table(x_all, "~/rimod/smallRNA/analysis/network_analysis/rLog_expression_network_genes.txt", sep="\t", col.names=NA, quote=F)


# ### Read ARACNE inferred network edges
# arac <- read.csv("~/rimod/smallRNA/analysis/network_analysis/ARACNE Inference 2--clustered default  edge.csv")
# arac.ints <- as.character(arac$name)
# arac.g1 <- as.character(sapply(arac.ints, function(x){strsplit(x, split=" ")[[1]][[1]]}))
# arac.g2 <- as.character(sapply(arac.ints, function(x){strsplit(x, split=" ")[[1]][[3]]}))
# arac.gg <- data.frame(gene1 = arac.g1, gene2 = arac.g2)
# 
# # Not considering directionality, find edges that are both found in the 
# # integrated approach and in the ab-initio ARACNe approach
# count <- 0
# cand.edges <- c()
# genes1 <- as.character(arac.gg$gene1)
# genes2 <- as.character(arac.gg$gene2)
# for (i in 1:nrow(arac.gg)){
#   g1 <- genes1[i]
#   g2 <- genes2[i]
#   
#   sub1 <- network[network$source == g1,]
#   sub2 <- network[network$target == g1,]
#   
#   sub1.count <- nrow(sub1[sub1$target == g2,])
#   sub2.count <- nrow(sub2[sub2$source == g2,])
#   
#   count <- count + sub1.count + sub2.count
# 
#   if (sub1.count > 0 || sub2.count > 0){
#     cand.edges <- c(cand.edges, i)
#   }
# }

# Create connections for PyPanda analysis

net.ppi <- network[network$interactionType == "ppi",]
net.tf <- network[!network$interactionType == "ppi",]

net.ppi.pp <- net.ppi[,c(1,2,4)]
net.ppi.tf <- net.tf[,c(1,2,4)]

# write.table(net.ppi.pp, "~/rimod/smallRNA/analysis/network_analysis/pypanda_analysis/network_ppi_data_lfc0.5_pval0.5_pypanda.txt", sep="\t", row.names = F, quote=F)
# write.table(net.ppi.tf, "~/rimod/smallRNA/analysis/network_analysis/pypanda_analysis/network_motif_data_lfc0.5_pval0.5_pypanda.txt", sep="\t", row.names = F, quote=F)


############################
### Compare to results of Paper "Weighted Protein Interaction Network Analysis of Frontotemporal Dementia" (2017)

pn_hubs <- c("COPS5", "ESR1", "HSP90AB1", "STUB1", "EGFR", "FN1", "HSP90AA1", "HSPA8", "PDCD6IP", "TP53", "VCP", "APP", "FSCN1", "GNBL21",
             "HDAC1", "HSPA4", "HTT", "PIN1", "VCAM1", "YWHAZ", "CDK2", "ELAVL1", "EP300", "MCM7", "PML", "RPS3", "TCP1", "TRIM32", "TUBA1A")

verts <- V(g)$name

hub_overlap <- pn_hubs[pn_hubs %in% verts]

##### Create network containing these hub genes
edges <- c()
for (i in 1:nrow(network)) {
  if (network[i,1] %in% hub_overlap || network[i,2] %in% hub_overlap){
    e <- c(as.character(network[i,1]), as.character(network[i,2]))
    edges <- c(edges, e)
  }
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

write_graph(g, file = "~/rimod/smallRNA/analysis/network_analysis/hub_gene_network.gml", format = "gml")
