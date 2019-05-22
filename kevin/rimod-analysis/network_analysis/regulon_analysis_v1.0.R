############################################################
## Regulon analysis based on TF-gene network
##
############################################################

# load libs
library(igraph)
library(biomaRt)
# Parse args
args <- commandArgs(trailingOnly = TRUE)

net_file <- args[1]
deg_file <- args[2]
lfc <- as.numeric(args[3])
outname <- args[4]

# FOR RESTING
#setwd("~/rimod/CAGE/cage_analysis/regulon_down_analysis_300418//mapt/")
#net_file <- "mapt_tf_target_mapping.txt"
#deg_file <- "~/rimod/CAGE/cage_analysis/regulon_analysis_250418/CAGE_deseq_analysis_2018-04-26_genewise/deseq_result_mapt.ndc_2018-04-26_14.22.04.txt"
#lfc <- -0.1
#utname <- "mapt"

# Define cutoffs for up-regulation
pval = 0.05
lfc_cut = lfc
jac.cutoff <- 0.1


# Load network
print(net_file)
network <- read.table(net_file, sep="\t", header = T, stringsAsFactors = F)
deg <- read.table(deg_file, sep="\t", header=T, row.names=1)

## Up-regulated TFs, merge with annotation
deg.sig <- deg[deg$padj <= pval,]
if (lfc_cut > 0){
  deg.sig <- deg.sig[deg.sig$log2FoldChange >= lfc_cut,]
} else {
  deg.sig <- deg.sig[deg.sig$log2FoldChange <= lfc_cut,]
}

## Use biomaRt to get all the HGNC names and ENSEMBL IDs
ensgenes <- rownames(deg.sig)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = ensgenes, mart = ensembl)

deg.sig <- merge(deg.sig, bm, by.x = "row.names", by.y="ensembl_gene_id")

# Get symbols of TFs in network
network$TF <- toupper(network$TF)
tfs <- network$TF
tfs <- tfs[!duplicated(tfs)]
# Get signfificantly up-regulated TFs
tfs.up <- c() # significant tfs 
for (tf in tfs) {
  if (grepl("::", tf)){
    tf1 <- strsplit(tf, split="::")[[1]][[1]]
    tf2 <- strsplit(tf, split="::")[[1]][[2]]
    if (tf1 %in% deg.sig$hgnc_symbol && tf2 %in% deg.sig$hgnc_symbol){
      tfs.up <- c(tfs.up, tf)
    }
  } else if(tf %in% deg.sig$hgnc_symbol){
    tfs.up <- c(tfs.up, tf)
  }
}

#################################################
### Creat network for all TFs
#################################################

##############################
### Define regulons
##############################
regs <- levels(as.factor(as.character(network$TF)))
regulons <- list()
for (i in 1:length(regs)) {
  reg <- regs[i]
  reg.targets <- as.character(network[network$TF == reg,]$ensembl_gene_id)
  reg.targets <- reg.targets[!duplicated(reg.targets)]    # remove duplicate targets
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
if (length(edges) > 0){
  # Create graph 
  g <- graph(edges = edges)
  E(g)$weight <- jacs
  
  # Calculate sizes of Regulons
  reg.size <- sapply(regulons, length)
  reg.size <- reg.size[names(reg.size) %in% V(g)$name]
  reg.size <- reg.size[match(V(g)$name, names(reg.size))]
  V(g)$reg.size <- reg.size
  
  write_graph(g, file = paste(outname, "_regulon_graph_cyto.gml", sep=""), format = "gml")
  save(regulons, file=paste(outname, "_regulons.RData", sep=""))
  saveRDS(regulons, file=paste(outname, "_regulons.rds", sep=""))
}else {
  print("No Regulon network found")
}

##########################################################
### Create regulons for sig. up-regulated TFs
##########################################################
if (length(tfs.up) > 0){
  # Subset network
  network.sig <- network[network$TF %in% tfs.up,]
  
  # Lower jac cutoff for smaller network
  jac.cutoff <- 0.01
  
  regs <- levels(as.factor(as.character(network.sig$TF)))
  regulons <- list()
  for (i in 1:length(regs)) {
    reg <- regs[i]
    reg.targets <- as.character(network.sig[network.sig$TF == reg,]$ensembl_gene_id)
    reg.targets <- reg.targets[!duplicated(reg.targets)]    # remove duplicate targets
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
  if (length(edges) > 0){
    # Create graph 
    g <- graph(edges = edges)
    E(g)$weight <- jacs
    
    # Calculate sizes of Regulons
    reg.size <- sapply(regulons, length)
    reg.size <- reg.size[names(reg.size) %in% V(g)$name]
    reg.size <- reg.size[match(V(g)$name, names(reg.size))]
    V(g)$reg.size <- reg.size
    
    # Get LFCs
    tf.deg <- deg.sig[deg.sig$hgnc_symbol %in% V(g)$name,]
    tf.deg <- tf.deg[!duplicated(tf.deg$hgnc_symbol),]
    tf.deg <- tf.deg[match(V(g)$name, tf.deg$hgnc_symbol),]
    V(g)$lfc <- tf.deg$log2FoldChange
    
    write_graph(g, file = paste(outname, "_sig_regulon_graph_cyto.gml", sep=""), format = "gml")
    save(regulons, file=paste(outname, "_sig_regulons.RData", sep=""))
    saveRDS(regulons, file=paste(outname, "_sig_regulons.rds", sep=""))
  }else {
    print("No Regulon network found")
  }
} else {
  print("No differentially expressed TFs found")
}
