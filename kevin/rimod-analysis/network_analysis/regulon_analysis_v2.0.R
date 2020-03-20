###################################
## miRNA and TF regulon analysis
###################################
library(biomaRt)
library(igraph)
library(stringr)
setwd("~/rimod/integrative_analysis/regulon_analysis/")


# Calculate Jaccard coefficient
jaccard <- function(a, b){
  int <- length(intersect(a,b))
  jac <- int / length(union(a, b))
  return(jac)
}
###
# Calculate percentage of element in a that
# are also in b
percent <- function(a, b){
  alen <- length(a)
  ilen <- length(intersect(a, b))
  pct <- ilen / alen
  return(pct)
}

####
# Generate miRNA and TF regulons
# for a given FTD-group
generateRegulon <- function(deg, mir.targets, tf, outname){
  
  ###
  # Generate miRNA regulons
  ###
  mir.genes <- as.character(unique(mir.targets$mirna))
  
  mir.regulons <- list()
  mir.lfc <- c()
  for (mir in mir.genes) {
    targets <- as.character(mir.targets[mir.targets$mirna == mir,]$targets)
    target.degs <- deg[deg$hgnc_symbol %in% targets,]
    lfc <- mean(target.degs$log2FoldChange)
    mir.lfc <- c(mir.lfc, lfc)
    mir.regulons$tmp <- targets
    names(mir.regulons)[length(mir.regulons)] <- mir
  }
  
  df <- data.frame(regulator = mir.genes, targetLfc = mir.lfc)
  
  # keep only miRNAs with DEGs among targets
  keep <- !is.na(df$targetLfc)
  df <- df[keep,]
  mir.regulons <- mir.regulons[keep]
  
  ##======================================================#
  
  ####
  # Create TF-gene regulons
  ####
  tfs <- as.character(tf$TF)
  tf.regulons <- list()
  tf.lfc <- c()
  
  for (t in tfs) {
    targets <- as.character(tf[tf$TF == t,]$Overlapping_Genes)
    targets <- as.character(str_split(targets, pattern=",", simplify = T))
    target.degs <- deg[deg$hgnc_symbol %in% targets,]
    lfc <- mean(target.degs$log2FoldChange)
    tf.lfc <- c(tf.lfc, lfc)
    tf.regulons$tmp <- targets
    names(tf.regulons)[length(tf.regulons)] <- t
  }
  tf.df <- data.frame(regulator = tfs, targetLfc = tf.lfc)
  
  # append regulons and dfs
  df <- rbind(df, tf.df)
  regulons <- c(mir.regulons, tf.regulons)
  
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
        if (jac > 0.04){
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
  
  # Calculate sizes of Regulons
  reg.size <- sapply(regulons, length)
  reg.size <- reg.size[names(reg.size) %in% V(g)$name]
  reg.size <- reg.size[match(V(g)$name, names(reg.size))]
  V(g)$reg.size <- reg.size
  # regulon fold change
  df <- df[match(V(g)$name, df$regulator),]
  V(g)$lfc <- df$targetLfc
  
  
  
  ###
  # Add membrane trafficking overlap
  ubi <- read.table("~/rimod/integrative_analysis/membrane_trafficking/membrane_trafficking_reactome.tsv", sep="\t", header=T)
  genes <- str_split(ubi$MoleculeName, pattern=" ", simplify = T)[,2]
  ubi$name <- genes
  

  keep <- names(regulons) %in% V(g)$name
  regulons <- regulons[keep]
  regulons <- regulons[match(V(g)$name, names(regulons))]
  
  pct_list <- c()
  inter_list <- c()
  for (r in regulons) {
    j <- percent(r, genes)
    inter <- length(intersect(r, genes))
    if (inter == 0){j <- 0}
    inter_list <- c(inter_list, inter)
    pct_list <- c(pct_list, j)
  }
  print(pct_list)
  V(g)$mtpct <- pct_list
  V(g)$mtinter <- inter_list
  
  
  write_graph(g, file = paste(outname, "_regulon_graph_cyto.gml", sep=""), format = "gml")
  return(regulons)
}
#=============================================================================================#


####
# MAPT
####
outname = "mapt"
# miRNA-targets
mir.targets <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/MAPT_miRNA_target_edge_table.txt", sep="\t", header=T)

# TF targets
tf.up <- read.table("~/rimod/integrative_analysis/tf_cage_rnaseq/MAPT_common_TFs_up.txt", sep="\t", header=T)
tf.down <- read.table("~/rimod/integrative_analysis/tf_cage_rnaseq/MAPT_common_TFs_down.txt", sep="\t", header=T)
tf <- rbind(tf.up, tf.down)

# DEGs
#deg <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_result_mapt.ndc_fro_2019-10-23_13.33.11.txt",
#                  sep="\t", header=T, row.names=1)
#deg <- deg[deg$padj <= 0.05,]
#ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
#bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values = rownames(deg), mart=ensembl)
#deg <- merge(deg,  bm, by.x="row.names", by.y="ensembl_gene_id")

deg <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/MAPT_cell_composition_filterd_DEGs.txt", sep="\t", header=T, row.names=2)
deg <- deg[,-1]


# subset TFs
tf <- tf[tf$TF %in% deg$hgnc_symbol,]

# Generate regulon
mapt.regulons <- generateRegulon(deg=deg, mir.targets=mir.targets, tf=tf, outname=outname)
#==================================================#


####
# FTD-GRN
#####
outname = "grn"
# miRNA-targets
mir.targets <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/GRN_miRNA_target_edge_table.txt", sep="\t", header=T)

# TF targets
tf_cutoff <- 20
tf.up <- read.table("~/rimod/integrative_analysis/tf_cage_rnaseq/GRN_common_TFs_up.txt", sep="\t", header=T)
tf.down <- read.table("~/rimod/integrative_analysis/tf_cage_rnaseq/GRN_common_TFs_down.txt", sep="\t", header=T)
tf <- rbind(tf.up, tf.down)

# DEGs
deg <- read.table("~/rimod/RNAseq/analysis/deconvolution/cell_type_specificity/GRN_cell_composition_filterd_DEGs.txt", sep="\t", header=T, row.names=2)
deg <- deg[,-1]


# subset TFs
tf <- tf[tf$TF %in% deg$hgnc_symbol,]

# Generate regulon
grn.regulons <- generateRegulon(deg=deg, mir.targets=mir.targets, tf=tf, outname=outname)
#==================================================#



####
# FTD-C9orf72
#####
### -> tried it, no network possible

####
# Check out overlap with membrane trafficking
####
###
# Add membrane trafficking overlap
memb <- read.table("~/rimod/integrative_analysis/membrane_trafficking/membrane_trafficking_reactome.tsv", sep="\t", header=T)
memb_genes <- str_split(memb$MoleculeName, pattern=" ", simplify = T)[,2]
memb$name <- memb_genes

# color palette
mypal <- brewer.pal(3, "Dark2")

# mapt overlap
mapt_ovl <- c()
mapt_pct <- c()
for (r in mapt.regulons) {
  ovl <- length(intersect(r, memb_genes))
  pct <- percent(r, memb_genes)
  mapt_ovl <- c(mapt_ovl, ovl)
  mapt_pct <- c(mapt_pct, pct)
}

df <- data.frame(Regulon = names(mapt.regulons), PctOvl = mapt_pct, Ovl = mapt_ovl)
df <- df[order(df$PctOvl),]
df <- df[order(df$Ovl),]
df$Regulon <- factor(df$Regulon, levels=as.character(df$Regulon))
keep <- df$Ovl > 3
df <- df[keep,]

p <- ggplot(df, aes(x=Regulon, y=Ovl)) + 
  geom_bar(stat = "identity", fill=mypal[1]) + 
  theme_minimal(base_size = 10) +
  coord_flip() +
  scale_fill_brewer(palette = "Dark2")
p
ggsave(filename = "mapt_regulon_membraneTrafficking_overlap.png", width=2, height=3.4)


# grn
grn_ovl <- c()
grn_pct <- c()
for (r in grn.regulons) {
  ovl <- length(intersect(r, memb_genes))
  pct <- percent(r, memb_genes)
  grn_ovl <- c(grn_ovl, ovl)
  grn_pct <- c(grn_pct, pct)
}

df <- data.frame(Regulon = names(grn.regulons), PctOvl = grn_pct, Ovl = grn_ovl)
df <- df[order(df$PctOvl),]
df <- df[order(df$Ovl),]
df$Regulon <- factor(df$Regulon, levels=as.character(df$Regulon))
keep <- df$Ovl > 3
df <- df[keep,]

p <- ggplot(df, aes(x=Regulon, y=Ovl)) + 
  geom_bar(stat = "identity", fill = mypal[2]) + 
  theme_minimal(base_size = 10) +
  coord_flip() +
  scale_fill_brewer(palette = "Dark2")
p
ggsave(filename = "grn_regulon_membraneTrafficking_overlap.png", width=2, height=3.4)



