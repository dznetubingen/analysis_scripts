##################################################
## Create PPI network from annotated CAGE peaks ##
##
## PPIs from IntAct and BIOGRID currently
## miRNA-interactions: not implemented yet
##
## Author:  Kevin Menden
## Date:    20.04.2018  
##################################################

# load libs
library(biomaRt)
library(igraph)


# parse command line arguments

args <- commandArgs(trailingOnly = TRUE)

annoted_peaks_file = args[1]   # the annotated peak file (in CSV format)
outname = args[2] # prefix for outfiles
# output directory (optional)
outdir = ""
if (length(args) == 3){
  outdir = args[3]
}

# Parse the peak file and extract transcripts
annot <- read.csv(annoted_peaks_file, row.names = 1)
colnames(annot)[1] = "peak_id"
transcripts <- as.character(annot$Nearest.PromoterID)
transcripts <- as.character(sapply(transcripts, function(x){strsplit(x, split="[.]")[[1]][[1]]}))
# get biomart mapping for DEGs
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol"), filters = "ensembl_transcript_id", values = transcripts, mart = ensembl)
ensGenes = as.character(bm$ensembl_gene_id)

##########
### Load IntAct PPIs and subset
print("Loading IntAct PPIs ...")
intact <- read.table("~/resources/ppi/IntAct/intact/hs.intact.sub.ensembl.txt", sep="\t", header = T, fill = T)
intact <- intact[,c(1,2,3,4,7,8,9,10),]
intact <- intact[intact$ensgA %in% ensGenes,]
intact <- intact[intact$ensgB %in% ensGenes,]

##########
### Load Biogrid PPIs
print("Loading BioGrid PPIs ...")
biog <- read.csv("~/resources/ppi/BIOGRID/BIOGRID-ALL-3.4.159.tab2.csv", row.names = 1)
biog <- biog[biog$Organism.Interactor.A == "9606",]
biog <- biog[biog$Organism.Interactor.B == "9606",]
# Mapt to ENSEMBL_GENE_ID
biogrid_genes <- c(as.character(biog$Entrez.Gene.Interactor.A), as.character(biog$Entrez.Gene.Interactor.B))
biogrid_genes <- biogrid_genes[!duplicated(biogrid_genes)]
bm <- getBM(attributes = c("ensembl_gene_id", "entrezgene"), filters = "entrezgene", values = biogrid_genes, mart = ensembl)
mdata <- merge(bm, biog, by.x="entrezgene", by.y="Entrez.Gene.Interactor.A")
colnames(mdata)[2] <- "ensembl_gene_id_A"
mdata <- merge(bm, mdata, by.x="entrezgene", by.y="Entrez.Gene.Interactor.B")
mdata <- mdata[,c(-1,-3)]
colnames(mdata)[1] <- "ensembl_gene_id_B"
biogrid_ints <- mdata[mdata$ensembl_gene_id_B %in% ensGenes,]
biogrid_ints <- biogrid_ints[biogrid_ints$ensembl_gene_id_A %in% ensGenes,]

############
### Load TF-gene network
network <- read.table("~/rimod/CAGE/cage_analysis/tf_enrichment_pipeline/c9orf72/c9orf72_tf_target_mapping.txt", sep="\t", header=T)

#######
### Create network the network
#######
print("Creating network ...")
edges <- c()
# Add PPIs from IntAct
for (i in 1:nrow(intact)){
  e <- c(as.character(intact[i,5]), as.character(intact[i,6]))
  edges <- c(edges, e)
}
# Add PPIs from BIOGRID
for (i in 1:nrow(biogrid_ints)) {
  e <- c(as.character(biogrid_ints[i,2]), as.character(biogrid_ints[i,1]))
  edges <- c(edges, e)
}
for (i in 1:nrow(network)){
  e <- c(as.character(network[i,3]), as.character(network[i,5]))
  edges <- c(edges, e)
}

g <- graph(edges=edges)

# Assign type
vnames <- V(g)$name
types = c()
for (i in 1:length(vnames)){
  v <- vnames[i]
  type = ""
  if (grepl("ENS",v)){
    type = "Gene"
  }
  else {
    type = "TF"
  }
  types <- c(types, type)
}
V(g)$type <- types

# Save the network as GML
file_name = paste(outdir, outname, "tf.gene.ppi_network.gml", sep="")
write_graph(g, file=file_name, format="gml")
print("Network creation succesfull.")