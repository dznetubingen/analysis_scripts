###############################################
## Create TF-target network
##
## This script creates a mapping from TF-gene 
## from the Results created by the tf-activity pipeline
###############################################

# load libs
library(biomaRt)
library(plyr)
library(igraph)


args <- commandArgs(trailingOnly = TRUE)

tf_targets_file = args[1]
outname = args[2]

targets <- read.table(tf_targets_file, sep="\t", header = T, stringsAsFactors = F)

# only include promotor-TSS targets
dim(targets)
targets <- targets[grepl("promoter-TSS", targets$Annotation),]

## Extract names of the transcripts
extractTranscriptNames = function(transc){
  newt <- as.character(sapply(transc, function(x){gsub("[(]", "", gsub("[)]", "", strsplit(x, split=" ")[[1]][[2]]))}))
  newt <- as.character(sapply(newt, function(x){strsplit(x, split="[.]")[[1]][[1]]}))
  return(newt)
}
transcripts <- extractTranscriptNames(targets$Annotation)
targets$ensembl_transcript_id <- transcripts
                                   
# Change Transcripts to Genes with biomaRt
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm <- getBM(values = transcripts, attributes = c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol"), filters = "ensembl_transcript_id", mart = ensembl)
bm <- bm[!duplicated(bm),]
mdata <- join(targets, bm, by=c("ensembl_transcript_id"), type="left", match="all")

df = mdata[,c(3,4,5,7,8,9)]
cols <- colnames(df)
cols[cols == "X0"] <- "no_hits"
colnames(df) <- cols

# Remove genes with no ensembl_gene_Id
df <- df[!is.na(df$ensembl_gene_id),]

### Compare with ChIP-seq results
if (length(args) == 3){
  encode_ovl_file = args[3]
  chip <- read.table(encode_ovl_file, sep="\t", header = T)
  chip <- chip[chip$Overlap >= 10,]
  df <- df[df$TF %in% chip$TF,]
}


#######
## Create igraph object
#######

edges <- c()
for (i in 1:nrow(df)){
  e <- c(df[i,3], df[i,5])
  edges <- c(edges, e)
}
g <- graph(edges=edges, directed=T)
l <- layout_with_fr(g)
plot(g, vertex.color = "gold", vertex.size = 5, edge.arrow.size = .5, layout = l, 
     vertex.label.cex = 0.6, vertex.label.color = "black")

### Assign labels

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
E(g)$hits <- df$no_hits

write.table(df, paste(outname, "_tf_target_mapping.txt", sep=""), sep="\t", quote=F, row.names=F)
write_graph(g, file =paste(outname, "_graph_cyto.gml", sep=""), format="gml")
