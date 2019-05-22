############################
# Angiogenesis GSVA plot
############################
library(GSVA)
library(GSEABase)
library(biomaRt)
library(pheatmap)
library(viridis)
library(limma)
library(stringr)

setwd("~/rimod/CAGE/cage_analysis/angiogenesis_gsva/")
region = "fro"

##########
## Data preparation
###########
# Load expression data
emat <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_2018-04-26_genewise/deseq_normalized_counts_2018-04-26_14.22.04.txt", sep = "\t", header = T, row.names=1)

md <- read.csv("~/rimod/files/FTD_Brain.csv", stringsAsFactors = FALSE)
md$SAMPLEID <- str_pad(md$SAMPLEID, width = 5, side = "left", pad = "0") # fill sample ids to 5 digits

# bring counts and md in similar format
cage.samples <- as.character(gsub("sample_","",colnames(emat)))
cage.samples <- as.character(sapply(cage.samples, function(x){strsplit(x, split=paste("_", region, sep=""))[[1]][[1]]}))

md <- md[md$SAMPLEID %in% cage.samples,]
md <- md[match(cage.samples, md$SAMPLEID),]
disease.codes <- c("FTD-C9", "FTD-MAPT", "FTD-GRN", "control")
keep <- md$DISEASE.CODE %in% disease.codes
md <- md[keep,]
emat <- emat[,keep]
md$DISEASE.CODE <- gsub("-", "_", md$DISEASE.CODE) # make disease code names safe
md$DISEASE.CODE[md$GENE == "P301L"] <- "FTD-MAPTP301l"

# Get gene symbols
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
bm <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = rownames(emat), mart = ensembl)

emat <- merge(emat, bm, by.x="row.names", by.y="ensembl_gene_id")
emat <- emat[!duplicated(emat$hgnc_symbol),]
rownames(emat) <- emat$hgnc_symbol
emat <- emat[,c(-1, -ncol(emat))]

# Load genesets
c2 <- getGmt("~/resources/genesets/c2.all.v6.1.symbols.gmt", geneIdType = SymbolIdentifier(), collectionType = BroadCollection(category = "c2"))
hallmark <- getGmt("~/resources/genesets/h.all.v6.1.symbols.gmt", geneIdType = SymbolIdentifier(), collectionType = BroadCollection(category = "h"))


#################
## Run GSVA and plot results
##################

gsva.res <- gsva(as.matrix(emat), hallmark)



# Limma analysis
G <- factor(md$DISEASE.CODE)
dm <- model.matrix(~ -1 + G)
colnames(dm) <- c("control", "C9", "GRN", "MAPT", "MAPTP301L")
cont = c("C9-control", "GRN-control", "MAPT-control", "(MAPT+GRN+C9)-control", "MAPTP301L-control")
fit = lmFit(gsva.res, design = dm)
cm <- makeContrasts(contrasts = cont, levels = dm)
cont.fit <- contrasts.fit(fit, contrasts = cm)
fit <- eBayes(cont.fit)
res <- decideTests(fit)
summary(res)

# Plotting
col_df <- data.frame(group = md$DISEASE.CODE)
rownames(col_df) <- colnames(gsva.res)

res <- topTable(fit, coef=2, number = Inf, p.value = 0.05)
keep <- row.names(res)
pheatmap(gsva.res[keep,], annotation_col = col_df, color = viridis(200, option="A"), show_colnames = F)


## Get genesets
angiogenesis <- hallmark[[40]]
