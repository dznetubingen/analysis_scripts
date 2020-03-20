#####################
# RNA-seq GSVA analysis for Angiogenesis verification
#####################
library(stringr)
library(biomaRt)
library(GSVA)
library(pheatmap)
library(viridis)
library(limma)
setwd("~/rimod/RNAseq/analysis/gsva_analysis_endothelial/")

# Load expression data
emat <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_vst_values_2019-10-23_13.33.11.txt", sep="\t", header=T, row.names=1)
rownames(emat) <- str_split(rownames(emat), pattern="[.]", simplify = T)[,1]
colnames(emat) <- gsub("X", "", colnames(emat))

# Load metadata
md <- read.table("~/rimod/RNAseq/rnaseq_frontal_md.txt", sep="\t", header=T)
md <- md[match(colnames(emat), md$ids),]

# Get gene symbols
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=rownames(emat), mart=ensembl)
emat <- merge(emat, bm, by.x="row.names", by.y="ensembl_gene_id")
emat <- emat[!duplicated(emat$hgnc_symbol),]
emat <- na.omit(emat)
emat <- emat[!emat$hgnc_symbol == "",]
rownames(emat) <- emat$hgnc_symbol
emat <- emat[, c(-1, -ncol(emat))]

# Load genesets
hallmark <- getGmt("~/resources/genesets/Human_Reactome_October_01_2019_symbol.gmt", geneIdType = SymbolIdentifier(), collectionType = BroadCollection(category = "h"))

# Run GSVA
gsva.res <- gsva(as.matrix(emat), hallmark)
df <- as.data.frame(gsva.res)

## Make heatmap
mutation <- as.character(md$mutated_gene)
mutation[md$mutation == "P301L"] <- "MAPT_P301L"
colanno <- data.frame(mutation)
rownames(colanno) <- colnames(df)


pheatmap(df, color=viridis(200, option="A"), annotation_col = colanno, show_rownames = F)


## Limma analysis
G <- factor(mutation)
dm <- model.matrix(~ -1 + G)
colnames(dm) <- c("C9orf72", "control", "GRN", "MAPT", "MAPT_P301L")
cont <- c("C9orf72-control", "GRN-control", "MAPT-control", "MAPT_P301L-control", "(MAPT+MAPT_P301L)-control")
fit <- lmFit(df, design=dm)
cm <- makeContrasts(contrasts = cont, levels=dm)
cont.fit <- contrasts.fit(fit, contrasts = cm)
fit = eBayes(cont.fit)
res <- decideTests(fit)
summary(res)

test <- topTable(fit, coef=4, number = Inf, p.value = 0.05)
test.upu <- test[test$logFC > 0,]
rownames(test) <- str_split(rownames(test), pattern="%", simplify = T)[,1]
