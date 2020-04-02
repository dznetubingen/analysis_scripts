######
# Human Proteomics analysis
#####
library(limma)
library(RColorBrewer)
library(pheatmap)
library(viridis)

setwd("~/dzne/rimod_package/Proteomics/")
out_dir = "/home/kevin/dzne/rimod_analysis/proteomics_analysis_010420/"

# load data
prot <- read.table("RIMOD_data_2017-04-27/swath_data_2017-04-27.tsv", sep="\t", header=T)
prot.ids <- prot[,1:2]
prot <- prot[,c(-1,-2)]
rownames(prot) <- prot.ids$GENE_SYMBOL

# Convert values to log2 scale
prot <- log2(prot)

# load metadata and order it according to data
design <- read.table("RIMOD_data_2017-04-27/metadata_2017-04-27.tsv", sep="\t", header=T)
design <- design[match(design$Sample.Label, colnames(prot)),]

# subset frontal and temporal
keep_front <- design$Brain.area == "frontal"
keep_temp <- design$Brain.area == "temporal"
front <- prot[, keep_front]
front.design <- design[keep_front,]
temp <- prot[, keep_temp]
temp.design <- design[keep_temp,]


####
# Frontal analysis
####

### Differential expression analysis
gender <- front.design$Gender
ph <- as.numeric(front.design$pH)
dc <- factor(as.character(front.design$Disease.Code))
dm.front <- model.matrix(~ -1 + dc + gender)
dm.front <- dm.front[,-5]
colnames(dm.front) <- c("FTD.C9", "FTD.GRN", "FTD.MAPT", "NDC", "genderM")

# Fit the model
fit <- lmFit(front, design = dm.front)
conts <- c("FTD.C9-NDC", "FTD.GRN-NDC", "FTD.MAPT-NDC")
cm <- makeContrasts(contrasts = conts, levels=dm.front)
cont.fit <- contrasts.fit(fit, contrasts = cm)
fit2 <- eBayes(cont.fit)
res <- decideTests(fit2)
summary(res)

# extract fold-changes
deg.c9 <- topTable(fit2, coef=1, number=Inf, p.value=1)
deg.grn <- topTable(fit2, coef=2, number=Inf, p.value=1)
deg.mapt <- topTable(fit2, coef=3, number=Inf, p.value=1)

deg.c9 <- data.frame(rownames(deg.c9), deg.c9$logFC)
deg.grn <- data.frame(rownames(deg.grn), deg.grn$logFC)
deg.mapt <- data.frame(rownames(deg.mapt), deg.mapt$logFC)

write.table(deg.c9, paste0(out_dir, "frontal_c9_lfcs.txt"), quote=F, row.name=F, col.names = F)
write.table(deg.grn, paste0(out_dir, "frontal_grn_lfcs.txt"), quote=F, row.name=F, col.names = F)
write.table(deg.mapt, paste0(out_dir, "frontal_mapt_lfcs.txt"), quote=F, row.name=F, col.names = F)

#=======================#


###
# Temporal analysis
####
gender <- temp.design$Gender
ph <- as.numeric(temp.design$pH)
dc <- factor(as.character(temp.design$Disease.Code))
dm.temp <- model.matrix(~ -1 + dc + gender)
dm.temp <- dm.temp[,-5]
colnames(dm.temp) <- c("FTD.C9", "FTD.GRN", "FTD.MAPT", "NDC", "genderM")

# Fit the model
fit <- lmFit(temp, design = dm.temp)
conts <- c("FTD.C9-NDC", "FTD.GRN-NDC", "FTD.MAPT-NDC")
cm <- makeContrasts(contrasts = conts, levels=dm.temp)
cont.fit <- contrasts.fit(fit, contrasts = cm)
fit2 <- eBayes(cont.fit)
res <- decideTests(fit2)
summary(res)

# extract fold-changes
deg.c9 <- topTable(fit2, coef=1, number=Inf, p.value=1)
deg.grn <- topTable(fit2, coef=2, number=Inf, p.value=1)
deg.mapt <- topTable(fit2, coef=3, number=Inf, p.value=1)

deg.c9 <- data.frame(rownames(deg.c9), deg.c9$logFC)
deg.grn <- data.frame(rownames(deg.grn), deg.grn$logFC)
deg.mapt <- data.frame(rownames(deg.mapt), deg.mapt$logFC)

write.table(deg.c9, paste0(out_dir, "temporal_c9_lfcs.txt"), quote=F, row.name=F, col.names = F)
write.table(deg.grn, paste0(out_dir, "temporal_grn_lfcs.txt"), quote=F, row.name=F, col.names = F)
write.table(deg.mapt, paste0(out_dir, "temporal_mapt_lfcs.txt"), quote=F, row.name=F, col.names = F)
