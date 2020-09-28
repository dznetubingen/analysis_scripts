setwd("~/rimod/CAGE/cage_analysis/tf_enrichment_analysis_050420/")

library(biomaRt)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

grn <- read.table("GRN_up_targets.txt")
nfkb <- grn[grn$TF == "NFKB2",]
nfkb_targets <- nfkb$TargetGene

ebf3 <- grn[grn$TF == "EBF3",]
ebf3_targets <- ebf3$TargetGene

# save TF targets
sp1 <- grn[grn$TF == "SP1",]$TargetGene
rela <- grn[grn$TF == "RELA",]$TargetGene
cebpd <- grn[grn$TF == "CEBPD",]$TargetGene
klf3 <- grn[grn$TF == "KLF3",]$TargetGene
write.table(sp1, "GRN_SP1_targets.txt", quote=F, col.names = F, row.names = F)
write.table(rela, "GRN_RELA_targets.txt", quote=F, col.names = F, row.names = F)
write.table(cebpd, "GRN_CEBPD_targets.txt", quote=F, col.names = F, row.names = F)
write.table(klf3, "GRN_KLF3_targets.txt", quote=F, col.names = F, row.names = F)




bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=nfkb_targets, mart=ensembl)
nfkb_targets_hgnc <- bm$hgnc_symbol
nfkb_targets_hgnc <- nfkb_targets_hgnc[!nfkb_targets_hgnc == ""]
write.table(nfkb_targets_hgnc, "GRN_NFKB2_targets.txt", quote=F, col.names = F, row.names = F)

bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=ebf3_targets, mart=ensembl)
ebf3_targets <- bm$hgnc_symbol
ebf3_targets <- ebf3_targets[!ebf3_targets == ""]
write.table(ebf3_targets, "GRN_EBF3_targets.txt", quote=F, col.names = F, row.names = F)


rfx2 <- read.table("MAPT_up_targets.txt")
rfx2 <- rfx2[rfx2$TF == "RFX2",]
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=rfx2$TargetGene, mart=ensembl)
rfx2 <- bm$hgnc_symbol
rfx2 <- rfx2[!rfx2 == ""]
write.table(rfx2, "MAPT_RFX2_targets.txt", quote=F, col.names = F, row.names=F)

# tead
mapt <- read.table("MAPT_up_targets.txt")
grn <- read.table("GRN_up_targets.txt")
mapt <- mapt[mapt$TF == "TEAD2",]
grn <- grn[grn$TF == "TEAD2",]
tead <- intersect(mapt$TargetGene, grn$TargetGene)
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=tead, mart=ensembl)
tead <- bm$hgnc_symbol
tead <- tead[!tead == ""]

write.table(tead, "TEAD2_MAPT_GRN_target_intersection.txt", quote=F, col.names = F, row.names = F)


###
# Down targets
####
tf_name = "EGR3"
tf <- read.table("GRN_down_targets.txt")
tf <- tf[tf$TF == tf_name,]
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=tf$TargetGene, mart=ensembl)
tf <- bm$hgnc_symbol
tf <- tf[!tf == ""]
write.table(tf, "GRN_EGR3_targets.txt", quote=F, col.names = F, row.names = F)

