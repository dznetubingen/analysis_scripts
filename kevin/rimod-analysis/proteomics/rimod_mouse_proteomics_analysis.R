###########
# RiMod Mouse Proteomics analysis
# extract results from collaborators
##########

setwd("~/rimod/Proteomics/mouse/")

# load data
keep_cols <- c(1,4,8,9,10,11,12,13)
mapt <- read.csv("P301L mice_proteomics.csv", sep=";")
mapt <- mapt[, keep_cols]
grn <- read.csv("GRN_KO_proteomics.csv", sep=";")
grn <- grn[, keep_cols]

dim(mapt)

# Make numbers readible by R
mapt$p.Tg.Con.1.5M <- as.numeric(gsub(",", ".", mapt$p.Tg.Con.1.5M))
mapt$p.Tg.Con.3M <- as.numeric(gsub(",", ".", mapt$p.Tg.Con.3M))
mapt$p.Tg.Con.6M <- as.numeric(gsub(",", ".", mapt$p.Tg.Con.6M))

grn$p.old <- as.numeric(gsub(",", ".", grn$p.old))
grn$p.semi <- as.numeric(gsub(",", ".", grn$p.semi))
grn$p.young <- as.numeric(gsub(",", ".", grn$p.young))

# Extract mapt sig genes
mapt.young <- na.omit(mapt[mapt$p.Tg.Con.1.5M < 0.05,])
mapt.middle <- na.omit(mapt[mapt$p.Tg.Con.3M < 0.05,])
mapt.old <- na.omit(mapt[mapt$p.Tg.Con.6M < 0.05,])

write.table(mapt.young$Gene.names, "MAPT_young_DEproteins.txt", quote=F, col.names = F, row.names = F)
write.table(mapt.middle$Gene.names, "MAPT_middle_DEproteins.txt", quote=F, col.names = F, row.names = F)
write.table(mapt.old$Gene.names, "MAPT_old_DEproteins.txt", quote=F, col.names = F, row.names = F)


grn.young <- na.omit(grn[grn$p.young < 0.05,])
grn.middle <- na.omit(grn[grn$p.semi < 0.05,])
grn.old <- na.omit(grn[grn$p.old < 0.05,])
write.table(grn.young$Gene.names, "GRN_young_DEproteins.txt", quote=F, col.names = F, row.names = F)
write.table(grn.middle$Gene.names, "GRN_middle_DEproteins.txt", quote=F, col.names = F, row.names = F)
write.table(grn.old$Gene.names, "GRN_old_DEproteins.txt", quote=F, col.names = F, row.names = F)
