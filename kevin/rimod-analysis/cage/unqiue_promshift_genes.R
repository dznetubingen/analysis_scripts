setwd("/Users/kevin/dzne/")

mapt.ori <- read.table("rimod_package/cage_analysis/promotor_shifting/frontal/mapt_promotor_shifting_genes_fro.txt", stringsAsFactors = F)$V1
grn.ori <- read.table("rimod_package/cage_analysis/promotor_shifting/frontal/grn_promotor_shifting_genes_fro.txt", stringsAsFactors = F)$V1
c9.ori<- read.table("rimod_package/cage_analysis/promotor_shifting/frontal/c9_promotor_shifting_genes_fro.txt", stringsAsFactors = F)$V1


# unique mapt
mapt <- mapt.ori[!mapt.ori %in% grn.ori]
mapt <- mapt[!mapt %in% c9.ori]

# unique grn
grn <- grn.ori[!grn.ori %in% mapt.ori]
grn <- grn[!grn %in% c9.ori]

# unique c9
c9 <- c9.ori[!c9.ori %in% mapt.ori]
c9 <- c9[!c9 %in% grn.ori]

write.table(mapt, "mapt_promshift_unique.txt", quote=F, col.names = F, row.names = F)
write.table(grn, "grn_promshift_unique.txt", quote=F, col.names = F, row.names = F)
write.table(c9, "c9_promshift_unique.txt", quote=F, col.names = F, row.names = F)
