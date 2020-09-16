###########
# Venn Diagram for methylation data
##########
library(VennDiagram)
setwd("~/rimod/paper_v2/figures/figure2/")


mapt <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_mapt.ndc_quant.txt", sep="\t", header=T)
grn <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_grn.ndc_quant.txt", sep="\t", header=T)
c9 <- read.table("~/rimod/Methylation/frontal_methylation_0818/DMPs_c9orf72.ndc_quant.txt", sep="\t", header=T)

mapt <- mapt[mapt$adj.P.Val <= 0.05,]
grn <- grn[grn$adj.P.Val <= 0.05,]
c9 <- c9[c9$adj.P.Val <= 0.05,]


dmp.list <- list("FTD-MAPT"=rownames(mapt),
                 "FTD-GRN"=rownames(grn),
                 "FTD-C9orf72"=rownames(c9))

# color palette
mypal <- brewer.pal(3, "Dark2")

#  make VENN plot
venn.plot <- venn.diagram(dmp.list, 
                          filename = "DMPs_overlap_venn.png",
                          cex = 1,
                          
                          # image
                          imagetype = "png",
                          resolution = 300,
                          height=800,
                          width=800,
                          
                          # circles
                          lty = 'blank',
                          fill = mypal,
                          
                          # names
                          cat.cex = 1,
                          cat.pos = c(-27, 27, 135),
                          cat.dist = c(0.055, 0.055, 0.055))