#######################
# Plots for Figure 2 of RiMod paper
#######################
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
setwd("~/rimod/paper/figures/figure2/")

#####
# RNA-seq plots
####

### Load RNA-seq results
rna.mapt <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_result_mapt.ndc_fro_2019-10-23_13.33.11.txt",
                       sep="\t", header=T, row.names = 1)
rna.grn <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_result_grn.ndc_fro_2019-10-23_13.33.11.txt",
                      sep="\t", header=T, row.names=1)
rna.c9 <- read.table("~/rimod/RNAseq/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_result_c9.ndc_fro_2019-10-23_13.33.11.txt",
                     sep="\t", header=T, row.names=1)


# define p-value cutoff
pval_cut <- 0.05

# pval cutting
rna.mapt <- rna.mapt[rna.mapt$padj <= pval_cut,]
rna.grn <- rna.grn[rna.grn$padj <= pval_cut,]
rna.c9 <- rna.c9[rna.c9$padj <= pval_cut,]

rna.list <- list("FTD-MAPT"=rownames(rna.mapt),
                 "FTD-GRN"=rownames(rna.grn),
                 "FTD-C9orf72"=rownames(rna.c9))

mypal <- brewer.pal(3, "Dark2")
# Plot Venn
venn.plot <- venn.diagram(rna.list, 
                          filename = "rnaseq_overlap_venn.png",
                          imagetype = "png",
                          resolution = 300,
                          fill = mypal,
                          cex = 1,
                          cat.cex = 1,
                          height=800,
                          width=800)



#==============================================================================#

#####
# CAGE-seq plots
#####
cage.mapt <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_fro_2019-10-23_13.30.53/deseq_result_mapt.ndc_2019-10-23_13.30.53.txt",
                       sep="\t", header=T, row.names = 1)
cage.grn <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_fro_2019-10-23_13.30.53/deseq_result_grn.ndc_2019-10-23_13.30.53.txt",
                      sep="\t", header=T, row.names=1)
cage.c9 <- read.table("~/rimod/CAGE/cage_analysis/CAGE_deseq_analysis_fro_2019-10-23_13.30.53/deseq_result_c9.ndc_2019-10-23_13.30.53.txt",
                     sep="\t", header=T, row.names=1)


# define p-value cutoff
pval_cut <- 0.05

# pval cutting
cage.mapt <- cage.mapt[cage.mapt$padj <= pval_cut,]
cage.grn <- cage.grn[cage.grn$padj <= pval_cut,]
cage.c9 <- cage.c9[cage.c9$padj <= pval_cut,]

cage.list <- list("FTD-MAPT"=rownames(cage.mapt),
                 "FTD-GRN"=rownames(cage.grn),
                 "FTD-C9orf72"=rownames(cage.c9))

# Plot Venn
venn.plot <- venn.diagram(cage.list, 
                          filename = "cage_overlap_venn.png",
                          imagetype = "png",
                          resolution = 300,
                          fill = mypal,
                          cex = 1,
                          cat.cex = 1,
                          height=800,
                          width=800)






#==============================================================================#


####
# CAGE-seq and RNA-seq overlap
####
# MAPT
mapt.list <- list("CAGE-seq"=rownames(cage.mapt),
                  "RNA-seq"=rownames(rna.mapt))
venn.plot <- venn.diagram(mapt.list, 
                          filename = "MAPT_rna_cage_overlap.png",
                          imagetype = "png",
                          resolution = 300,
                          fill = c(mypal[1], mypal[2]),
                          cex = 2,
                          cat.cex = 2)

# GRN
grn.list <- list("CAGE-seq"=rownames(cage.grn),
                  "RNA-seq"=rownames(rna.grn))
venn.plot <- venn.diagram(grn.list, 
                          filename = "GRN_rna_cage_overlap.png",
                          imagetype = "png",
                          resolution = 300,
                          fill = c(mypal[1], mypal[2]),
                          cex = 2,
                          cat.cex = 2)

# C9orf72
c9.list <- list("CAGE-seq"=rownames(cage.c9),
                 "RNA-seq"=rownames(rna.c9))
venn.plot <- venn.diagram(c9.list, 
                          filename = "C9orf72_rna_cage_overlap.png",
                          imagetype = "png",
                          resolution = 300,
                          fill = c(mypal[1], mypal[2]),
                          cex = 2,
                          cat.cex = 2)

###
# Barplot comparision
###
group <- c("FTD-MAPT","FTD-MAPT","FTD-MAPT", "FTD-GRN","FTD-GRN","FTD-GRN", "FTD-C9orf72", "FTD-C9orf72", "FTD-C9orf72")
name <- c("RNA-seq", "CAGE-seq", "Intersect", "RNA-seq", "CAGE-seq", "Intersect", "RNA-seq", "CAGE-seq", "Intersect")
size <- c(nrow(rna.mapt), nrow(cage.mapt), length(intersect(rownames(rna.mapt), rownames(cage.mapt))),
          nrow(rna.grn), nrow(cage.grn), length(intersect(rownames(rna.grn), rownames(cage.grn))),
          nrow(rna.c9), nrow(cage.c9), length(intersect(rownames(rna.c9), rownames(cage.c9))))

df <- data.frame(DC = group, Method = name, Size = size)

p <- ggplot(df, aes(x=DC, y=Size, fill=Method)) +
  geom_bar(stat = "identity", position="dodge2")
p
ggsave(filename = "cage_rnaseq_barplot.png")







