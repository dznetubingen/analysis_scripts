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
                          filename = "cageseq_overlap_venn.png",
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
df$DC <- factor(group, levels = c("FTD-MAPT", "FTD-GRN", "FTD-C9orf72"))

p <- ggplot(df, aes(x=Method, y=Size, fill=DC)) +
  geom_bar(stat = "identity", position="dodge2") +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(. ~ DC, scales = "free") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(size = 12, angle = 30))
p

png("cage_rnaseq_barplot.png", height=300, width=500)
p
dev.off()


#==============================================================#

####
# String analysis plotting
####
setwd("~/rimod/paper/figures/figure2/")
cutoff <- 12
# MAPT
rea.mapt <- read.table("~/rimod/RNAseq/analysis/pathways_analysis/stringdb/mapt_REACT.tsv", sep="\t", header=F, stringsAsFactors = F)
rea.mapt <- rea.mapt[, 1:6]
colnames(rea.mapt) <- c("react", "name", "es", "direction", "size", "fdr")
# change second pathway name
rea.mapt$name[2] <- rea.mapt$react[2]
rea.mapt <- rea.mapt[1:cutoff,]

rea.mapt$negLog <- -log10(rea.mapt$fdr)
rea.mapt <- rea.mapt[order(rea.mapt$negLog),]
rea.mapt$name <- factor(rea.mapt$name, levels = rea.mapt$name)
# make barplot
p <- ggplot(rea.mapt, aes(x=name, y=negLog)) + 
  geom_bar(stat = "identity", fill = mypal[1]) + 
  theme_minimal(base_size = 15) +
  coord_flip() + 
  scale_fill_continuous()
p

png("rnaseq_mapt_string_reactome_barplot.png", height=400, width=400)
p
dev.off()


# GRN
rea.mapt <- read.table("~/rimod/RNAseq/analysis/pathways_analysis/stringdb/grn.enrichment.RCTM.tsv", sep="\t", header=F, stringsAsFactors = F)
rea.mapt <- rea.mapt[, 1:6]
colnames(rea.mapt) <- c("react", "name", "es", "direction", "size", "fdr")
# change second pathway name
rea.mapt <- rea.mapt[1:cutoff,]
rea.mapt$name[5] <- rea.mapt$react[5]
rea.mapt$name[4] <- rea.mapt$react[4]


rea.mapt$negLog <- -log10(rea.mapt$fdr)
rea.mapt <- rea.mapt[order(rea.mapt$negLog),]
rea.mapt$name <- factor(rea.mapt$name, levels = rea.mapt$name)
# make barplot
p <- ggplot(rea.mapt, aes(x=name, y=negLog)) + 
  geom_bar(stat = "identity", fill = mypal[2]) + 
  theme_minimal(base_size = 15) +
  coord_flip() + 
  scale_fill_continuous()
p

png("rnaseq_grn_string_reactome_barplot.png", height=400, width=400)
p
dev.off()

#============================================#

####
# Cell composition plot
####

fracs <- read.table("~/rimod/RNAseq/analysis/deconvolution/cdn_predictions.txt", sep="\t", header=T)
colnames(fracs)[1] <- "sample"
fracs$sample <- gsub("X", "", fracs$sample)
fracs$sample <- str_split(fracs$sample, pattern="_", simplify = T)[,1]
fracs$Neurons <- fracs$InNeurons + fracs$ExNeurons
fracs <- fracs[, c(-3, -9)]



# get design matrix
md <- read.csv("~/rimod/files/FTD_Brain.csv")
md <- md[md$REGION == "frontal",]
md$sample <- str_split(md$GIVENSAMPLENAME, pattern="_", simplify = T)[,1]
md <- md[md$sample %in% fracs$sample,]
md <- md[match(fracs$sample, md$sample),]
fracs$group <- md$DISEASE.CODE

# remove sporadic
fracs <- fracs[!fracs$group == "Sporadic-TDP",]


# melting
fracs <- melt(fracs)

p <- ggplot(fracs, aes(x=sample, y=value, fill=variable)) + 
  geom_bar(stat = "identity", colour="white") + 
  facet_grid(cols = vars(group), scales = "free_x", space="free_x") + 
  theme_minimal(base_size = 15) + 
  theme(axis.text.x = element_blank())
p

png("deconvolution_plot_rnaseq.png", height=300, width=500)
p
dev.off()
