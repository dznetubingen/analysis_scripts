#######################
# Plots for Figure 2 of RiMod paper
#######################
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
library(reshape2)
library(stringr)
library(extrafont)
loadfonts()
setwd("/Users/kevin/dzne/rimod_analysis/figure2/")

base_path = "/Users/kevin/dzne/"

# Setup color palette
# color palette for NDC, C9orf72, GRN, MAPT and MAPT-P301L
mypal_long <- c()
# Color palette for all groups incl. NDC
mypal <- c()
# color palette only for disease groups
mypal_short <- c("#67e08a", "#db6e1a", "#7570B3")

#####
# RNA-seq plots
####

### Load RNA-seq results
rna.mapt <- read.table("/Users/kevin/dzne/rimod_package/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_result_mapt.ndc_fro_2019-10-23_13.33.11.txt",
                       sep="\t", header=T, row.names = 1)
rna.grn <- read.table("/Users/kevin/dzne/rimod_package/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_result_grn.ndc_fro_2019-10-23_13.33.11.txt",
                      sep="\t", header=T, row.names=1)
rna.c9 <- read.table("/Users/kevin/dzne/rimod_package/analysis/RNAseq_analysis_fro_2019-10-23_13.33.11/deseq_result_c9.ndc_fro_2019-10-23_13.33.11.txt",
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
                          fill = mypal_short,
                          
                          # names
                          cat.cex = 1,
                          cat.pos = c(-27, 27, 135),
                          cat.dist = c(0.055, 0.055, 0.055),
                          
                          # fonts
                          cat.fontfamily = c("Arial", "Arial", "Arial"),
                          fontfamily = c("Arial", "Arial", "Arial", "Arial", "Arial", "Arial", "Arial"))


#####
# smRNA-seq plots
######

srna.mapt <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/analysis_0719/deseq_result_mapt.ndc_frontal_smRNAseq.txt",
                       sep="\t", header=T, row.names = 1)
srna.grn <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/analysis_0719/deseq_result_grn.ndc_frontal_smRNAseq.txt",
                      sep="\t", header=T, row.names=1)
srna.c9 <- read.table("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/analysis_0719/deseq_result_c9.ndc_frontal_smRNAseq.txt",
                     sep="\t", header=T, row.names=1)



# define p-value cutoff
pval_cut <- 0.05

# pval cutting
srna.mapt <- srna.mapt[srna.mapt$padj <= pval_cut,]
srna.grn <- srna.grn[srna.grn$padj <= pval_cut,]
srna.c9 <- srna.c9[srna.c9$padj <= pval_cut,]

srna.list <- list("FTD-MAPT"=rownames(srna.mapt),
                 "FTD-GRN"=rownames(srna.grn),
                 "FTD-C9orf72"=rownames(srna.c9))

mypal <- brewer.pal(3, "Dark2")
# Plot Venn
venn.plot <- venn.diagram(srna.list, 
                          filename = "smrnaseq_overlap_venn.png",
                          cex = 1,
                          
                          # image
                          imagetype = "png",
                          resolution = 300,
                          height=800,
                          width=800,
                          
                          # circles
                          lty = 'blank',
                          fill = mypal_short,
                          
                          # names
                          cat.cex = 1,
                          cat.pos = c(-27, 27, 135),
                          cat.dist = c(0.055, 0.055, 0.055),
                          
                          # fonts
                          cat.fontfamily = c("Arial", "Arial", "Arial"),
                          fontfamily = c("Arial", "Arial", "Arial", "Arial", "Arial", "Arial", "Arial"))

#==============================================================================#



#==============================================================================#

#####
# CAGE-seq plots
#####
cage.mapt <- read.table("/Users/kevin/dzne/rimod_package/cage_analysis/CAGE_deseq_analysis_fro_2019-10-23_13.30.53/deseq_result_mapt.ndc_2019-10-23_13.30.53.txt",
                       sep="\t", header=T, row.names = 1)
cage.grn <- read.table("/Users/kevin/dzne/rimod_package/cage_analysis/CAGE_deseq_analysis_fro_2019-10-23_13.30.53/deseq_result_grn.ndc_2019-10-23_13.30.53.txt",
                      sep="\t", header=T, row.names=1)
cage.c9 <- read.table("/Users/kevin/dzne/rimod_package/cage_analysis/CAGE_deseq_analysis_fro_2019-10-23_13.30.53/deseq_result_c9.ndc_2019-10-23_13.30.53.txt",
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
                          fill = mypal_short,
                          
                          # names
                          cat.cex = 1,
                          cat.pos = c(-27, 27, 135),
                          cat.dist = c(0.055, 0.055, 0.055),
                          
                          # fonts
                          cat.fontfamily = c("Arial", "Arial", "Arial"),
                          fontfamily = c("Arial", "Arial", "Arial", "Arial", "Arial", "Arial", "Arial"))






#==============================================================================#


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
pval_line <- -log10(0.05)
cutoff <- 12
# MAPT
rea.mapt <- read.table("/Users/kevin/dzne/rimod_package/analysis/pathways_analysis/stringdb/mapt_REACT.tsv", sep="\t", header=F, stringsAsFactors = F)
rea.mapt <- rea.mapt[, 1:6]
colnames(rea.mapt) <- c("react", "name", "es", "direction", "size", "fdr")
# change second pathway name
rea.mapt$name[2] <- rea.mapt$react[2]
rea.mapt$name[12] <- rea.mapt$react[12]
rea.mapt <- rea.mapt[1:cutoff,]

rea.mapt$negLog <- -log10(rea.mapt$fdr)
rea.mapt <- rea.mapt[order(rea.mapt$negLog),]
rea.mapt$name <- factor(rea.mapt$name, levels = rea.mapt$name)
# make barplot
p <- ggplot(rea.mapt, aes(x=name, y=negLog)) + 
  geom_bar(stat = "identity", fill = mypal_short[1]) + 
  theme_minimal(base_size = 15) + 
  xlab("") +
  ylim(0, 20) +
  geom_hline(yintercept = pval_line, linetype="dotted") +
  coord_flip() + 
  scale_fill_continuous() +
  theme(axis.text.y = element_text(angle = 30))
p

png("rnaseq_mapt_string_reactome_barplot.png", height=300, width=450)
p
dev.off()


# GRN
rea.mapt <- read.table("/Users/kevin/dzne/rimod_package/analysis/pathways_analysis/stringdb/grn.enrichment.RCTM.tsv", sep="\t", header=F, stringsAsFactors = F)
rea.mapt <- rea.mapt[, 1:6]
colnames(rea.mapt) <- c("react", "name", "es", "direction", "size", "fdr")
# change second pathway name
rea.mapt <- rea.mapt[1:cutoff,]
rea.mapt$name[5] <- rea.mapt$react[5]
rea.mapt$name[4] <- rea.mapt$react[4]
rea.mapt$name[10] <- rea.mapt$react[10]
rea.mapt$name[9] <- rea.mapt$react[9]


rea.mapt$negLog <- -log10(rea.mapt$fdr)
rea.mapt <- rea.mapt[order(rea.mapt$negLog),]
rea.mapt$name <- factor(rea.mapt$name, levels = rea.mapt$name)
# make barplot
p <- ggplot(rea.mapt, aes(x=name, y=negLog)) + 
  geom_bar(stat = "identity", fill = mypal_short[2]) + 
  theme_minimal(base_size = 15) +
  xlab("") +
  geom_hline(yintercept = pval_line, linetype="dotted") +
  ylim(0, 20) +
  coord_flip() + 
  scale_fill_continuous()+
  theme(axis.text.y = element_text(angle = 30))
p

png("rnaseq_grn_string_reactome_barplot.png", height=300, width=450)
p
dev.off()

#============================================#


#####
# miRNA-target pathway analysis
#####

pval_line <- -log10(0.05)

# mapt
mir.reac <- read.csv("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/target_mrna_correlation_analysis_0819/reactome_pathway_enrichment/mapt_downTargets.csv")
mir.reac <- mir.reac[, c(1,2,3,4,5)]
mir.reac <- mir.reac[mir.reac$source == "REAC",]
mir.reac <- mir.reac[1:5,]
mir.reac <- mir.reac[order(mir.reac$negative_log10_of_adjusted_p_value),]
mir.reac$term_id <- factor(mir.reac$term_id, levels = mir.reac$term_id)

p <- ggplot(mir.reac, aes(x=term_id, y=negative_log10_of_adjusted_p_value)) + 
  geom_bar(stat = "identity", fill = mypal[1]) + 
  theme_minimal(base_size = 15) +
  geom_hline(yintercept = pval_line, linetype="dotted") +
  xlab("") +
  ylab("") +
  ylim(0, 9) +
  coord_flip() + 
  scale_fill_continuous()
p

png("mapt_mirTargets_reactome.png", height=200, width=400)
p
dev.off()


# grn
mir.reac <- read.csv("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/target_mrna_correlation_analysis_0819/reactome_pathway_enrichment/grn_downTargets.csv")
mir.reac <- mir.reac[, c(1,2,3,4,5)]
mir.reac <- mir.reac[mir.reac$source == "REAC",]
mir.reac <- mir.reac[order(mir.reac$negative_log10_of_adjusted_p_value),]
mir.reac$term_id <- factor(mir.reac$term_id, levels = mir.reac$term_id)

p <- ggplot(mir.reac, aes(x=term_id, y=negative_log10_of_adjusted_p_value)) + 
  geom_bar(stat = "identity", fill = mypal[2]) + 
  theme_minimal(base_size = 15) +
  xlab("") +
  ylim(0, 9) +
  geom_hline(yintercept = pval_line, linetype="dotted") +
  coord_flip() + 
  scale_fill_continuous()
p
png("grn_mirTargets_reactome.png", height=200, width=400)
p
dev.off()


# c9orf72
mir.reac <- read.csv("/Users/kevin/dzne/rimod_package/smRNAseq/analysis/target_mrna_correlation_analysis_0819/reactome_pathway_enrichment/c9orf72_downTargets.csv")
mir.reac <- mir.reac[, c(1,2,3,4,5)]
mir.reac <- mir.reac[mir.reac$source == "REAC",]
mir.reac <- mir.reac[order(mir.reac$negative_log10_of_adjusted_p_value),]
mir.reac$term_name <- factor(mir.reac$term_name, levels = mir.reac$term_name)

p <- ggplot(mir.reac, aes(x=term_name, y=negative_log10_of_adjusted_p_value)) + 
  geom_bar(stat = "identity", fill = mypal[3]) + 
  theme_minimal(base_size = 15) +
  geom_hline(yintercept = pval_line, linetype="dotted") +
  coord_flip() + 
  scale_fill_continuous()
p

png("c9orf72_mirTargets_reactome.png", height=400, width=400)
p
dev.off()



####
# Cell composition plot
####

fracs <- read.table("/Users/kevin/dzne/rimod_package/analysis/deconvolution/cdn_predictions.txt", sep="\t", header=T)
colnames(fracs)[1] <- "sample"
fracs$sample <- gsub("X", "", fracs$sample)
fracs$sample <- str_split(fracs$sample, pattern="_", simplify = T)[,1]
fracs$Neurons <- fracs$InNeurons + fracs$ExNeurons
fracs <- fracs[, c(-3, -9)]



# get design matrix
md <- read.csv("/Users/kevin/dzne/rimod_package/files/FTD_Brain.csv")
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

# endothelial cells
end <- fracs[fracs$variable == "Endothelial",]
p <- ggplot(end, aes(x=sample, y=value, fill=group)) + 
  geom_bar(stat = "identity", colour="white") + 
  facet_grid(cols = vars(group), scales = "free_x", space="free_x") + 
  theme_minimal(base_size = 15) + 
  theme(axis.text.x = element_blank())
p

p <- ggplot(end, aes(x=group, y=value, fill=group)) + 
  geom_boxplot() + 
  theme_minimal(base_size = 15) + 
  theme(axis.text.x = element_blank())
p

####
# Alternative splicing plots
####
as.mapt <- read.table("/Users/kevin/dzne/rimod_package/as_analysis/mapt_AS_genes_dPSI_0.2.txt", sep="\t", header=T, stringsAsFactors = F)$x
as.grn <- read.table("/Users/kevin/dzne/rimod_package/as_analysis/grn_AS_genes_dPSI_0.2.txt", sep="\t", header=T, stringsAsFactors = F)$x
as.c9 <- read.table("/Users/kevin/dzne/rimod_package/as_analysis/c9orf72_AS_genes_dPSI_0.2.txt", sep="\t", header=T, stringsAsFactors = F)$x

altsplice.list <- list("FTD-MAPT"=as.mapt,
                  "FTD-GRN"=as.grn,
                  "FTD-C9orf72"=as.c9)

# Plot Venn
venn.plot <- venn.diagram(altsplice.list, 
                          filename = "AS_overlap_venn.png",
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
                          cat.dist = c(0.055, 0.055, 0.055),
                          
                          # fonts
                          cat.fontfamily = c("Arial", "Arial", "Arial"),
                          fontfamily = c("Arial", "Arial", "Arial", "Arial", "Arial", "Arial", "Arial"))

test <- intersect(as.mapt, intersect(as.grn, as.c9))
write.table(test, "test.txt", row.names=F, quote=F, col.names=F)
