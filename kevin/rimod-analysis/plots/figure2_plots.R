#######################
# Plots for Figure 2 of RiMod paper
#######################
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
library(reshape2)
library(stringr)
setwd("~/rimod/paper/figures/figure2/")

#==============================================================#

####
# String analysis plotting
####
pval_line <- -log10(0.05)
setwd("~/rimod/paper/figures/figure2/")
cutoff <- 12
# MAPT
rea.mapt <- read.table("~/rimod/RNAseq/analysis/pathways_analysis/stringdb/cc_filtered_string/mapt.enrichment.RCTM.tsv", sep="\t", header=F, stringsAsFactors = F)
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
  geom_bar(stat = "identity", fill = mypal[1]) + 
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
rea.mapt <- read.table("~/rimod/RNAseq/analysis/pathways_analysis/stringdb/cc_filtered_string/grn_ccF.enrichment.RCTM.tsv", sep="\t", header=F, stringsAsFactors = F)
rea.mapt <- rea.mapt[, 1:6]
colnames(rea.mapt) <- c("react", "name", "es", "direction", "size", "fdr")
# change second pathway name
rea.mapt <- rea.mapt[1:cutoff,]
rea.mapt$name[3] <- rea.mapt$react[3]
#rea.mapt$name[4] <- rea.mapt$react[4]
rea.mapt$name[10] <- rea.mapt$react[10]
rea.mapt$name[12] <- rea.mapt$react[12]


rea.mapt$negLog <- -log10(rea.mapt$fdr)
rea.mapt <- rea.mapt[order(rea.mapt$negLog),]
rea.mapt$name <- factor(rea.mapt$name, levels = rea.mapt$name)
# make barplot
p <- ggplot(rea.mapt, aes(x=name, y=negLog)) + 
  geom_bar(stat = "identity", fill = mypal[2]) + 
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
mir.reac <- read.csv("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/reactome_pathway_enrichment/mapt_downTargets.csv")
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
mir.reac <- read.csv("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/reactome_pathway_enrichment/grn_downTargets.csv")
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
mir.reac <- read.csv("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/reactome_pathway_enrichment/c9orf72_downTargets.csv")
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

fracs <- read.table("~/rimod/RNAseq/analysis/deconvolution/cdn_predictions.txt", sep="\t", header=T)
colnames(fracs)[1] <- "sample"
fracs$sample <- gsub("X", "", fracs$sample)
fracs$sample <- str_split(fracs$sample, pattern="_", simplify = T)[,1]
fracs$Neurons <- fracs$InNeurons + fracs$ExNeurons
fracs <- fracs[, c(-3, -9)]



# get design matrix
#md <- read.csv("/Users/kevin/dzne/rimod_package/files/FTD_Brain.csv")
md <- read.csv("~/rimod/files/FTD_Brain_corrected.csv")
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
as.mapt <- read.table("~/rimod/RNAseq/as_analysis/majiq/mapt_AS_genes_dPSI_0.2.txt", sep="\t", header=T, stringsAsFactors = F)$x
as.grn <- read.table("~/rimod/RNAseq/as_analysis/majiq/grn_AS_genes_dPSI_0.2.txt", sep="\t", header=T, stringsAsFactors = F)$x
as.c9 <- read.table("~/rimod/RNAseq/as_analysis/majiq/c9orf72_AS_genes_dPSI_0.2.txt", sep="\t", header=T, stringsAsFactors = F)$x

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
                          cat.dist = c(0.055, 0.055, 0.055))

test <- intersect(as.mapt, intersect(as.grn, as.c9))
write.table(test, "test.txt", row.names=F, quote=F, col.names=F)
