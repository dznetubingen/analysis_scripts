#################
# RiMod miRNA experiment mimic/inhibitor
# MICROGLIA
##############
library(ggplot2)
library(DESeq2)

setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/mirna_mimics/microglia_080920/")

cts <- read.csv("results/salmon/salmon_merged_gene_counts.csv", check.names = F, row.names = 1)

# make the metadata file
groups <- colnames(cts)
groups[grepl("MIR-150-5P-inhi", groups)] <- "Inhibitor150" 
groups[grepl("MIR-150-5P-mimic", groups)] <- "Mimic150"
groups[grepl("MIR-19B-3P-mimic", groups)] <- "Mimic19b"
groups[grepl("MIR-19B-3P-inhi", groups)] <- "Inhibitor19b"
groups[grepl("193a-3PID-inhi", groups)] <- "Inhibitor193a"
groups[grepl("193a-3PID-mimic", groups)] <- "Mimic193a"
groups[grepl("Let7-C-inhi", groups)] <- "InhibitorLet7"
groups[grepl("LetmiR-control-mimic", groups)] <- "MimicLetmir"
groups[grepl("Negcontrol-A-inhi", groups)] <- "NegInhibitor"
groups[grepl("Negcontrol-A-mimic", groups)] <- "NegMimic"

md <- data.frame(sample = colnames(cts), group = groups)
write.table(md, "microglia_mirna_mimic_metadata.txt", sep="\t", quote=F)


cts <- round(cts)
dds <- DESeqDataSetFromMatrix(cts,
                              colData = md,
                              design = ~ group)

dds <- estimateSizeFactors(dds)
#keep <- rowSums(counts(dds, normalized=TRUE) >= 10) >= 3
keep <- rowSums(counts(dds, normalized=TRUE)) > 10
dds <- dds[keep,]


# run DESeq2
dds <- DESeq(dds)


#== Extract results ==#
# 19B Mimic
res.19b <- results(dds, c("group", "Mimic19b", "NegMimic"))
res.19b <- na.omit(res.19b)
deg.19b <- res.19b[res.19b$padj <= 0.05,]

deg.19b.up <- deg.19b[deg.19b$log2FoldChange > 0,]
deg.19b.down <- deg.19b[deg.19b$log2FoldChange < 0,]
write.table(rownames(deg.19b.down), "mimic_19B_down.txt", row.names = F, quote=F)
write.table(rownames(deg.19b.up), "mimic_19B_up.txt", row.names = F, quote=F)

write.table(res.19b, "mimic_micro_19B_res.txt", sep="\t", quote=F)

# 19B Inhibitor
res.inhib.19b <- results(dds, c("group", "Inhibitor19b", "NegInhibitor"))
res.inhib.19b <- na.omit(res.inhib.19b)
deg.inhib.19b <- res.inhib.19b[res.inhib.19b$padj <= 0.05,]

deg.inhib.19b.up <- deg.inhib.19b[deg.inhib.19b$log2FoldChange > 0,]
deg.inhib.19b.down <- deg.inhib.19b[deg.inhib.19b$log2FoldChange < 0,]
write.table(rownames(deg.inhib.19b.down), "inhib_19B_down.txt", row.names = F, quote=F)
write.table(rownames(deg.inhib.19b.up), "inhib_19B_up.txt", row.names = F, quote=F)

write.table(res.inhib.19b, "inhib_micro_19B_res.txt", sep="\t", quote=F)

# 150 Mimic
res.150 <- results(dds, c("group", "Mimic150", "NegMimic"))
res.150 <- na.omit(res.150)
deg.150 <- res.150[res.150$padj <= 0.05,]

deg.150.up <- deg.150[deg.150$log2FoldChange > 0,]
deg.150.down <- deg.150[deg.150$log2FoldChange < 0,]
write.table(rownames(deg.150.down), "mimic_150_down.txt", row.names = F, quote=F)
write.table(rownames(deg.150.up), "mimic_150_up.txt", row.names = F, quote=F)
write.table(res.150, "mimic_micro_150_res.txt", sep="\t", quote=F)


# 150 Inhibitor
res.inhib.150 <- results(dds, c("group", "Inhibitor150", "NegInhibitor"))
res.inhib.150 <- na.omit(res.inhib.150)
deg.inhib.150 <- res.inhib.150[res.inhib.150$padj <= 0.05,]

deg.inhib.150.up <- deg.inhib.150[deg.inhib.150$log2FoldChange > 0,]
deg.inhib.150.down <- deg.inhib.150[deg.inhib.150$log2FoldChange < 0,]
write.table(rownames(deg.inhib.150.down), "inhib_150_down.txt", row.names = F, quote=F)
write.table(rownames(deg.inhib.150.up), "inhib_150_up.txt", row.names = F, quote=F)

write.table(res.inhib.150, "inhib_micro_150_res.txt", sep="\t", quote=F)



# 193a Mimic
res.193a <- results(dds, c("group", "Mimic193a", "MimicLetmir"))
res.193a <- na.omit(res.193a)
deg.193a <- res.193a[res.193a$padj <= 0.05,]

deg.193a.up <- deg.193a[deg.193a$log2FoldChange > 0,]
deg.193a.down <- deg.193a[deg.193a$log2FoldChange < 0,]
write.table(rownames(deg.193a.down), "mimic_193a_down.txt", row.names = F, quote=F)
write.table(rownames(deg.193a.up), "mimic_193a_up.txt", row.names = F, quote=F)

write.table(res.193a, "mimic_micro_193a_res.txt", sep="\t", quote=F)


# 193a Inhibitor
res.inhib.193a <- results(dds, c("group", "Inhibitor193a", "InhibitorLet7"))
res.inhib.193a <- na.omit(res.inhib.193a)
deg.inhib.193a <- res.inhib.193a[res.inhib.193a$padj <= 0.05,]

deg.inhib.193a.up <- deg.inhib.193a[deg.inhib.193a$log2FoldChange > 0,]
deg.inhib.193a.down <- deg.inhib.193a[deg.inhib.193a$log2FoldChange < 0,]
write.table(rownames(deg.inhib.193a.down), "inhib_193a_down.txt", row.names = F, quote=F)
write.table(rownames(deg.inhib.193a.up), "inhib_193a_up.txt", row.names = F, quote=F)

write.table(res.inhib.193a, "inhib_micro_193a_res.txt", sep="\t", quote=F)


# Count transformation and PCA

vst.mat <- varianceStabilizingTransformation(dds)

plotPCA(vst.mat, intgroup="group", ntop=2000)



######
# Make enrichment plots
# (depends on the manually download of gProfiler stuff)
###
library(ggplot2)

make_enrichment_plot <- function(res, cut=20, source="GO:BP", title=""){
  res <- res[, 1:8]
  res <- res[res$source == source,]
  if (nrow(res) >= cut){
    res <- res[1:cut,]  
  }
  res <- res[order(res$adjusted_p_value, decreasing = T),]
  res$term_name <- factor(as.character(res$term_name), levels=as.character(res$term_name))
  
  p <- ggplot(data=res, aes(x=term_name, y=negative_log10_of_adjusted_p_value, color=adjusted_p_value)) + 
    geom_point(aes(size=term_size)) +
    theme_minimal() +
    coord_flip() +
    scale_color_gradient(low="red", high="blue") +
    theme(axis.text.y = element_text(size=14), axis.text.x = element_text(size = 14)) +
    labs(title = title, x="", y="-log10 adj. P-value")
  return(p)
}

setwd("gprofiler_enrichments/")
cutoff <- 10
width = 10 
height = 4
###
# Mir-150
# Mir 150 mimic up
df <- read.csv("mir150_mimic_up.csv")
p <- make_enrichment_plot(df, cut=cutoff, title="miR-150 Mimic Up")
ggsave("mir150_mimic_up.png", width=width, height=height)

  # mir-150 mimic down
  df <- read.csv("mir150_mimic_down.csv")
  p <- make_enrichment_plot(df, cut=cutoff, title="miR-150 Mimic Down")
  ggsave("mir150_mimic_down.png", width=width, height=height)


# Mir 150 inhib up
df <- read.csv("mir150_inhib_up.csv")
p <- make_enrichment_plot(df, cut=cutoff, title="miR-150 inhib Up")
ggsave("mir150_inhib_up.png", width=width, height=height)

# mir-150 inhib down
df <- read.csv("mir150_inhib_down.csv")
p <- make_enrichment_plot(df, cut=cutoff, title="miR-150 inhib Down")
ggsave("mir150_inhib_down.png", width=width, height=height)

#===== end miR-150 ========#


####
# miR-19b
####

# mimic up
df <- read.csv("mir19b_mimic_up.csv")
p <- make_enrichment_plot(df, cut=cutoff, title="miR-19b Mimic up")
ggsave("miR19b_mimic_up.png", width=width, height=height)

# mimic down
df <- read.csv("mir19b_mimic_down.csv")
p <- make_enrichment_plot(df, cut=cutoff, title="miR-19b Mimic down")
ggsave("miR19b_mimic_down.png", width=width, height=height)

# inhib up
df <- read.csv("mir19b_inhib_up.csv")
p <- make_enrichment_plot(df, cut=cutoff, title="miR-19b inhib up")
ggsave("miR19b_inhib_up.png", width=width, height=height)

# inhib down
df <- read.csv("mir19b_inhib_down.csv")
p <- make_enrichment_plot(df, cut=cutoff, title="miR-19b inhib down")
ggsave("miR19b_inhib_down.png", width=width, height=height)

#======= end miR-19b =======#


####
# miR-193a
####

# mimic up
df <- read.csv("mir193a_mimic_up.csv")
p <- make_enrichment_plot(df, cut=cutoff, title="miR-193a mimic up")
ggsave("miR193a_mimic_up.png", width=width, height=height)

# mimic down
df <- read.csv("mir193a_mimic_down.csv")
p <- make_enrichment_plot(df, cut=cutoff, title="miR-193a mimic down")
ggsave("miR193a_mimic_down.png", width=width, height=height)



##########
# Target comparison
#########
targets <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/MAPT_miRNA_target_edge_table.txt", sep="\t", header=T)


# miR-150
mir150 <- targets[grepl("hsa-miR-150-", targets$mirna),]
mir150_down <- rownames(deg.150.down)

length(intersect(mir150$targets, mir150_down))
# miR-193
mir193 <- targets[grepl("hsa-miR-193a", targets$mirna),]
mir193_down <- rownames(deg.193a.down)

length(intersect(mir193$targets, mir193_down))

# mir-19b
targets <- as.character(read.table("~/rimod/smallRNA/frontal/target_selection/targets_miR-19b-3p.txt", sep="\t")$V1)
mir19_down <- rownames(deg.19b.down)
length(intersect(mir19_down, targets))
