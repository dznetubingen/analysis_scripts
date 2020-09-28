###########
# Make metadata templates for neuron RNA-seq experiments
###########

setwd("~/rimod/daisy_data/")

####
# iPSC microglia experiments
####

cts <- read.csv("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/mirna_mimics/microglia_080920/results/salmon/salmon_merged_gene_counts.csv", check.names = F, row.names = 1)

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

md <- data.frame(SampleID = colnames(cts), Group = groups, Name = colnames(cts))

write.table(md, "microglia_miRNA_experiments_metadata.txt", sep="\t", quote=F, row.names = F)

#================ end microglia experiments ======================#


####
# iPSC neurons miRNA experiments metadata
####

md <- read.table("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/mirna_mimics/rimod_mirna_neurons_md.txt", sep="\t", header=T)
md$Name <- md$sample
colnames(md) <- c("SampleID", "Group", "Name")

write.table(md, "neurons_miRNA_experiments_metadata.txt", sep="\t", quote=F, row.names=F)

#===================== end neuron miRNA experiments ===========#


####
# iPSC neurons experiment


cts <- read.table("~/rimod/Neurons/rimod_neurons/frontal_lengthScaeldTPM_counts.txt", sep="\t", header=T, row.names = 1, check.names = F)


md <- read.table("~/rimod/Neurons/rimod_neurons/iPSCNeurons_mRNAseq_metadata.txt", sep="\t", header=T, stringsAsFactors = F)
md <- md[md$sample %in% colnames(cts),]
md <- md[match(colnames(cts), md$sample),]
colnames(md) <- c("SampleID", "Mutation", "Group", "Gender", "Flowcell")

md$Name <- md$SampleID

write.table(md, "iPSC_neurons_rimod_metadata.txt", sep="\t", quote=F, row.names = F)
#============= end iPSC neurons ============#