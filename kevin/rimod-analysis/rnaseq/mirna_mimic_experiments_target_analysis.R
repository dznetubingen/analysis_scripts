###############
# Evaluation of miRNA-target overlap between predicted and 'actual' targets
# according to the experiments
###############

setwd("/media/kevin/89a56127-927e-42c0-80de-e8a834dc81e8/mirna_mimics/")

targets <- read.table("~/rimod/smallRNA/frontal/analysis/target_mrna_correlation_analysis_0819/MAPT_miRNA_target_edge_table.txt", sep="\t", header=T, stringsAsFactors = F)

###
# hsa-miR-150-5p
###
mir150 <- targets[targets$mirna == "hsa-miR-150-5p",]

neur150 <- read.table("mimic_neurons_150_res.txt")
mic150 <- read.table("microglia_080920/mimic_micro_150_res.txt")

neur150 <- na.omit(neur150[mir150$targets,])
mic150 <- na.omit(mic150[mir150$targets,])


neur.mimic.down <- read.table("mimic_150_down.txt")
neur_inter <- intersect(mir150$targets, neur.mimic.down$V1)

mic.mimic.down <- read.table("microglia_080920/mimic_150_down.txt")
mic_inter <- intersect(mir150$targets, mic.mimic.down$V1)

mic.inhib.up <- read.table("microglia_080920/inhib_150_up.txt")
mic_inhib_inter <- intersect(mic.inhib.up$V1, mir150$targets)
#==== end hsa-miR-150-5p ====#



###
# hsa-miR-19b-3p
####
mir19 <- read.table("~/rimod/smallRNA/frontal/target_selection/targets_miR-19b-3p.txt")
neur19 <- read.table("mimic_neurons_19B_res.txt")
mic19 <- read.table("microglia_080920/mimic_micro_19B_res.txt")

neur.mimic.down <- read.table("mimic_19B_down.txt")
mic.mimic.down <- read.table("microglia_080920/mimic_19B_down.txt")

neur.inter <- intersect(neur.mimic.down$V1, mir19$V1)
mic.inter <- intersect(mic.mimic.down$V1, mir19$V1)

mir19.inhib <- read.table("inhib_19B_up.txt")

neur19 <- na.omit(neur19[mir19$V1,])
mic19 <- na.omit(mic19[mir19$V1,])
#==== end hsa-miR-19b-3p ===#

###
# hsa-miR-193a-3p
####
mir193 <- targets[targets$mirna == "hsa-miR-193a-3p",]
mic193 <- read.table("microglia_080920/mimic_micro_193a_res.txt")

mic193 <- na.omit(mic193[mir193$targets,])

mic.mimic.down <- read.table("microglia_080920/mimic_193a_down.txt")
mic.inhib.up <- read.table("microglia_080920/inhib_193a_up.txt")

mic.inter <- intersect(mic.mimic.down$V1, mir193$targets)
mic.inhib.inter <- intersect(mic.inhib.up$V1, mir193$targets)
#=== end hsa-miR-193a-3p ===#