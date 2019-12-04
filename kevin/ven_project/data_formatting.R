##########
# Format sc/sn RNA-seq datasets to be easily parseable with Scanpy
# ROSMAP and Lake et al. data is already in good format
# Only have to get the count table for L5 VENs
#########
library(Matrix)

setwd("~/ven_project/data/")

load("L5_VEN/data/FI_layer5_count_data.rda")

write.table(FI_layer5_count_data, "L5_VEN/l5_ven_count_table.txt", sep="\t", quote=F, col.names = NA)
