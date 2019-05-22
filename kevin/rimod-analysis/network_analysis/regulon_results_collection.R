##############################
## Collection and combination of all regulon analysis results
##############################
setwd("~/rimod/CAGE/cage_analysis/regulon_all_frontal/")


####
# MAPT regulons
####
# Up-regulated
mapt_up_bg <- readRDS("regulon_analysis_250418/mapt/mapt_sig_regulons.rds")
mapt_up_degbg <- readRDS("regulon_analysis_degbg/mapt/mapt_sig_regulons.rds")

mapt_up <- intersect(names(mapt_up_bg), names(mapt_up_degbg))
mapt_up_reg <- mapt_up_bg[names(mapt_up_bg) %in% mapt_up]

# Down-regulated
mapt_down_bg <- readRDS("regulon_down_analysis_300418/mapt/mapt_sig_regulons.rds")
mapt_down_degbg <- readRDS("regulon_analysis_degbg/mapt_down/mapt_down_regulons.rds")
# no intersection
mapt_down <- intersect(names(mapt_down_bg), names(mapt_down_degbg))


####
# GRN regulons
####
# Up-regulated
grn_up_bg <- readRDS("regulon_analysis_250418/grn/grn_sig_regulons.rds")
grn_up_degbg <- readRDS("regulon_analysis_degbg/grn/grn_sig_regulons.rds")
grn_up <- intersect(names(grn_up_bg), names(grn_up_degbg))

# Down-regulated
grn_down_bg <- readRDS("regulon_down_analysis_300418/grn/grn_sig_regulons.rds")
grn_down_degbg <- readRDS("regulon_analysis_degbg/grn_down/grn_sig_regulons.rds")
grn_down <- intersect(names(grn_down_bg), names(grn_down_degbg))


####
# C9orf72 regulons
####
# Up-regulated
c9_up_bg <- readRDS("regulon_analysis_250418/c9orf72/c9orf72_regulons.rds")
c9_up_degbg <- readRDS("regulon_analysis_degbg/c9orf72/c9orf72_regulons.rds")
c9_up <- intersect(names(c9_up_bg), names(c9_up_degbg))

# Down-regulated
c9_down_bg <- readRDS("regulon_down_analysis_300418/c9orf72/c9orf72_regulons.rds")
d9_down_degbg <- readRDS("regulon_analysis_degbg/c9orf72_down/c9orf72_down_regulons.rds")
c9_down <- intersect(names(c9_down_bg), names(c9_down_bg))


# Create dataframe for all mutations
max_len <- max(length(mapt_up), length(mapt_down), length(grn_up), length(grn_down))
# Fill up
mapt_up <- c(mapt_up, rep("", max_len-length(mapt_up)))
mapt_down <- c(mapt_down, rep("", max_len-length(mapt_down)))
grn_up <- c(grn_up, rep("", max_len-length(grn_up)))
grn_down <- c(grn_down, rep("", max_len-length(grn_down)))

df <- data.frame(mapt_up, grn_up, mapt_down, grn_down)

write.table(df, "high_confidence_TFs_MAPT_GRN.txt", sep="\t", quote = F, row.names = F)


# Test to intersect with chipseq data
# Mapt
mup <- read.table("regulon_analysis_degbg/mapt/results/encode/extended_peaks_seqs.merged.encodeIntersect.txt", sep="\t", header = T)
mdown <- read.table("regulon_analysis_degbg/mapt_down/results/encode/extended_peaks_seqs.merged.encodeIntersect.txt", sep="\t", header = T)
# Test for upTFs
mup_test <- mup[mup$TF %in% mapt_up,]
mdown_test <- mdown[mdown$TF %in% mapt_up,]
