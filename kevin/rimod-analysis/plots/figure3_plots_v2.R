#####
# Visualizatin of HumanBase modules
####
setwd("~/rimod/RNAseq/analysis/human_base/")


###
# Extract subset of modules
###
subsetModules <- function(df, n_mods=3){
  mods <- as.character(levels(df$CLUSTER_NAME))
  
  # Initial df
  mod <- df[df$CLUSTER_NAME == "M1",]
  mod <- mod[order(mod$TERM_Q_VALUE),]
  mod <- mod[1:n_mods,]
  mod.df <- mod
  
  for (m in mods) {
    mod <- df[df$CLUSTER_NAME == m,]
    mod <- mod[order(mod$TERM_Q_VALUE),]
    mod <- mod[1:n_mods,]
    mod.df <- rbind(mod.df, mod)
  }
  mod.df <- mod.df[c(-1:-n_mods),] # remove initial
  return(na.omit(mod.df))
}


# GRN down
df <- read.table("rnaseq_grn_filtered_down_enrichment.txt", sep="\t", header=T)
grn.down <- subsetModules(df)

# GRN up
df <- read.table("rnaseq_grn_filtered_up_enrichment.txt", sep="\t", header=T)
grn.up <- subsetModules(df)


# MAPT down
df <- read.table("rnaseq_mapt_filtered_down_enrichment.txt", sep="\t", header=T)
mapt.down <- subsetModules(df)

# MAPT up
df <- read.table("rnaseq_mapt_filtered_up_enrichment.txt", sep="\t", header=T)
mapt.up <- subsetModules(df)
