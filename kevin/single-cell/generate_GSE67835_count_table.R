############################
# Merge GSE67835 count data
############################
setwd("~/deepcell_project/revision1_september19/rosmap_deconvolution/training_data/GSE67835_RAW/")

files = list.files(pattern = "*.csv")

# starting dataframe
df <- read.table(files[1], sep="\t", row.names = 1)
colnames(df)[1] <- gsub(".csv", "", files[1])

for (i in 2:length(files)) {
  f <- files[i]
  fname <- gsub(".csv", "", f)
  tmp <- read.table(f, sep="\t", row.names=1)
  colnames(tmp)[1] <- fname 
  df <- merge(df, tmp, by="row.names")
  rownames(df) <- df[,1]
  df <- df[,-1]
}

write.table(df, "../GSE67835_count_table.txt", sep="\t", quote=F, col.names = NA)
