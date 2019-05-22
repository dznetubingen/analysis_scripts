#########################
## Generate PPI files from STRING
#############################

library(biomaRt)
# Load PPI information
ints <- read.table("~/resources/ppi/9606.protein.actions.v10.5.txt", header=T, fill = T)
md <- read.table("~/resources/ppi/9606.protein.aliases.v10.5.txt", fill=T)

ints <- ints[!is.na(ints$score),]
ints <- ints[ints$score >= 700,]
ints$item_id_a <- gsub("9606.", "", ints$item_id_a)
ints$item_id_b <- gsub("9606.", "", ints$item_id_b)
ppi <- read.table("~/resources/ppi/9606.protein.links.detailed.v10.5.txt", header=T)
ppi <- ppi[ppi$combined_score >= 700,] # only consider interactions with high scores
ppi$protein1 <- gsub("9606.", "", ppi$protein1)
ppi$protein2 <- gsub("9606.", "", ppi$protein2)

# Get corresponding symbols
ensembl <- useMart("ensembl", "hsapiens_gene_ensembl")
bm <- getBM(mart = ensembl, attributes = c("ensembl_peptide_id", "hgnc_symbol"), values = ints$item_id_a, filters="ensembl_peptide_id")
ints <- merge(ints, bm, by.x="item_id_a", by.y="ensembl_peptide_id")
colnames(ints)[8] <- "hgnc_symbol_a"
bm <- getBM(mart = ensembl, attributes = c("ensembl_peptide_id", "hgnc_symbol"), values = ints$item_id_b, filters="ensembl_peptide_id")
ints <- merge(ints, bm, by.x = "item_id_b", by.y="ensembl_peptide_id")
colnames(ints[9]) <- "hgnc_symbol_b"
ints <- ints[!duplicated(ints),]
# Save the table for later usage!
write.table(ints, "~/resources/ppi/adjusted_int_table_900_symbols.txt", sep="\t", quote=F, row.names = F)
# Load earlier created interaction table
