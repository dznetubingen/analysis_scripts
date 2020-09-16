md <- read.csv("~/rimod/files/FTD_Brain_corrected.csv")

md <- md[!duplicated(md$SAMPLEID),]
table(md$DISEASE.CODE)

mapt <- md[md$DISEASE.CODE == "FTD-MAPT",]
grn <- md[md$DISEASE.CODE == "FTD-GRN",]
c9 <- md[md$DISEASE.CODE == "FTD-C9",]
control <- md[md$DISEASE.CODE == "control",]


groups <- c("FTD-MAPT", "FTD-GRN", "FTD-C9", "control")


for (g in groups) {
  print(g)
  df <- md[md$DISEASE.CODE == g,]
  print("age")
  print(mean(na.omit(df$AGE)))
  print(sd(na.omit(df$AGE)))
  
  print("PMI")
  print(mean(na.omit(as.numeric(as.character(na.omit(df$PMD.MIN.))))))
  print(sd(na.omit(as.numeric(as.character(na.omit(df$PMD.MIN.))))))
  
  
  print("ph")
  print(mean(na.omit(as.numeric(as.character(na.omit(df$PH))))))
  print(sd(na.omit(as.numeric(as.character(na.omit(df$PH))))))
  
  print("RIN")
  print(mean(as.numeric(as.character(na.omit(df$RIN)))))
  print(sd(as.numeric(as.character(na.omit(df$RIN)))))
  
  print("males")
  print(table(df$GENDER))
  }
