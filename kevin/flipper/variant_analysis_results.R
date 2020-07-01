###########
# Flipper calculate SNP ratios
###########
setwd("~/flipper/")
library(ggplot2)

# Tursiops
tref <- 3326478
t1 <- 5369212
t2 <- 5123166
tursiops <- mean(c(t1 / tref, t2 / tref))
print(tursiops)

# Stenella
sref <- 7095215
s1 <- 10757121
s2 <- 10665600
stenella <- mean(c(s1 / sref, s2 / sref))
print(stenella)

# Grampus
gref <- 4009395
g1 <- 4927483
grampus <- g1 / gref
print(grampus)


df <- data.frame(SNPs = c(tref, t1, t2, sref, s1, s2, gref, g1), 
                 Species = c("Tursiops", "Tursiops", "Tursiops", "Stenella", "Stenella", "Stenella", "Grampus", "Grampus"),
                 Type = c("Reference", "NonReference", "NonReference", "Reference", "NonReference", "NonReference", "Reference", "NonReference"),
                 Individual = c("TRef", "T1", "T2", "SRef", "S1", "S2", "GRef", "G1"))

ggplot(df, aes(x=Individual, y = SNPs, fill = Species)) +
  geom_bar(stat="identity") +
  theme_minimal()

ggsave("barplot_SNP_numbers.png", width=5, height=5)
