####
# Flipper population size plot
####
library(ggplot2)
source("~/gitRepos/analysis_scripts/kevin/flipper/psmc_plotting.R")

setwd("~/flipper/psmc_analysis/")

# Tursiops
df.tt <- psmc.result("tursiops/diploid_tt_203_t20.psmc", mu=1.5e-8, g = 21.1)
df.tt$Ne <- df.tt$Ne / 10000
df.tt$Species <- rep("T. truncatus", nrow(df.tt))

# Stenella
df.s <- psmc.result("stenella/diploid_sc_214.psmc", mu=1.5e-8, g = 22.5)
df.s$Ne <- df.s$Ne / 10000
df.s$Species <- rep("S. coerulleoalba", nrow(df.s))

# Grampus
df.g <- psmc.result("grampus/diploid_gg_225.psmc", mu=1.5e-8, g = 19.6)
df.g$Ne <- df.g$Ne / 10000
df.g$Species <- rep("G. griseus", nrow(df.g))

df <- rbind(df.tt, df.s, df.g)

ggplot(df, aes(x=YearsAgo, y=Ne, color=Species)) +
  geom_line(size=1) +
  scale_x_log10() + annotation_logticks(side="b") +
  theme_minimal() +
  xlab("Years ago") + ylab("Effective populationo size (10^4)")

ggsave("psmc_plot.png", width=10, height = 4)
