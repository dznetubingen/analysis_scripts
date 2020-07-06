####
# Flipper population size plot
####
source("~/gitRepos/analysis_scripts/kevin/flipper/psmc_plotting.R")

setwd("~/flipper/psmc_analysis/")

# Tursiops
df.tt <- psmc.result("tursiops/diploid_tt_203_t20.psmc", mu=1.5e-8, g = 21.1)
#df.tt$Ne <- df.tt$Ne / 10000
df.tt$Species <- rep("T. truncatus", nrow(df.tt))

# Stenella
df.s <- psmc.result("stenella/diploid_sc_214.psmc", mu=1.5e-8, g = 19.6)
#df.s$Ne <- df.s$Ne / 10000
df.s$Species <- rep("S. coerulleoalba", nrow(df.s))

df <- rbind(df.tt, df.s)

ggplot(df, aes(x=YearsAgo, y=Ne, color=Species)) +
  geom_line(size=2) +
  scale_x_log10() + annotation_logticks(side="b") +
  theme_minimal() +
  xlab("Years ago") + ylab("Effective populationo size (10^4)")

ggsave("psmc_plot.png", width=10, height = 4)
