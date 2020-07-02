####
# Flipper population size plot
####
source("~/gitRepos/analysis_scripts/kevin/flipper/psmc_plotting.R")

setwd("~/flipper/psmc_analysis/")

test$Ne <- test$Ne / 1000

df.tt <- psmc.result("tursiops/diploid_tt_203.psmc", mu=1.5e-8, g = 21.1)
df.tt$Ne <- df.tt$Ne / 10000

ggplot(df.tt, aes(x=YearsAgo, y=Ne)) +
  geom_line() +
  scale_x_log10() + annotation_logticks(side="b") +
  theme_minimal()
