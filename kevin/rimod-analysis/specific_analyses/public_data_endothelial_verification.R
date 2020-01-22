####################
# Analysis of deconvolution results of public data set (Nature Medicine, swarup, 2019)
###################
library(ggplot2)
library(reshape2)

setwd("~/rimod/public_data/mRNAseq_syn7820574/")


# Load cortex data
fracs <- read.table("scaden_predictions_cortex.txt", sep="\t", header=T, row.names = 1)
md <- read.table("metadata_HumanFTD_Cortex.txt", sep="\t", header=T)

md <- md[match(rownames(fracs), md$SampleID),]
fracs$dc <- md$Diagnosis

tau <- fracs[fracs$dc == "FTD-Tau",]
tauneg <- fracs[fracs$dc == "FTD-TauNeg",]
control <- fracs[fracs$dc == "Control",]


####
# Calculate average increase/decrase
####
mean_fun <- median
fracs <- fracs[, -1]
tau <- fracs[fracs$dc == "FTD-Tau",]
tauneg <- fracs[fracs$dc == "FTD-TauNeg",]
control <- fracs[fracs$dc == "Control",]

# calculate averages
tau <- tau[, -ncol(tau)]
tau <- apply(tau, 2, mean_fun)
tauneg <- tauneg[, -ncol(tauneg)]
tauneg <- apply(tauneg, 2, mean_fun)
control <- control[, -ncol(control)]
control <- apply(control, 2, mean_fun)


# Tau
tau_list <- c()
for (i in 1:length(control)) {
  ct.mean <- control[i]
  diff <- tau[i] - ct.mean
  pct <- diff / ct.mean
  print(pct)
  tau_list <- c(tau_list, pct)
}

# TauNeg
tauneg_list <- c()
for (i in 1:length(control)) {
  ct.mean <- control[i]
  diff <- tauneg[i] - ct.mean
  pct <- diff / ct.mean
  print(pct)
  tauneg_list <- c(tauneg_list, pct)
}


# make data frame
tau <- t(data.frame(tau_list))
tauneg <- t(data.frame(tauneg_list))
df <- data.frame(rbind(tau, tauneg))
df$group <- c("Tau", "TauNeg")


# plotting
df <- melt(df)

ggplot(data=df, aes(x=variable, y=value, fill=group)) +
  geom_bar(stat='identity', position = position_dodge()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, size = 12)) +
  xlab("") + 
  ylab("Percentage difference to control")
