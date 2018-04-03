library(ggplot2)
library(tidyr)
library(plyr)

# Load Data
fig6A = read.csv("/Volumes/labdata/mary/Manuscripts_In_Progress (teruel1@stanford.edu)/Xenopus_Paper/final_analysis_scripts/Figure 6/Fig6A.csv")
colnames(fig6A)[3:ncol(fig6A)] = c(0:4800)

# Reshape Data
df_cells = gather(fig6A, "Time","FRET", 3:(ncol(fig6A)))
df_cells$Time = as.numeric(df_cells$Time)
df_cells$Condition = factor(df_cells$Condition, levels=c('Low','Mid','High'))

# Plot Figure 6A
ggplot(df_cells, aes(Time, FRET, group = Cell)) + geom_line(alpha=.5, size=0.5) + 
  theme_classic() + facet_grid(Condition ~ .) + theme(legend.position="none")

# AUC Histogram
df_hist = ddply(df_cells, .(Condition, Cell), summarize, AUC = sum(FRET))

ggplot(df_hist, aes(AUC)) + geom_histogram() + theme_classic() + 
  facet_grid(Condition ~ .)
