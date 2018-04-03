library(ggplot2)
library(tidyr)
library(plyr)

# Load Data
fig6CD = as.data.frame(read.csv("/Volumes/labdata/mary/Manuscripts_In_Progress (teruel1@stanford.edu)/Xenopus_Paper/final_analysis_scripts/Figure 6/Fig6CD.csv"))
colnames(fig6CD)[3:ncol(fig6CD)] = seq(10, 78, 2)

# Reformat
df_cells = gather(fig6CD, "Time","FRET", 5:ncol(fig6CD))
df_cells$Time = as.numeric(df_cells$Time)

# Figure 6A Timecourse
ggplot(df_cells[df_cells$Cell %in% sample(1:nrow(cells_temp), 300),], aes(Time, FRET, group = Cell)) + geom_line(size=1) + 
  ylim(1,2.5) + theme_classic() + facet_grid(EGF ~ .) + theme(legend.position="none")

# Figure 6A Histograms
AUC = ddply(df_cells, .(Cell, EGF), summarize, AUC = sum(FRET))
ggplot(AUC, aes(AUC)) + geom_histogram(bins = 80) + facet_grid(EGF~.) + 
  theme_classic() + geom_vline(xintercept = 45)

# Figure 6B
percent_response = ddply(AUC, .(EGF), summarize, Percent = sum((AUC >= 45)) / length(AUC))
ggplot(percent_response, aes(x = log2(EGF), y = Percent)) + geom_point(size = 5) +
  theme_classic()
