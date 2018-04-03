library(ggplot2)

fig7eMek = as.data.frame(read.csv("/Volumes/labdata/mary/Manuscripts_In_Progress (teruel1@stanford.edu)/Xenopus_Paper/final_analysis_scripts/Figure 7/Figure 7E/Fig7eMek.csv"))

auc = apply(fig7eMek[,5:38], 1, function(x) sum(x) - min(x))

num = round(length(auc)*.15)
high = as.numeric(auc[order(auc)[length(auc)-num]])
low = as.numeric(auc[order(auc)[num]])

fig7eMek$AUC = "Mid"
fig7eMek$AUC[which(auc>=high)] = "High"
fig7eMek$AUC[which(auc<=low)] = "Low"

ggplot(fig7eMek[fig7eMek$AUC %in% c("High","Low"),], aes(AUC, MEK_Norm)) + 
  geom_boxplot(notch = T, outlier.shape = NA) + theme_classic() + ylim(0.2, 0.45)
