library(ggplot2)

fig7eErk = as.data.frame(read.csv("/Volumes/labdata/mary/Manuscripts_In_Progress (teruel1@stanford.edu)/Xenopus_Paper/final_analysis_scripts/Figure 7/Figure 7E/Fig7eErk.csv"))

auc = apply(fig7eErk[,5:46], 1, function(x) sum(x) - min(x))

num = round(length(auc)*.15)
high = as.numeric(auc[order(auc)[length(auc)-num]])
low = as.numeric(auc[order(auc)[num]])

fig7eErk$AUC = "Mid"
fig7eErk$AUC[which(auc>=high)] = "High"
fig7eErk$AUC[which(auc<=low)] = "Low"

ggplot(fig7eErk[fig7eErk$AUC %in% c("High","Low"),], aes(AUC, ERK_Norm)) + 
  geom_boxplot(notch = T, outlier.shape = NA) + theme_classic() + ylim(1.5, 3.5)
