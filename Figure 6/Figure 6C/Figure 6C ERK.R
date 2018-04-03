library(ggplot2)

fig6cErk = as.data.frame(read.csv("/Users/mary/GitHub/Kovary_MSB_2018/Figure 6/Figure 6C/Fig6cErk.csv"))

auc = apply(fig6cErk[,5:46], 1, function(x) sum(x) - min(x))

num = round(length(auc)*.15)
high = as.numeric(auc[order(auc)[length(auc)-num]])
low = as.numeric(auc[order(auc)[num]])

fig6cErk$AUC = "Mid"
fig6cErk$AUC[which(auc>=high)] = "High"
fig6cErk$AUC[which(auc<=low)] = "Low"

ggplot(fig6cErk[fig6cErk$AUC %in% c("High","Low"),], aes(AUC, ERK_Norm)) + 
  geom_boxplot(notch = T, outlier.shape = NA) + theme_classic() + ylim(1.5, 3.5)

