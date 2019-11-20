data = read.table('/ysm-gpfs/pi/gerstein/xk4/TFmembership/method_comparison/kmeans_prob_dist3.txt',header = T)
data = t(data)
library(Hmisc)
correlation = rcorr(data,type = 'pearson')
kmeans_corr = correlation$r
write.table(kmeans_corr,'kmeans_corr.txt',row.names = F,sep = '\t')

