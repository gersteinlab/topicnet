args <- commandArgs(TRUE)
k <- as.integer(args[1])

data = read.table("/ysm-gpfs/pi/gerstein/xk4/TFmembership/edge_allencode.mat.txt",row.names = 1,header = TRUE)
kcluster = kmeans(data,k,iter.max = 1000,nstart = 20)

kmeans_center = kcluster$centers

library(pdist)
distance = as.matrix(pdist(data,kmeans_center))

#Convert distance to probability
#dist =  1.0/distance
#cluster_prob = t(apply(dist,1,function(x) x/sum(x)))
#write.table(cluster_prob,file="kmeans_prob_dist1.txt",row.names = F,sep = '\t')

#Correlation
#data = t(cluster_prob)
library(Hmisc)
correlation = rcorr(t(distance),type = 'pearson')
kmeans_corr = correlation$r
write.table(kmeans_corr,'kmeans_corr.txt',row.names = F,sep = '\t')

