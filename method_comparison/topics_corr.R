topics = read.table('/ysm-gpfs/pi/gerstein/xk4/TFmembership/LDAmodels/50_topics.txt',header = T)
topics_corr = cor(topics,method = 'pearson')
write.table(topics_corr,'topics_corr.txt',row.names = F,sep = '\t')

