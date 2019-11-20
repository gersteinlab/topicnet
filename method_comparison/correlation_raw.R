data = read.table("/ysm-gpfs/pi/gerstein/xk4/TFmembership/edge_allencode.mat.txt",row.names = 1,header = TRUE)

library(Hmisc)
data = t(data)
correlation = rcorr(data,type = 'pearson')
raw_corr = correlation$r

filename = 'raw_corr.txt'
write.table(raw_corr,filename,row.names = F, sep = '\t')


