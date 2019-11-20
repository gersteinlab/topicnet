nmf_data = read.table('NMF_50.txt',header = T)
nmf_data = t(nmf_data)
library(Hmisc)
correlation = rcorr(nmf_data,type = 'pearson')
nmf_corr = correlation$r
filename = 'nmf_corr.txt'
write.table(nmf_corr,filename,row.names = F, sep = '\t')
