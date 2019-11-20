args <- commandArgs(TRUE)
rank <- as.integer(args[1])
method <- 'brunet'
seed = 'random'

library(NMF)

data = read.table("/ysm-gpfs/pi/gerstein/xk4/TFmembership/edge_allencode.mat.txt",row.names = 1,header = TRUE)
data = t(data)
nmf_model = nmf(data,rank,method,seed)
#nmf_data = coef(nmf_model)
#nmf_data = t(nmf_data)

#filename = paste('NMF_',rank,'.txt',sep = '')
#write.table(nmf_data, filename, sep = '\t', row.names = F)

#nmf_basis = basis(nmf_model)
#nmf_basis = t(nmf_basis)
#filename = paste('NMF_basis_',rank,'.txt',sep = '')
#write.table(nmf_basis,filename,sep = '\t',row.names = F)

filename = paste('nmf_',rank,'.rdata',sep = '')

