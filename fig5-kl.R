# KL divergence between time points

load("TFmem_data.fig6-20190520.RData")

library(gplots)
library(entropy)

### Pairwise KL divergence
pair.mat <- function(m, f, symm=FALSE){ # pairwise distance or similarity between rows of m given the distance function
  n = nrow(m)
  d = matrix(0L, nrow=n, ncol=n)
  for(i in seq(1,n)){
    for(j in seq(i,n)){
      vi = m[i,]
      vj = m[j,]
      if(symm) {
        sym = (f(vi, vj) + f(vj, vi))/2
        d[i, j] = sym
        d[j, i] = sym
      }
      else {
        d[i, j] = f(vi, vj)
        d[j, i] = f(vj, vi)
      }
    }
  }
  rownames(d) = rownames(m)
  colnames(d) = rownames(m)
  return(d)
}

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

get_upper_tri<-function(cormat){
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}

topics.mcf7.esr1.kl = pair.mat(topics.mcf7.esr1, KL.plugin)
topics.mcf7.esr1.kl.symm = get_upper_tri((topics.mcf7.esr1.kl + t(topics.mcf7.esr1.kl))/2)
diag(topics.mcf7.esr1.kl.symm)=NA
colnames(topics.mcf7.esr1.kl.symm) = paste(time,"min",sep="")
rownames(topics.mcf7.esr1.kl.symm) = paste(time,"min",sep="")
topics.mcf7.esr1.kl.symm = topics.mcf7.esr1.kl.symm[,rev(paste(time,"min",sep=""))]

pdf(paste("/Users/tianxiao/Documents/Lab/TF membership/fig6/esr1-klsymm.pdf", sep=""), height = 10, width = 10)
heatmap.2(topics.mcf7.esr1.kl.symm, scale="none", symm=TRUE, trace="none", col = colorRampPalette(c("steelblue","gray95"))(32), 
          main="ESR1", Rowv = FALSE,
          distfun = function(x) as.dist(x)) # cluster by the sqrt distance matrix
dev.off()

