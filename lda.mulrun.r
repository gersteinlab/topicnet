


library(topicmodels)
#library(ldatuning)
#library(parallel)

args<-commandArgs(TRUE)

ntopic <- as.integer(args[1])
fin = args[2]
fout=args[3]
nn = as.integer(args[4])




data = read.table(fin, sep="\t", header=T, row.names=1, stringsAsFactors=F)

data.matrix=as.matrix(data)


num.ds = nrow(data)

out.list=list()

for (idx in 1:nn){

  tmp.lda=LDA(data.matrix,k=ntopic, method="Gibbs")
  out.list[[idx]]=tmp.lda
  
}


save(out.list, file=paste(fout,".rdata",sep=""))
