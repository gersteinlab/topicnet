args<-commandArgs(TRUE)

n_start <- as.integer(args[1])
n_end <- as.integer(args[2])
step <- as.integer(args[3])

library(topicmodels)
library(ldatuning)

data <- read.table('/ysm-gpfs/pi/gerstein/xk4/TFmembership/edge_allencode.mat.txt')
data.matrix <- unname(data)

print("Finding topics")
find.topics <- FindTopicsNumber(data.matrix, topics = seq(n_start,n_end,step), metrics = c("Arun2010", "CaoJuan2009", "Griffiths2004"))

print("Write to file")
find.topics <- data.frame(find.topics)
filename <- paste("find.topics",".",n_start,"-",n_end,".txt",sep="")
write.table(find.topics,filename,sep='\t',row.names = F)

FindTopicsNumber_plot(find.topics)


