# General analysis of the time-course data

### lda.t50.t2t: topic to gene matrix
### refseq.groseq.vec: groseq measure

load("TFmem_data.fig6.RData")
load("med.fgsea.rdata")
load("fig6.groseq.expr.rdata")
load("TFmem_data.fig6-20190520.RData")

c2all = read.table("c2all.0.25.0.txt", header = TRUE)

require(fgsea)
require(gplots)
require(tsne)

########## Time course change
### Plot time course change
t_order = colnames(refseq.groseq.vec)

ESR1.list = paste(c("SRX469393","SRX469394","SRX469395","SRX469396","SRX469397","SRX469398"),"MCF.7:human:ESR1", sep=".")
time = factor(c("0","2","5","10","40","160"), levels=c("0","2","5","10","40","160"))
topics.mcf7.esr1 = topics.mcf7[ESR1.list,]

pdf('fig6/fig6.heatmap.pdf')
heatmap.2(t(topics.mcf7.esr1), Colv = FALSE, labCol = paste(time,"min",sep=""), trace="none")
dev.off()

plot.topic = paste("X",c(3,4,31,43),sep="")
plot.topic.id = c(3,4,31,43)
n.topic = length(plot.topic)

#pdf('/Users/tianxiao/Documents/Lab/TF membership/fig6/plot.timecourse/fig6.timecourse.pdf',width=8, height=6)

#par(mfrow=c(2,1), mar = c(2,4,2,4))
#mp=matplot(topics.mcf7.esr1[,plot.topic], type="b", xlim=c(0.8,6), ylim=c(0,0.2), 
#           xlab="time",ylab="module proportion", pch=1, xaxt="n", lty=2)
#ax=axis(1, at=1:6, labels=time)
#legend("topleft", legend=plot.topic, col=1:n.topic, pch=1, cex=0.8)

esr1.groseq = refseq.groseq.vec["ESR1", t_order]
esr1.groseq.vec = c(esr1.groseq[,1], 0, 0, esr1.groseq[,2], esr1.groseq[,3], esr1.groseq[,4])
ylim = c(max(esr1.groseq), min(esr1.groseq))
#barplot(esr1.groseq.vec, col="purple", xaxt="n", ylim=ylim)
#dev.off()

### Timecourse plotted separately
pdf('fig6/plot.timecourse/plot.timecourse.sep.pdf', width=6, height=10)
par(mfrow=c(4,1))
for(i in plot.topic.id){
  plot(as.numeric(time), topics.mcf7.esr1[,i], type="b", col="blue", ylim=c(0,0.2),
       xlab="time", ylab="module proportion", main = paste("X",i,sep=""), pch=1, xaxt="n")
  axis(1, at=1:6, labels=time)
}
dev.off()



########## Annotation using top genes
### Gene set enrichment
for(i in plot.topic.id){
  cat(i,"\n")
  write(as.character(head(c2all[c2all$tid==i,]$pathway, n=20)), file=paste("c2_topic",i,".txt",sep=""))
}


### Top weighted genes for the topics
get_top_genes <- function(mat, i, n=20){
  v = mat[,i]
  v = v[names(v) %in% rownames(refseq.groseq.vec)]
  return(names(v)[order(v, decreasing = TRUE)][1:n])
}


top20.genes.list = list()
top1000.genes.list = list()
#for(i in seq(1,length(plot.topic.id))){
for(i in seq(1,50)){
    #i = plot.topic.id[j]
  top20.genes.i = get_top_genes(lda.t50.t2t, i, n=20)
  top20.genes.list[[i]] = top20.genes.i
  
  top1000.genes.i = get_top_genes(lda.t50.t2t, i, n=1000)
  top1000.genes.list[[i]] = top1000.genes.i
  #pdf(paste('/Users/tianxiao/Documents/Lab/TF membership/fig6/gro-seq/groseq_topic',i,'.pdf',sep=''))
  #heatmap.2(as.matrix(refseq.groseq.vec[top.genes.i,t_order]), Rowv=FALSE, Colv=FALSE, dendrogram = "none", trace="none", cexCol = 0.8, cexRow = 0.8)
  #dev.off()
}
names(top20.genes.list) = paste("X", seq(1,50), sep="")
names(top1000.genes.list) = paste("X", seq(1,50), sep="")



########## GRO-seq analysis
gene.var = apply(refseq.groseq.vec, FUN=var, 1)
refseq.groseq.var = refseq.groseq.vec[order(gene.var, decreasing = TRUE)[1:1000],]
#refseq.groseq.var = refseq.groseq.vec[(gene.var>0.1),]

heatmap.2(as.matrix(refseq.groseq.var), Colv=FALSE, scale = "row", col = c,
          trace="none", cexCol = 0.8, cexRow = 0.8)

### GRO-seq correlation with ERS1
gro.esr1 = refseq.groseq.vec["ESR1",]
n.gene = nrow(refseq.groseq.vec)
corr.with.esr1 = rep(0, n.gene)
names(corr.with.esr1) = rownames(refseq.groseq.vec)

for(i in seq(1, n.gene)){
  corr.with.esr1[i] = cor(t(gro.esr1), t(refseq.groseq.vec[i,]))
}

corr.with.esr1 = corr.with.esr1[!is.na(corr.with.esr1)]

corr.with.esr1.pos = corr.with.esr1[corr.with.esr1>0.8]
write(names(corr.with.esr1.pos), file="corr.list/corr_list_pos_0.8_all.txt")

corr.with.esr1.neg = corr.with.esr1[corr.with.esr1< -0.8]
write(names(corr.with.esr1.neg), file="corr.list/corr_list_neg_0.8_all.txt")

#corr.with.esr1.top20.list = list()
corr.with.esr1.top1000.list = list()

topic_plot=c(3,4,31,43)
for(i in topic_plot){
  #gset = top20.genes.list[[i]]
  #corr.with.esr1.top20.list[[i]] = corr.with.esr1[gset]
  gset = top1000.genes.list[[i]]
  write(names(corr.with.esr1.pos[gset[gset %in% names(corr.with.esr1.pos)]]), file=paste("corr.list/corr_list_pos_0.8_topic",i,".txt", sep=""))
  write(names(corr.with.esr1.neg[gset[gset %in% names(corr.with.esr1.neg)]]), file=paste("corr.list/corr_list_neg_0.8_topic",i,".txt", sep=""))
}

#names(corr.with.esr1.top20.list) = names(top20.genes.list)
#names(corr.with.esr1.top200.list) = names(top200.genes.list)
#corr.with.esr1.top.abs.list = lapply(corr.with.esr1.top.list, abs)

#corr.with.esr1.top.vec = unlist(corr.with.esr1.top.list)

#corr.with.esr1.top.abs.high.list = lapply(corr.with.esr1.top.list, function(v){names(v[v>0.7])}) # highly correlated
#cat(unlist(corr.with.esr1.top.abs.high.list))


### Plotted separately
mean.corr = mean(corr.with.esr1[!is.na(corr.with.esr1)])
topic.all.corr = mean(unlist(lapply(corr.with.esr1.top200.list, function(x){mean(abs(x)[!is.na(x)])})))
  
pdf('fig6/plot.timecourse/top.genes.corr.sep.pdf', width=4, height=14)
par(mfrow=c(4,1), mar=c(2,6,2,6), xpd=T)
col.panel = c("darkblue", "lightblue")
for(t in plot.topic){
  v = corr.with.esr1.top20.list[[t]]
  v.sign = ((v>0)+1)
  barplot(abs(v), col = col.panel[v.sign],  horiz = T, las=1)
  lines(x=c(mean.corr,mean.corr), y=c(-3,26), lty=2, col="red")
  lines(x=c(topic.all.corr,topic.all.corr), y=c(-3,26), lty=2, col="blue")
}
dev.off()

########## ##########
### Plotted together
pdf('fig6/plot.timecourse/top.gene.corr.pdf',width=10, height=6)
bar.col = c(rep("red",20), rep("green",20), rep("blue",20), rep("orange",20))
bp = barplot(unlist(corr.with.esr1.top20.list), xaxt='n', col = bar.col, ylim=c(-1.2,1),space=0.8)
abline(v=bp[c(20,40,60),]+0.8, lty=2)
text(x = bp[c(10,30,50,70)+0.4,], y =-1.1, labels=plot.topic)
dev.off()

### Boxplot of 100 top genes
top50.genes.list = list()
for(j in seq(1,length(plot.topic.id))){
  i = plot.topic.id[j]
  top.genes.i = get_top_genes(lda.t50.t2t, i, n = 50)
  top50.genes.list[[j]] = top.genes.i
}
names(top50.genes.list) = plot.topic


corr.with.esr1.top50.list = list()
for(i in seq(1, length(plot.topic))){
  gset = plot.topic[[i]]
  corr.with.esr1.top50.list[[i]] = corr.with.esr1[gset]
}

names(corr.with.esr1.top50.list) = names(top.genes.list)
corr.with.esr1.top50.list[["rand"]] = corr.with.esr1[sample(rownames(refseq.groseq.vec), 50)]
corr.with.esr1.top.abs.list = lapply(corr.with.esr1.top50.list, abs)

pdf('fig6/plot.timecourse/top50.gene.corr.pdf',width=10, height=6)
boxplot(corr.with.esr1.top.abs.list)
dev.off()


corr.with.esr1.top50.list[["rand"]] = NULL
corr.with.esr1.top50.abs.high.list = lapply(corr.with.esr1.top50.list, function(v){names(v[v>0.6])}) # highly correlated
corr.with.esr1.top50.abs.high.all = unlist(corr.with.esr1.top50.abs.high.list)

corr.with.esr1.top50.abs.high.all = corr.with.esr1.top50.abs.high.all[!is.na(corr.with.esr1.top50.abs.high.all)]

