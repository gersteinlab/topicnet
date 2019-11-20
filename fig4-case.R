# Rewiring case study: GM12878 vs K562
# Identify topic rewiring between GM12878 and K562

load("case_TFmem_data.fig4-190904.RData")

library(reshape2)
library(entropy)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(gplots)
library(ggplot2)

top.genes.list = list()

# Set up tfs and cells to study
tf='EP300'
#tf='ZBTB33'
cell1='GM12878'
cell2='K562'

doc.row1 = paste(cell1, tf, sep=":")
doc.row2 = paste(cell2, tf, sep=":")

topics1 = topics.50[doc.row1,]
topics2 = topics.50[doc.row2,]

data1 = data[doc.row1,]
data2 = data[doc.row2,]

tf.topic50.case = tf.topic50.cell.klsymm.df[(tf.topic50.cell.klsymm.df$`Cell-1` == cell1)*
                                              (tf.topic50.cell.klsymm.df$`Cell-2` == cell2)>0,]

rew.diff = topics1 - topics2
top.rew.topics = order(abs(rew.diff), decreasing = TRUE)[1:4]

###### Plot top contributed genes for all topics
get_top_genes <- function(mat, i, n=20){
  v = mat[,i]
  return(names(v)[order(v, decreasing = TRUE)][1:n])
}

col1 = rgb(0, 0, 255, max = 255, alpha = 255/2, names = "blue50")
col2 = rgb(255, 0, 0, max = 255, alpha = 255/2, names = "red50")
col.list = c(col1,col2)
names(col.list) = c(cell1, cell2)

for(i in seq(1,length(top.rew.topics))) {
  t = top.rew.topics[i]
  diff = rew.diff[t]
  if(diff>0){
    color = col.list[cell1]
    target = data1
  }else{
    color = col.list[cell2]
    target = data2
  }
  
  name = paste(tf,cell1,cell2,sep="-")
  
  pdf(paste("case.gm-k562/top.genes/", name, "_rew_", i, "_topic_", t, ".pdf", sep=""), 
      height=8, width=4.5)
  par(mai=c(2,2,2,2))
  top.genes = get_top_genes(lda.t50.t2t, t, n=20)
  top.genes.val = lda.t50.t2t[rev(top.genes), t]
  xlim = c(min(top.genes.val)*0.75, max(top.genes.val))
  
  target.id = which(target[top.genes]==1)
  nontarget.id = seq(1,20)[!seq(1,20) %in% target.id]
  bp = barplot(top.genes.val, space=0.2, horiz=TRUE, xlim=xlim, xpd=FALSE, xlab='contribution',
          las=2, main=paste("Topic",t), cex.axis=0.75, col=color, yaxt='n')
  axis(2, at = rev(bp[target.id,]), labels = top.genes[target.id], col.axis = color, las=1, cex.axis=0.75)
  axis(2, at = rev(bp[nontarget.id,]), labels = top.genes[nontarget.id], col.axis = "black", las=1, cex.axis=0.75)
  
  scale = max(top.genes.val)
  legend(scale, bp[1,], legend=c(cell1,cell2), fill=col.list, cex=0.75, xpd=TRUE)
  dev.off()
  
  top.genes.list[[paste("topic_", t, "_", name,sep="")]] = top.genes
}

###### Plot circular topic distribution plot
###### NEW VERSION: only plot topic proportions

for(i in c(seq(1,5), seq(nrow(tf.topic50.case)-4,nrow(tf.topic50.case)))) {
  r = tf.topic50.case[i,]
  
  
  tf = as.character(r[1])
  cell1 = as.character(r[2])
  cell2 = as.character(r[3])
  
  cell.tf1 = paste(cell1, tf, sep=":")
  cell.tf2 = paste(cell2, tf, sep=":")
  
  data.cell1 = data[cell.tf1,]
  data.cell2 = data[cell.tf2,]
  cmp = data.cell1 - data.cell2
  c1.gain = names(cmp)[cmp==1]
  c2.gain = names(cmp)[cmp==-1]
  
  c1.gain.mat=lda.t50.t2t[c1.gain,]
  c2.gain.mat=lda.t50.t2t[c2.gain,]
  
  ### Record gain/loss targets
  tf.list = list(cell1=cell1, cell2=cell2, c1.gain=c1.gain, c2.gain=c2.gain, 
                 c1.gain.mat=c1.gain.mat,
                 c2.gain.mat=c2.gain.mat)
  id = paste(paste(tf, cell1, sep="_"), cell2, sep="_")
  gainloss.gene.list[[id]] = tf.list
  
  topic50.cell.tf1 = topics.50[cell.tf1, ]
  topic50.cell.tf2 = topics.50[cell.tf2, ]
  
  c1.gain.mat.mean = colMeans(c1.gain.mat)
  c2.gain.mat.mean = colMeans(c2.gain.mat)
  
  ### Circular
  ymax = max(c(c1.gain.mat.mean, c2.gain.mat.mean))
  #ymax = signif(max(c(c1.gain.mat.mean, c2.gain.mat.mean)), digits=1)
  
  dat.break = data.frame(id=seq(51,55), prop=rep(NA,5))
  
  dat.plot1 = data.frame(id=seq(1,50), prop=as.numeric(topic50.cell.tf1))
  dat.plot2 = data.frame(id=seq(1,50), prop=as.numeric(topic50.cell.tf2))
  
  dat.plot1 = rbind(dat.plot1, dat.break)
  dat.plot2 = rbind(dat.plot2, dat.break)
  
  label_data = data.frame(id = dat.plot1$id, angle= 90 - 360 * (dat.plot1$id-0.5) /length(dat.plot1$id), label=dat.plot1$id)
  label_data$label[51:55] = NA
  
  p = ggplot(dat.plot1, aes(x=id, y=mean)) +       
    geom_bar(data=dat.plot1, aes(x=id, y=prop), stat="identity", fill=alpha("blue", 0.5)) +
    geom_bar(data=dat.plot2, aes(x=id, y=prop), stat="identity", fill=alpha("red", 0.5)) +
    #geom_line(data=dat.plot1, aes(x=id, y=prop), color=alpha("blue", 0.5)) +
    #geom_line(data=dat.plot2, aes(x=id, y=prop), color=alpha("red", 0.5)) +
    scale_y_continuous(limits=c(-0.85, 1), breaks=seq(0, 1, length.out=5)) + 
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5), 
      plot.margin = unit(rep(1,4), "cm")) +
    coord_polar(start = 0) +
    ggtitle(tf)
  
  p = p + geom_text(data=data.frame(x = rep(0, 5), y = seq(0, 1, length.out=5)), aes(x,y), 
                    label = seq(0, 1, length.out=5), 
                    hjust=1.125, size=4)
  
  p = p + geom_text(data=label_data, aes(x=id, y=-0.04, label=label), 
                    size=4, angle= label_data$angle, 
                    inherit.aes = FALSE, hjust = 0.75 )
  
  ggsave(filename=paste("case.gm-k562/topic.distr/plot_",i, "_", id, "_barplot.pdf", sep=""), plot=p, width=8, height=8)
  
}

###### OLDER VERSION: "mean contrubution" included
for(i in seq(1,5)) {
#for(i in seq(nrow(tf.topic50.case)-4,nrow(tf.topic50.case))) {
    r = tf.topic50.case[i,]
  
  
  tf = as.character(r[1])
  cell1 = as.character(r[2])
  cell2 = as.character(r[3])
  
  cell.tf1 = paste(cell1, tf, sep=":")
  cell.tf2 = paste(cell2, tf, sep=":")
  
  data.cell1 = data[cell.tf1,]
  data.cell2 = data[cell.tf2,]
  cmp = data.cell1 - data.cell2
  c1.gain = names(cmp)[cmp==1]
  c2.gain = names(cmp)[cmp==-1]
  
  c1.gain.mat=lda.t50.t2t[c1.gain,]
  c2.gain.mat=lda.t50.t2t[c2.gain,]
  
  ### Record gain/loss targets
  tf.list = list(cell1=cell1, cell2=cell2, c1.gain=c1.gain, c2.gain=c2.gain, 
                 c1.gain.mat=c1.gain.mat,
                 c2.gain.mat=c2.gain.mat)
  id = paste(paste(tf, cell1, sep="_"), cell2, sep="_")
  gainloss.gene.list[[id]] = tf.list
  
  topic50.cell.tf1 = topics.50[cell.tf1, ]
  topic50.cell.tf2 = topics.50[cell.tf2, ]
  
  c1.gain.mat.mean = colMeans(c1.gain.mat)
  c2.gain.mat.mean = colMeans(c2.gain.mat)
  
  ### Circular
  ymax = max(c(c1.gain.mat.mean, c2.gain.mat.mean))
  #ymax = signif(max(c(c1.gain.mat.mean, c2.gain.mat.mean)), digits=1)
  
  dat.break = data.frame(id=seq(51,55), mean=rep(NA,5), prop=rep(NA,5))
  
  dat.plot1 = data.frame(id=seq(1,50), mean=c1.gain.mat.mean, prop=as.numeric(topic50.cell.tf1)*-0.85*ymax-ymax*0.15)
  dat.plot2 = data.frame(id=seq(1,50), mean=c2.gain.mat.mean, prop=as.numeric(topic50.cell.tf2)*-0.85*ymax-ymax*0.15)
  
  dat.plot1 = rbind(dat.plot1, dat.break)
  dat.plot2 = rbind(dat.plot2, dat.break)
  
  label_data = data.frame(id = dat.plot1$id, angle= 90 - 360 * (dat.plot1$id-0.5) /length(dat.plot1$id), label=dat.plot1$id)
  label_data$label[51:55] = NA
  
  p = ggplot(dat.plot1, aes(x=id, y=mean)) +       
    geom_bar(data=dat.plot1, aes(x=id, y=mean), stat="identity", fill=alpha("blue", 0.5)) +
    geom_bar(data=dat.plot2, aes(x=id, y=mean), stat="identity", fill=alpha("red", 0.5)) +
    geom_line(data=dat.plot1, aes(x=id, y=prop), color=alpha("blue", 0.5)) +
    geom_line(data=dat.plot2, aes(x=id, y=prop), color=alpha("red", 0.5)) +
    scale_y_continuous(limits=c(-ymax, ymax), breaks=seq(0, signif(ymax, digits=1), length.out=4)) + 
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5), 
      plot.margin = unit(rep(1,4), "cm")) +
    coord_polar(start = 0) +
    ggtitle(tf)
  
  p = p + geom_text(data=data.frame(x = rep(0, 4), y = seq(0, signif(ymax, digits=1), length.out=4)), aes(x,y), 
                    label = signif(seq(0, signif(ymax, digits=1), length.out=4), digits=1), 
                    hjust=1.125, size=4)
  
  p = p + geom_text(data=label_data, aes(x=id, y=-ymax*0.04, label=label), 
                    size=4, angle= label_data$angle, 
                    inherit.aes = FALSE, hjust = 0.75 )
  
  
  ggsave(filename=paste("case.gm-k562/topic.distr/plot_",i, "_", id, "_barplot.pdf", sep=""), plot=p, width=8, height=8)
 
}
