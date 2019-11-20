# Summarize rewiring values of TFs across cell lines

load("alltf.rewiring.rdata")
load("TFmem_data.fig4-190116.RData")
load("med.fgsea.rdata")

library(reshape2)
library(entropy)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(gplots)
library(ggplot2)

### Use GM12878 as ref (disabled)
#cell.ref = "GM12878"

#cell.count = table(tf.doc$cell)
#cell.plot = c("HeLa-S3","HepG2","K562")
#cell.filtered = cell.plot

#cell.plot = names(which(cell.count>50))
#cell.filtered = cell.plot
#cell.plot = cell.plot[-(cell.plot==cell.ref)]
#tf.doc.ref = tf.doc[tf.doc$tf %in% tf.doc[tf.doc$cell == cell.ref,]$tf,]

### Pairwise comparison
cell.plot = c("HeLa-S3","HepG2","K562", "MCF-7", "GM12878")

tf.doc.ref = tf.doc[tf.doc$cell %in% cell.plot,]

tf.shared.count = table(tf.doc.ref[tf.doc.ref$cell %in% cell.filtered,]$tf)
#tf.shared = names(tf.shared.count)[tf.shared.count==4]
tf.shared = names(tf.shared.count)[tf.shared.count>=3]

comp.list = combn(cell.plot, 2)

rew.mat = matrix(0L, nrow=ncol(comp.list), ncol=length(tf.shared))
rownames(rew.mat) = paste(comp.list[1,],comp.list[2,], sep=" vs ")
colnames(rew.mat) = tf.shared

for(i in seq(1,ncol(comp.list))){
  cell = comp.list[1,i]
  cell.ref = comp.list[2,i]
  for(tf in tf.shared){
    plot.id = paste(cell, tf, sep=":")
    ref.id = paste(cell.ref, tf, sep=":")
    
    if(!((plot.id %in% rownames(topics.50)) && (ref.id %in% rownames(topics.50)))){
      rew.mat[i, tf] = NA
      next
    }
    
    topic50.plot = topics.50[plot.id,]
    topic50.ref = topics.50[ref.id,]
    
    rew.mat[i, tf] = (KL.plugin(topic50.plot, topic50.ref) + KL.plugin(topic50.ref, topic50.plot))/2
  }
}

rew.mat.raw = rew.mat

#rew.mat = rew.mat.raw[,which(apply(rew.mat.raw, 2, function(x){sum(is.na(x))})<nrow(rew.mat.raw)-2)]
#rew.mat = rew.mat.raw[,which(apply(rew.mat.raw, 2, function(x){sum(is.na(x))})==0)]

rew.mat.mean = apply(rew.mat, 2, function(x){mean(x, na.rm=TRUE)})
rew.mat = rbind(rew.mat, rew.mat.mean)
rownames(rew.mat)[nrow(rew.mat)] = "mean"

rew.order = order(rew.mat.mean)
rew.mat = rew.mat[,rew.order]

pdf("tf.klsymm/summary/top.rew.full.long.pdf", height=5,width=20)
heatmap.2(rew.mat, col = colorRampPalette(c("blue3","white","red", "red2","red3","red4"))(32), 
          trace = "none", cexCol = 0.8, cexRow = 1,
          Rowv = FALSE, Colv = FALSE, rowsep = nrow(rew.mat)-1, na.color = "grey")
dev.off()

