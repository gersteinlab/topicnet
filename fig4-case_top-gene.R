# Rewiring case study: GM12878 vs K562
# Top genes of rewired topics

load("case_TFmem_data.fig4-190904.RData")


library(reshape2)
library(entropy)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(gplots)
library(ggplot2)

t = 3

# Set up tfs and cells to study
tf="ZBTB33"
cell1 = "GM12878"
cell2 = "K562"

dat1 = data[paste(cell1, tf, sep=":"),]
dat2 = data[paste(cell2, tf, sep=":"),]

dat.diff = dat1 - dat2

gene.names = colnames(data)
  
if(rew.diff[,t]>0){ 
  # choose the cell line to compare expression 
  # check cell line-specific genes in the cell line enriched with the selected topic
  cell.gain = cell1
  gain.gene = gene.names[dat.diff==1]
}else{
  cell.gain = cell2
  gain.gene = gene.names[dat.diff==-1]
}

top.genes.t = get_top_genes(lda.t50.t2t, t, n=20)


### Write to csv
top.genes.t.write = top.genes.t[top.genes.t %in% colnames(data)]
dat1.gene = data[paste(cell1, tf, sep=":"),top.genes.t.write]
dat2.gene = data[paste(cell2, tf, sep=":"),top.genes.t.write]

data.write = data.frame(gene=as.character(top.genes.t.write), 
                        gm=as.numeric(dat1.gene), k=as.numeric(dat2.gene))
  

write.csv(data.write, quote=FALSE,
          paste("top.genes/top.genes-topic",t,"-",tf,".csv", sep=""))


### Draw barplot
top.genes.val = lda.t50.t2t[rev(top.genes.t), t]
xlim = c(min(top.genes.val)*0.75, max(top.genes.val))

pdf(paste("top.genes/top.genes-topic", t, ".pdf", sep=""),
    height=8, width=6)
par(mai=c(2,2,2,2))
bp = barplot(top.genes.val, space=0.2, horiz=TRUE, xlim=xlim, xpd=FALSE, xlab='contribution',
             las=2, main=paste("Topic",t), cex.axis=0.75, col='darkblue')

dev.off()
