# Rewiring case study: GM12878 vs K562
# Gene expression change of the top-weighted genes in rewired topics

load("case_TFmem_data.fig4-190904.RData")

library(reshape2)
library(entropy)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(preprocessCore)

# Set up tfs and cells to study
tf = 'ZBTB33'
#tf = 'EP300'
cell1='GM12878'
cell2='K562'

doc.row1 = paste(cell1, tf, sep=":")
doc.row2 = paste(cell2, tf, sep=":")
topics1 = topics.50[doc.row1,]
topics2 = topics.50[doc.row2,]

rew.diff = topics1 - topics2

exp.raw = read.table("data/expression/activeGene.4chengfei.txt", header=TRUE)
exp.clean = exp.raw[complete.cases(exp.raw),c("tname","expr","cell")]

mat.exp = dcast(exp.clean, tname~cell, mean, value.var = "expr", drop=TRUE)
mat.exp = mat.exp[complete.cases(mat.exp),]
rownames(mat.exp) = mat.exp$tname
mat.exp = mat.exp[-1]
mat.exp = as.matrix(mat.exp)

mat.exp.qm = normalize.quantiles(mat.exp,copy=TRUE)
rownames(mat.exp.qm) = rownames(mat.exp)
colnames(mat.exp.qm) = colnames(mat.exp)

mat.exp.log = log(mat.exp.qm,2)

#exp.k562 = exp.clean[exp.clean$cell=="k562",c("tname","expr")]
#exp.gm = exp.clean[exp.clean$cell=="gm",c("tname","expr")]

rew.topic = list()
rew.topic[["ZBTB33"]] = c(49,16,31,17)
rew.topic[["EP300"]] = c(34,10,36,5)


name = paste(tf,cell1,cell2,sep="-")
rew.topic.tf = rew.topic[[tf]]

gene.names = colnames(data)

for(i in seq(1,length(rew.topic.tf))) {
  t = rew.topic.tf[i]
  dat1 = data[paste(cell1, tf, sep=":"),]
  dat2 = data[paste(cell2, tf, sep=":"),]
  
  dat.diff = dat1 - dat2
  
  
  if(rew.diff[,t]>0){ 
    # choose the cell line to compare expression 
    # check cell line-specific genes in the cell line enriched with the selected topic
    cell.gain = cell1
    gain.gene = gene.names[dat.diff==1]
  }else{
    cell.gain = cell2
    gain.gene = gene.names[dat.diff==-1]
  }
  
  top.genes.t = get_top_genes(lda.t50.t2t, t, n=100)
  exp.top.gene = mat.exp.log[rownames(mat.exp.log) %in% top.genes.t,]
  exp1 = exp.top.gene[,"gm"]
  exp2 = exp.top.gene[,"k562"]
  
  pdf(paste("case.gm-k562/top.genes.expr/", name, "_rew_", i, "_topic_", t,"_top100_exp.pdf", sep=""),
      width=6, height=6)
  boxplot(exp1, exp2, names = c(cell1, cell2), col=col.list[c(cell1, cell2)], 
          main=paste(tf,"-Topic ",t,sep=""),
          sub=paste("Enriched in", cell.gain))
  dev.off()
  
  
  exp1.gain.gene = exp1[names(exp1) %in% gain.gene]
  exp2.gain.gene = exp2[names(exp2) %in% gain.gene]
  
  pdf(paste("case.gm-k562/top.genes.expr/", name, "_rew_", i, "_topic_", t,"_top100_", cell.gain, "-gain_exp.pdf", sep=""),
      width=6, height=6)
  boxplot(exp1.gain.gene, exp2.gain.gene, names = c(cell1, cell2), col=col.list[c(cell1, cell2)], 
          main=paste(tf,"-Topic ",t,sep=""),
          sub=paste("Enriched in", cell.gain))
  dev.off()
}


####### Cell cycle-related genes
cc.gene.list = c("UHRF2", "MAU2", "NCAPG", "USP39", "SPDL1", "CENPT", "CDC20", "CDC5L", "SEPT7", "GSG2", "SMC4")

exp1.cc = mat.exp.log[rownames(mat.exp) %in% cc.gene.list,"gm"]
exp2.cc = mat.exp.log[rownames(mat.exp) %in% cc.gene.list,"k562"]

dat1.cc = data["GM12878:ZBTB33", colnames(data) %in% cc.gene.list]
dat2.cc = data["K562:ZBTB33", colnames(data) %in% cc.gene.list]

dat.diff.cc = dat1.cc - dat2.cc
cc.gm.gain.list = names(dat.diff.cc)[which(dat.diff.cc==1)]

exp1.gain.cc = mat.exp.log[rownames(mat.exp) %in% cc.gm.gain.list,"gm"]
exp2.gain.cc = mat.exp.log[rownames(mat.exp) %in% cc.gm.gain.list,"k562"]

pdf(paste("case.gm-k562/top.genes.expr/ZBTB33-GM12878-K562_rew_1_topic_49_top100_GM12878-gain_cell-cycle_exp.pdf", sep=""),
    width=6, height=6)
boxplot(exp1.gain.cc, exp2.gain.cc, names = c(cell1, cell2), col=col.list[c(cell1, cell2)], 
        main="ZBTB33-Topic 49 (cell cycle-related)",
        sub=paste("Enriched in GM12878"))
dev.off()



