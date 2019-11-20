# Calculate KL divergence for TFs in a given cell line

load("alltf.rewiring.rdata")
load("TFmem_data.fig4-190115.RData")

library(reshape2)
library(entropy)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(gplots)

cell.count = table(tf.doc$cell)
cell.filtered = names(cell.count)[cell.count>50]

cell.topic50.tf.list = list() # tf: topic50 values in relative cell-tf documents
cell.topic50.tf.kl.list = list() # tf: KL divergence based on topic50 values between cell types
cell.topic50.tf.klsymm.list = list() # tf: symmetric KL divergence matrix (mean of two directions)

for(cell in cell.filtered){
  cell.doc.id = which(tf.doc$cell == cell)
  n.doc = length(cell.doc.id)
  cell.topic50 = topics.50[cell.doc.id,]
  cell.tf = tf.doc[cell.doc.id,]$tf
  rownames(cell.topic50) = cell.tf
  
  cell.topic50.tf.list[[cell]] = cell.topic50
  
  cell.topic50.tf.kl = pair.mat(cell.topic50, KL.plugin)
  cell.topic50.tf.kl.list[[cell]] = cell.topic50.tf.kl
  
  cell.topic50.tf.klsymm = (cell.topic50.tf.kl + t(cell.topic50.tf.kl))/2
  cell.topic50.tf.klsymm.list[[cell]] = cell.topic50.tf.klsymm
}

### Plot the KL divergence for each cell
for(cell in names(cell.topic50.tf.klsymm.list)) {
  cell.topic50.tf.leu.klsymm = cell.topic50.tf.klsymm.list[[cell]]
  pdf(paste("cell.klsymm/", cell, "-klsymm.pdf", sep=""), height = 10, width = 10)
  heatmap.2(cell.topic50.tf.leu.klsymm, scale="none", symm=TRUE, trace="none", col = colorRampPalette(c("red","yellow"))(32), main=cell,
            distfun = function(x) as.dist(x)) # cluster by the sqrt distance matrix
  dev.off()
} 
