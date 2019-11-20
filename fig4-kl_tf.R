# Calculate KL divergence for different cells given a TF

load("alltf.rewiring.rdata")
load("TFmem_data.fig4-190116.RData")
load("med.fgsea.rdata")

library("reshape2")
library(entropy)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(gplots)
library(ggplot2)

tf.count = table(tf.doc$tf)
tf.filtered = names(tf.count)[tf.count>=2]

######## Rewiring of all TFs (Cell 1: Cell 2: TF: KL divergence)
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

tf.topic50.list = list() # tf: topic50 values in relative cell-tf documents
tf.topic50.cell.kl.list = list() # tf: KL divergence based on topic50 values between cell types
tf.topic50.cell.klsymm.list = list() # tf: symmetric KL divergence matrix (mean of two directions)

for(tf in tf.filtered){
  tf.doc.id = which(tf.doc$tf == tf)
  n.doc = length(tf.doc.id)
  tf.topic50 = topics.50[tf.doc.id,]
  tf.cell = tf.doc[tf.doc.id,]$cell
  rownames(tf.topic50) = tf.cell
  
  if(nrow(tf.topic50)<2) next # only calculate for TFs with >=2 cell line data
  tf.topic50.list[[tf]] = tf.topic50
  
  tf.topic50.cell.kl = pair.mat(tf.topic50, KL.plugin)
  tf.topic50.cell.kl.list[[tf]] = tf.topic50.cell.kl
  
  tf.topic50.cell.klsymm = (tf.topic50.cell.kl + t(tf.topic50.cell.kl))/2
  tf.topic50.cell.klsymm.list[[tf]] = tf.topic50.cell.klsymm
}

### Plot the KL divergence for each tf
for(tf in names(tf.topic50.cell.klsymm.list)) {
  tf.topic50.cell.klsymm = tf.topic50.cell.klsymm.list[[tf]]
  pdf(paste("/Users/tianxiao/Documents/Lab/TF membership/fig4/tf.klsymm/", tf, "-klsymm.pdf", sep=""), height = 10, width = 10)
  heatmap.2(tf.topic50.cell.klsymm, scale="none", symm=TRUE, trace="none", col = colorRampPalette(c("red","yellow"))(32), main=tf,
            distfun = function(x) as.dist(x)) # cluster by the sqrt distance matrix
  dev.off()
}

########## Mean of KL for each TF
tf.topic50.cell.klsymm.mean = unlist(lapply(tf.topic50.cell.klsymm.list, function(m) {n = nrow(m); sum(m)/(n^2-n)}))
tf.topic50.cell.klsymm.mean.bot.tf = sort(tf.topic50.cell.klsymm.mean)[1:10]
tf.topic50.cell.klsymm.mean.top.tf = sort(tf.topic50.cell.klsymm.mean, decreasing = TRUE)[1:10]

tf.topic50.cell.klsymm.mean.tf = c(tf.topic50.cell.klsymm.mean.top.list, tf.topic50.cell.klsymm.mean.bot.list)



########## Top/botton rewiring pairs (TF: cell1 vs cell2)
tf.topic50.cell.klsymm.df.list = lapply(tf.topic50.cell.klsymm.list, function(x){
  xx=x;xx[lower.tri(xx)]=0;xm=melt(xx);colnames(xm)=c("Cell-1", "Cell-2","Rewire");return(xm)
})

tf.topic50.cell.klsymm.raw = do.call("rbind", tf.topic50.cell.klsymm.df.list) # Data frame (TF, Cell1, Cell2, rewire)
tf.topic50.cell.klsymm.raw$TF = unlist(lapply(strsplit(rownames(tf.topic50.cell.klsymm.raw), "[.]"), function(x){x[1]}))
tf.topic50.cell.klsymm.raw = tf.topic50.cell.klsymm.raw[,c("TF","Cell-1", "Cell-2","Rewire")]

tf.topic50.cell.klsymm.raw$`Cell-1` = as.character(tf.topic50.cell.klsymm.raw$`Cell-1`)
tf.topic50.cell.klsymm.raw$`Cell-2` = as.character(tf.topic50.cell.klsymm.raw$`Cell-2`)

tf.topic50.cell.klsymm.df = tf.topic50.cell.klsymm.raw[tf.topic50.cell.klsymm.raw$Rewire!=0,]
tf.topic50.cell.klsymm.df = tf.topic50.cell.klsymm.df[order(tf.topic50.cell.klsymm.df$Rewire, decreasing=TRUE),]

tf.topic50.cell.klsymm.df.top.pair = tf.topic50.cell.klsymm.df[1:10,]
tf.topic50.cell.klsymm.df.bot.pair = tf.topic50.cell.klsymm.df[(nrow(tf.topic50.cell.klsymm.df)-9):nrow(tf.topic50.cell.klsymm.df),]

tf.topic50.cell.klsymm.df.pair = rbind(tf.topic50.cell.klsymm.df.top.pair, tf.topic50.cell.klsymm.df.bot.pair)


########## Plot histogram of mean rewiring for TFs
tf.topic50.cell.klsymm.mean.df = data.frame(x=names(tf.topic50.cell.klsymm.mean),
                                            y=tf.topic50.cell.klsymm.mean)


pdf("rewire/average_rewiring.pdf", width = 8, height = 8)
angles= c(360/53*seq(53,27)+90, 360/53*seq(52,27)+90)
p = ggplot(tf.topic50.cell.klsymm.mean.df, aes(x=x, y=y, fill=y)) +
  geom_histogram(binwidth=1, stat='identity') +theme_light() +
  scale_fill_gradient2(low='blue', mid="cyan",high='red',midpoint = 0.5) +
  theme(axis.title.y=element_text(angle=0)) + 
  theme(axis.text.x=element_text(angle=angles, vjust = 0, hjust=0))+
  theme(panel.border = element_blank())+
  coord_polar() + theme(axis.text = element_text(size = 10, angle=0))+ aes(x=reorder(x, y))+
  xlab("TF")+ylab("Cell specifity")
plot(p)
dev.off()

########## Plot table for top rewiring events
tf.topic50.cell.klsymm.df.top20.pair = tf.topic50.cell.klsymm.df[1:20,]
tf.topic50.cell.klsymm.df.bot20.pair = tf.topic50.cell.klsymm.df[(nrow(tf.topic50.cell.klsymm.df)-19):nrow(tf.topic50.cell.klsymm.df),]

# Color panels
tf.all = unique(c(tf.topic50.cell.klsymm.df.top20.pair$TF, tf.topic50.cell.klsymm.df.bot20.pair$TF))
cell.all = unique(c(tf.topic50.cell.klsymm.df.top20.pair$`Cell-1`, tf.topic50.cell.klsymm.df.bot20.pair$`Cell-1`,
                    tf.topic50.cell.klsymm.df.top20.pair$`Cell-2`, tf.topic50.cell.klsymm.df.bot20.pair$`Cell-2`))

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))

col.tf = col_vector[1:length(tf.all)]
names(col.tf) = tf.all

col.cell = col_vector[1:length(cell.all)]
names(col.cell) = cell.all

top.fill = cbind(col.tf[tf.topic50.cell.klsymm.df.top20.pair$TF], 
                 col.cell[tf.topic50.cell.klsymm.df.top20.pair$`Cell-1`], 
                 col.cell[tf.topic50.cell.klsymm.df.top20.pair$`Cell-2`])
#bot.fill = cbind(col.tf[tf.rew.bot$TF], col.cell[tf.rew.bot$`Cell-1`], col.cell[tf.rew.bot$`Cell-2`])

plotGrid <- function(dat, dat.fill) {
  tt1 = ttheme_default(base_size=8, core=list(bg_params = list(fill=dat.fill[,1])))
  p1 = tableGrob(dat$TF, cols = "TF", rows = NULL, theme = tt1)
  p1$widths = unit(1.5, "cm")
  
  tt2 = ttheme_default(base_size=8, core=list(bg_params = list(fill=dat.fill[,2])))
  p2 = tableGrob(dat$`Cell-1`, cols = "Cell-1", rows = NULL, theme = tt2)
  p2$widths = unit(4.5, "cm")
  
  tt3 = ttheme_default(base_size=8, core=list(bg_params = list(fill=dat.fill[,3])))
  p3 = tableGrob(dat$`Cell-2`, cols = "Cell-2", rows = NULL, theme = tt3)
  p3$widths = unit(4.5, "cm")
  
  tt4 = ttheme_default(base_size=8)
  p4 = tableGrob(round(dat$Rewire,4), cols = "Rewire", rows = NULL, theme = tt4)
  p4$widths = unit(2, "cm")
  
  grid.arrange(gtable_combine(p1,p2,p3,p4, along=1))
}

pdf("rewire/KL_table.top.pdf", width = 8, height = 12)
plotGrid(tf.topic50.cell.klsymm.df.top20.pair, top.fill)
dev.off()


