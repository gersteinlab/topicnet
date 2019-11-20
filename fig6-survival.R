setwd('/Users/tianxiao/Documents/Lab/TF membership/fig5')

load("TFmem_data.fig5-20190208.RData")

library(reshape2)
library(entropy)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(preprocessCore)
library(survival)
library(survminer)
library(fgsea)
library(MASS)

### Set up cancer type
ct = 'BRCA'

### Load data
tf.leu.top = tf.topic50.cell.klsymm.df.leu.top$TF

dat.exp.raw = read.table("TCGA/"+ct+"/expression/expression_matrix.txt", sep = "\t",
                         check.names=FALSE, row.names = 1, header = TRUE)
dat.exp.raw = dat.exp.raw[rowSums(dat.exp.raw[-1])!=0,]
                           
dat.exp.norm = as.data.frame(normalize.quantiles(as.matrix(dat.exp.raw[-1]), copy=TRUE))
rownames(dat.exp.norm) = rownames(dat.exp.raw)
colnames(dat.exp.norm) = colnames(dat.exp.raw)[-1]
dat.exp.norm = cbind(dat.exp.raw$Gene, dat.exp.norm)
colnames(dat.exp.norm)[1] = "Gene"
dat.exp.norm$Gene = as.character(dat.exp.norm$Gene)
dat.exp.norm = dat.exp.norm[,!duplicated(colnames(dat.exp.norm))]
  
dat.survival = read.table("TCGA/"+ct+"/survival/survival.txt", sep = "\t", row.names = 2, header = TRUE)
dat.survival = dat.survival[row.names(dat.survival) %in% colnames(dat.exp.norm),]
dat.survival = dat.survival[!is.na(dat.survival$vital_status),]

sample.id = rownames(dat.survival)
dat.survival = data.frame(days=dat.survival$days_to_death, 
                          days_to_last_follow_up=dat.survival$days_to_last_follow_up, 
                          vital_status=dat.survival$vital_status)
rownames(dat.survival) = sample.id

for(i in seq(1, nrow(dat.survival))) {
  if(dat.survival$vital_status[i]=="alive"){dat.survival$days[i] = dat.survival$days_to_last_follow_up[i]}
}

dat.survival = dat.survival[-2]
dat.survival = dat.survival[!is.na(dat.survival$days),]
dat.survival$vital_status = ifelse(dat.survival$vital_status == 'alive', 0,1)

###### Topic weight * gene expression
dat.exp.sub = dat.exp.norm[dat.exp.norm$Gene %in% rownames(lda.t50.t2t),]

gene.dup = names(which(table(dat.exp.sub$Gene)>1))

for(gd in gene.dup){
  dat.gd = dat.exp.norm[dat.exp.norm$Gene==gd,-1]
  dat.exp.sub = dat.exp.sub[dat.exp.sub$Gene!=gd,]
  dat.gd.mean = colMeans(dat.gd)
  dat.exp.sub = rbind(dat.exp.sub, c(0,dat.gd.mean))
  #dat.exp.sub[nrow(dat.exp.sub),-1] = as.numeric(dat.exp.sub[nrow(dat.exp.sub),-1])
  dat.exp.sub[nrow(dat.exp.sub),1] = gd
}

rownames(dat.exp.sub) = dat.exp.sub$Gene
dat.exp.sub = dat.exp.sub[-1]
lda.t50.t2t.sub = lda.t50.t2t[rownames(dat.exp.sub),]

dat.exp.sub.rank = t(apply(dat.exp.sub, 1, function(x){v=order(x);return(v/length(x))}))
colnames(dat.exp.sub.rank) = colnames(dat.exp.sub)

t50.val = t(t(lda.t50.t2t.sub) %*% as.matrix(dat.exp.sub.rank))
t50.val = t50.val[rownames(dat.survival),]
colnames(t50.val) = paste("X", seq(1,50), sep="")
dat.leu = data.frame(vital_status=dat.survival$vital_status,
                     days=dat.survival$days, t50.val)

#status = c(1,2)
#names(status) = c("alive", "dead")
#dat.leu$status = status[dat.leu$vital_status]

res.cox = coxph(as.formula(paste('Surv(days, vital_status) ~ ',
                                 paste(colnames(t50.val), collapse="+"))), data=dat.leu)

cox.coef = summary(res.cox)$coefficients
cox.topic = rownames(cox.coef)[cox.coef[,5]<0.05]

for( t in cox.topic ){
  #res.cox2 = coxph(Surv(days_to_death, status) ~ X21, data=dat.leu)
  dat.t = dat.leu
  km = kmeans(dat.t[,t], 3)
  
  dat.t$label = 0
  dat.t$label[(dat.t[,t]>quantile(dat.t[,t], 0.5))] = 1
  dat.t$label[(dat.t[,t]<quantile(dat.t[,t], 0.5))] = 2
  dat.t = dat.t[dat.t$label>0,]
  
  #dat.t$label = (dat.t[,t]>median(dat.t[,t]))+1
  #dat.t$label = (dat.t[,t]>25)+1
  #dat.t$label = (dat.t[,t]>quantile(dat.t[,t], 0.5))+1
  val.t = dat.t[,t]
  mean.label = c(mean(val.t[dat.t$label==1]), mean(val.t[dat.t$label==2]))
  dat.t$cluster = round(mean.label[dat.t$label],2)
  
  #dat.t$cluster = round(km$centers[km$cluster],2)
  
  fit = survfit(Surv(days, vital_status) ~ cluster, data=dat.t)
  #pdf(paste("TCGA/figures/survival/survival-",t,".pdf", sep=""), width=10, height=10, onefile = FALSE)
  gg=ggsurvplot(fit, data = dat.t, pval = TRUE, , conf.int = TRUE,
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                legend=c(0.1,0.1))
  gg
  ggsave(paste("TCGA/"+ct+"/figures/survival-",t,".pdf", sep=""), print(gg), width=10, height=10)
  #dev.off()
}



### Component analysis of significant topics
gene.set.kegg <- gmtPathways("gene_set/c2.cp.kegg.v6.2.symbols.gmt")
gene.set.cancer <- gmtPathways("gene_set/c4.cm.v6.2.symbols.gmt")
gene.set.onco <- gmtPathways("gene_set/c6.all.v6.2.symbols.gmt")

plot.gsea <- function(gsea, pathways, ranks, path=NULL, n=10){
  if(length(ranks)==0){return(0)}
  if(!is.null(path)){pdf(path, width = 12, height = 10)}
  topPathwaysUp <- gsea[ES > 0][head(order(pval), n=n), pathway]
  #topPathwaysDown <- gsea[ES < 0][head(order(pval), n=n), pathway]
  #topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  topPathways <- topPathwaysUp
  plotGseaTable(pathways[topPathways], ranks, gsea, 
                gseaParam = 0.5)
}


for( i in seq(1,50)){
  pdf(paste("TCGA/topic.genes/gsea/oncosig/oncosig_",i,".pdf", sep=""), height = 10, width=6)
  gsea.cancer = fgsea(gene.set.onco, lda.t50.t2t[,i], nperm=1000)
  plot.gsea(gsea.cancer, gene.set.onco, lda.t50.t2t[,i])
  dev.off()
}

for( i in seq(1,50)){
  pdf(paste("TCGA/topic.genes/gsea/cancermodule/cm_",i,".pdf", sep=""), height = 10, width=6)
  gsea.cancer = fgsea(gene.set.cancer, lda.t50.t2t[,i], nperm=1000)
  plot.gsea(gsea.cancer, gene.set.cancer, lda.t50.t2t[,i])
  dev.off()
}


### TFs related to significant topics
topics.50.k562 = topics.50[grep("K562", rownames(topics.50)),]
topics.50.gm = topics.50[grep("GM12878", rownames(topics.50)),]

for(i in seq(1,50)){
  topics.50.k562.t21 = topics.50.k562[i]
  topics.50.k562.top = topics.50.k562[order(topics.50.k562.t21, decreasing = TRUE)[1:20],]
  topics.50.k562.t21.top = topics.50.k562.top[,i]
  names(topics.50.k562.t21.top) = rownames(topics.50.k562.top)

  pdf(paste("TCGA/TF/TF_enrichment-",i,".pdf", sep=""), width=10, height=10, onefile = FALSE)
  par(mai=c(2,2,2,2))
  barplot(t(rev(topics.50.k562.t21.top)), horiz = TRUE, las=2, main=colnames(topics.50)[i])
  dev.off()
}


