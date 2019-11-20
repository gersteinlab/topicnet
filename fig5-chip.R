# GSEA of the ChIP-seq data

load("TFmem_data.fig6-20190520.RData")
load("TFmem_data.fig6.RData")
load("med.fgsea.rdata")
load("fig6.groseq.expr.rdata")

library(fgsea)

ESR1.list = paste(c("SRX469393","SRX469394","SRX469395","SRX469396","SRX469397","SRX469398"),"MCF.7:human:ESR1", sep=".")
time = factor(c("0","2","5","10","40","160"), levels=c("0","2","5","10","40","160"))

chip_raw = read.table("ts.mcf7.19147t.csv", sep=",", row.names = 1)
colnames(chip_raw) = rownames(lda.t50.t2t)
# Columns: genes (same as rownames of lda.t50.t2t)
# Rows: TFs (same as colnames of topics.mcf7)

### Overall count of targets
chip.mcf7 = chip_raw[ESR1.list,]
pdf("chip-seq/target_count.pdf", width=8, height=6)
plot(as.numeric(time), rowSums(chip.mcf7), type="b", xlab="time",ylab="Target count", pch=1, xaxt="n", col="blue")
axis(1, at=1:6, labels=time)
dev.off()



### GSEA
gene.set.kegg <- gmtPathways("c2.cp.kegg.v6.2.symbols.gmt")
gene.set.reactome <- gmtPathways("c2.cp.reactome.v6.2.symbols.gmt")
gene.set.bp <- gmtPathways("c5.bp.v6.2.symbols.gmt")


plot.gsea <- function(gsea, pathways, ranks, path=NULL, n=10){
  if(length(ranks)==0){return(0)}
  if(!is.null(path)){pdf(path, width = 12, height = 10)}
  topPathwaysUp <- gsea[ES > 0][head(order(pval), n=n), pathway]
  topPathwaysDown <- gsea[ES < 0][head(order(pval), n=n), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plotGseaTable(pathways[topPathways], ranks, gsea, 
                gseaParam = 0.5)
}

for(i in c(3, 4, 43)){
  gene.imp = lda.t50.t2t[,i]

  gsea = fgsea(gene.set.kegg, gene.imp, nperm=1000)
  plot.gsea(gsea, gene.set.kegg, gene.imp)
}
