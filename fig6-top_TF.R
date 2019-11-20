# Top TFs in the topics

load("TFmem_data.fig4-190227.RData")

cell.count = table(tf.doc$cell)
cell.filtered = names(cell.count)[cell.count>50]

for(cell in cell.filtered){
  topic50.cell = topics.50[grep(cell,rownames(topics.50)),]
  tf.list.cell = rownames(topic50.cell)
  for(t in colnames(topic50.cell)){
    topic50.v = topic50.cell[,t]
    topic50.v.top10 = topic50.v[order(topic50.v, decreasing=TRUE)[1:10]]
    names(topic50.v.top10) = tf.list.cell[order(topic50.v, decreasing=TRUE)[1:10]]
    topic50.v.top10 = sort(topic50.v.top10)
    pdf(paste("TF.enrichment/", cell,"/", cell,"-",t,".pdf", sep=""))
    par(mai=c(2,2,2,2))
    barplot(topic50.v.top10, horiz = T, las=2)
    dev.off()
  }
}

tf.list = rownames(topics.50)
for(t in colnames(topics.50)){
  topic50.v = topics.50[,t]
  topic50.v.top10 = topic50.v[order(topic50.v, decreasing=TRUE)[1:10]]
  names(topic50.v.top10) = tf.list[order(topic50.v, decreasing=TRUE)[1:10]]
  topic50.v.top10 = sort(topic50.v.top10)
  pdf(paste("TF.enrichment/all/ALL-",t,".pdf", sep=""))
  par(mai=c(2,2,2,2))
  barplot(topic50.v.top10, horiz = T, las=2)
  dev.off()
}
