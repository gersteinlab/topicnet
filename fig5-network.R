# Plot gene networks from stringdb

target.0m = colnames(chip.mcf7)[which(chip.mcf7[1,]==1)]
target.10m = colnames(chip.mcf7)[which(chip.mcf7[4,]==1)]
target.40m = colnames(chip.mcf7)[which(chip.mcf7[5,]==1)]
target.160m = colnames(chip.mcf7)[which(chip.mcf7[6,]==1)]

t3.weight = lda.t50.t2t[,3]
t3.top.genes = names(sort(t3.weight, decreasing=TRUE)[1:500])
t3.top.genes.target.0min = t3.top.genes[t3.top.genes %in% target.0m]

t4.weight = lda.t50.t2t[,4]
t4.top.genes = names(sort(t4.weight, decreasing=TRUE)[1:500])
t4.top.genes.target.10min = t4.top.genes[t4.top.genes %in% target.10m]
t4.top.genes.target.40min = t4.top.genes[t4.top.genes %in% target.40m]
t4.top.genes.target.160min = t4.top.genes[t4.top.genes %in% target.160m]


### Plot network
library("igraph")
library("plotrix")

## Function for plotting an elliptical node
myellipse <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/30 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  draw.ellipse(x=coords[,1], y=coords[,2],
               a = vertex.size, b=vertex.size/2, col=vertex.color, border="deepskyblue4")
}

# Register the shape with igraph
add_shape("ellipse", clip=shapes("circle")$clip, plot=myellipse)


# Load gene list
imp.gene.symbol = t3.top.genes.target.0min
string.network = read.table("plot.network/t3_0min.tsv", sep="\t")

string.network.edges = cbind(as.character(string.network$V1), as.character(string.network$V2))

string.network.nodes = data.frame(gene=unique(as.character(c(string.network.edges[,1], string.network.edges[,2]))))
string.network.nodes$label = rep(0, length(string.network.nodes))

string.network.nodes$label[string.network.nodes$gene %in% imp.gene.symbol] = 1
string.network.nodes$gene = as.character(string.network.nodes$gene)

graph = graph.data.frame(string.network.edges, vertices = string.network.nodes, directed = F)
V(graph)$color = ifelse(V(graph)$label==1,"deepskyblue3","lightblue1")
V(graph)$label.cex = 0.5

pdf("plot.network/t3_0min.pdf", width=10, height=10)
plot(graph, vertex.label=V(graph)$name, vertex.shape="ellipse", vertex.size=2)
dev.off()
