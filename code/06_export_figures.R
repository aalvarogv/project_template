library(igraph)
library(ggplot2)

g <- readRDS("../../RESULTS/processed/network_graph.rds")

pdf("../../RESULTS/figures/network_plot.pdf")
plot(g, vertex.size=5, vertex.label=NA)
dev.off()
