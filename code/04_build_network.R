library(igraph)
library(tidyverse)

edges <- read.csv("../../RESULTS/raw/string_edges.csv")

g <- graph_from_data_frame(edges, directed = FALSE)

# eliminar genes sin conexiones
g <- delete_vertices(g, V(g)[degree(g) == 0])

# métricas
metrics <- data.frame(
  gene = V(g)$name,
  degree = degree(g),
  betweenness = betweenness(g),
  closeness = closeness(g)
)

write.csv(metrics, "../../RESULTS/processed/network_metrics.csv", row.names = FALSE)

# guardar objeto de red
saveRDS(g, "../../RESULTS/processed/network_graph.rds")
