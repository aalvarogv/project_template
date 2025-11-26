library(igraph)
library(tidyverse)


base_dir <- dirname(getwd())

edges_path   <- file.path(base_dir, "RESULTS", "raw", "string_edges.csv")
chembl_path  <- file.path(base_dir, "RESULTS", "raw", "chembl_targets.csv")
out_dir      <- file.path(base_dir, "RESULTS", "processed")
metrics_path <- file.path(out_dir, "network_metrics.csv")
graph_path   <- file.path(out_dir, "network_graph.rds")
summary_path <- file.path(out_dir, "network_summary.csv")

cat("Leyendo aristas desde:", edges_path, "\n")

if (!file.exists(edges_path)) {
  stop("No se encuentra el fichero de aristas: ", edges_path)
}

edges <- read.csv(edges_path, stringsAsFactors = FALSE)

# Nos aseguramos de que las columnas gene1 y gene2 existen
if (!all(c("gene1", "gene2") %in% colnames(edges))) {
  stop("El fichero de aristas debe tener columnas 'gene1' y 'gene2'.")
}

# Construir grafo no dirigido a partir de las aristas
g <- graph_from_data_frame(
  d = edges[, c("gene1", "gene2")],
  directed = FALSE
)

cat("Nodos iniciales:", vcount(g), "\n")
cat("Aristas iniciales:", ecount(g), "\n")

# Eliminar genes sin conexiones 
g <- delete_vertices(g, V(g)[degree(g) == 0])

cat("Nodos tras eliminar aislados:", vcount(g), "\n")
cat("Aristas tras eliminar aislados:", ecount(g), "\n")


# Integrar información de ChEMBL 

is_chembl_target <- rep(FALSE, vcount(g))
names(is_chembl_target) <- V(g)$name

if (file.exists(chembl_path)) {
  cat("Integrando información de ChEMBL desde:", chembl_path, "\n")
  chembl <- read.csv(chembl_path, stringsAsFactors = FALSE)
  
  # genes con chembl_id no NA → diana conocida
  chembl_genes <- chembl %>%
    filter(!is.na(chembl_id) & chembl_id != "") %>%
    pull(gene) %>%
    unique()
  
  is_chembl_target[intersect(names(is_chembl_target), chembl_genes)] <- TRUE
}

V(g)$is_chembl_target <- is_chembl_target

# Cálculo de métricas de red
metrics <- tibble(
  gene      = V(g)$name,
  degree    = degree(g),
  betweenness = betweenness(g, normalized = TRUE),
  closeness   = closeness(g, normalized = TRUE),
  clustering_coef = transitivity(g, type = "local", isolates = "zero"),
  is_chembl_target = V(g)$is_chembl_target
)

# Crear directorio de salida si no existe
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

write.csv(metrics, metrics_path, row.names = FALSE)
cat("Métricas de red guardadas en:", metrics_path, "\n")

# Pequeño resumen global de la red
net_summary <- tibble(
  n_nodes   = vcount(g),
  n_edges   = ecount(g),
  density   = edge_density(g),
  transitivity_global = transitivity(g, type = "global"),
  n_chembl_targets = sum(V(g)$is_chembl_target)
)

write.csv(net_summary, summary_path, row.names = FALSE)
cat("Resumen de red guardado en:", summary_path, "\n")

# Guardar objeto de red
saveRDS(g, graph_path)
cat("Grafo guardado en:", graph_path, "\n")

# Detección de comunidades (Algoritmo de Louvain)
# set.seed para reproducibilidad
set.seed(123) 
comm <- cluster_louvain(g)

# Añadir la pertenencia a módulo como atributo del nodo
V(g)$module <- membership(comm)

# Guardar tabla de asignación de módulos
modules_df <- data.frame(
  gene = V(g)$name,
  module = as.numeric(V(g)$module) # ID numérico del cluster
)

modules_path <- file.path(out_dir, "network_modules.csv")
write.csv(modules_df, modules_path, row.names = FALSE)
cat("Módulos guardados en:", modules_path, "\n")

# Actualizar el guardado del grafo para que incluya la info de módulos
saveRDS(g, graph_path)