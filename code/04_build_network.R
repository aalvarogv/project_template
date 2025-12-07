#!/usr/bin/env Rscript
library(igraph)
library(tidyverse)

# ===== 1. Rutas =====
base_dir <- dirname(getwd())
edges_path   <- file.path(base_dir, "results", "raw", "string_edges.csv")
chembl_path  <- file.path(base_dir, "results", "raw", "chembl_targets.csv")
out_dir      <- file.path(base_dir, "results", "processed")
metrics_path <- file.path(out_dir, "network_metrics.csv")
graph_path   <- file.path(out_dir, "network_graph.rds")
summary_path <- file.path(out_dir, "network_summary.csv")
modules_path <- file.path(out_dir, "network_modules.csv")

cat("Leyendo aristas desde:", edges_path, "\n")
if (!file.exists(edges_path)) {
  stop("No se encuentra el fichero de aristas: ", edges_path)
}

edges <- read.csv(edges_path, stringsAsFactors = FALSE)

# ===== 2. Validaciones =====
req_cols <- c("gene1","gene2","score")
if (!all(req_cols %in% names(edges))) {
  stop("El CSV de aristas debe contener columnas: ",
  paste(req_cols, collapse=", "))
}

# Limpiar datos básicos
edges <- edges %>%
  filter(!is.na(gene1), !is.na(gene2), gene1 != "", gene2 != "") %>%
  filter(gene1 != gene2) %>%                       # sin bucles
  mutate(score = as.numeric(score)) %>%
  filter(!is.na(score), score > 0)                 # score válido

# Peso (fortaleza) y coste (para distancias)
eps <- 1e-6
edges <- edges %>%
  mutate(weight = pmin(pmax(score, eps), 1.0),
         cost   = 1/weight)


# ===== 3. Grafo simple, no dirigido =====
g <- graph_from_data_frame(
  d = edges[, c("gene1","gene2","weight","cost")],
  directed = FALSE
)

# Combinar aristas múltiples conservando el mayor peso y el menor coste
g <- igraph::simplify(
  g,
  remove.multiple = TRUE,
  remove.loops    = TRUE,
  edge.attr.comb  = list(
    weight   = "max",   # para score STRING
    cost     = "min",   # para distancias (1/weight)
    .default = "first"  # resto de atributos (si los hubiera)
  )
)

cat("Nodos iniciales:", vcount(g), "\n")
cat("Aristas iniciales:", ecount(g), "\n")

# Eliminar aislados
g <- delete_vertices(g, V(g)[degree(g) == 0])

cat("Nodos tras eliminar aislados:", vcount(g), "\n")
cat("Aristas tras eliminar aislados:", ecount(g), "\n")


# ===== 4. Integrar ChEMBL =====
V(g)$is_chembl_target <- FALSE
if (file.exists(chembl_path)) {
  cat("Integrando información de ChEMBL desde:", chembl_path, "\n")
  chembl <- read.csv(chembl_path, stringsAsFactors = FALSE)
  if (all(c("gene","chembl_id") %in% names(chembl))) {
    chembl_genes <- chembl %>%
      filter(!is.na(chembl_id), chembl_id != "") %>%
      mutate(gene = toupper(trimws(gene))) %>%
      pull(gene) %>% unique()
    map_names <- toupper(trimws(V(g)$name))
    V(g)$is_chembl_target[map_names %in% chembl_genes] <- TRUE
  }
}


# ===== 5. Métricas (no ponderadas y ponderadas) =====
# Strength (grado ponderado por 'weight')
deg              <- igraph::degree(g, mode = "all")
strg             <- igraph::strength(g, mode = "all", weights = igraph::E(g)$weight)

# Betweenness/closeness no ponderadas
btw_unw          <- igraph::betweenness(g, directed = FALSE, normalized = TRUE)
cls_unw          <- igraph::closeness(g, mode = "all", normalized = TRUE)

# Betweenness/closeness ponderadas usando 'cost' (distancias)
btw_w            <- igraph::betweenness(g, directed = FALSE, normalized = TRUE, weights = E(g)$cost)
cls_w            <- igraph::closeness(g, mode = "all", normalized = TRUE, weights = E(g)$cost)

# Closeness armónica (mejor en redes desconexas)
harmonic_closeness <- harmonic_centrality(g, mode = "all", weights = E(g)$cost, normalized = TRUE)

# Clustering local (ponderado si hay pesos)
clust_local      <- igraph::transitivity(g, type = "local", isolates = "zero", weights = E(g)$weight)

# Atributo de nodo
V(g)$degree            <- deg
V(g)$strength          <- strg
V(g)$betweenness       <- btw_unw
V(g)$betweenness_w     <- btw_w
V(g)$closeness         <- cls_unw
V(g)$closeness_w       <- cls_w
V(g)$closeness_harmonic<- harmonic_closeness
V(g)$clust_coef        <- clust_local


# ===== 6. Comunidades (Louvain, ponderado) =====
set.seed(123)
comm <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
V(g)$module <- membership(comm)


# ===== 7. Resumen global =====
comp <- tryCatch(
  igraph::components(g),
  error = function(e) igraph::clusters(g)  # fallback versiones antiguas / masking
)
giant_id <- which.max(comp$csize)
gc_nodes <- igraph::V(g)[comp$membership == giant_id]
g_giant  <- igraph::induced_subgraph(g, gc_nodes)

# --- helper compatible para longitud media de camino ---
avg_path_len <- function(g, weights = NULL) {
  f <- tryCatch(get("mean_distance", asNamespace("igraph")), error = function(e) NULL)
  if (!is.null(f)) {
    return(f(g, directed = FALSE, weights = weights))
  }
  return(NA_real_)
}

avg_path_unw <- if (igraph::vcount(g_giant) > 1) avg_path_len(g_giant, weights = NULL) else NA_real_
avg_path_w   <- if (igraph::vcount(g_giant) > 1) avg_path_len(g_giant, weights = igraph::E(g_giant)$cost) else NA_real_

assort_deg   <- igraph::assortativity_degree(g, directed = FALSE)
modularity_w <- igraph::modularity(comm, weights = igraph::E(g)$weight)

net_summary <- tibble(
  n_nodes   = igraph::vcount(g),
  n_edges   = igraph::ecount(g),
  density   = igraph::edge_density(g),
  transitivity_global = igraph::transitivity(g, type = "global", weights = igraph::E(g)$weight),
  n_components         = comp$no,
  giant_component_size = length(gc_nodes),
  avg_path_len_unw_giant = avg_path_unw,
  avg_path_len_w_giant   = avg_path_w,
  assortativity_degree = assort_deg,
  n_modules           = length(unique(igraph::V(g)$module)),
  modularity_louvain  = modularity_w,
  n_chembl_targets = sum(igraph::V(g)$is_chembl_target)
)


# ===== 8. Salidas =====
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# Métricas de nodo
metrics <- tibble(
  gene               = V(g)$name,
  degree             = V(g)$degree,
  strength           = V(g)$strength,
  betweenness        = V(g)$betweenness,
  betweenness_w      = V(g)$betweenness_w,
  closeness          = V(g)$closeness,
  closeness_w        = V(g)$closeness_w,
  closeness_harmonic = V(g)$closeness_harmonic,
  clustering_coef    = V(g)$clust_coef,
  module             = V(g)$module,
  is_chembl_target   = V(g)$is_chembl_target
) %>% arrange(desc(betweenness))

write.csv(metrics, metrics_path, row.names = FALSE)
write.csv(net_summary, summary_path, row.names = FALSE)

modules_df <- tibble(
  gene   = V(g)$name,
  module = as.integer(V(g)$module)
)
write.csv(modules_df, modules_path, row.names = FALSE)

saveRDS(g, graph_path)

cat("Métricas de red guardadas en:", metrics_path, "\n")
cat("Resumen de red guardado en:", summary_path, "\n")
cat("Módulos guardados en:", modules_path, "\n")
cat("Grafo guardado en:", graph_path, "\n")