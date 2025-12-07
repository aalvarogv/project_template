#!/usr/bin/env Rscript
library(igraph)
library(ggplot2)
library(tidyverse) 


# ===== 1. Rutas =====
base_dir <- dirname(getwd())
proc_dir <- file.path(base_dir, "results", "processed")
fig_dir  <- file.path(base_dir, "results", "figures")

# Comprobación rápida
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Carga de datos procesados
graph_path   <- file.path(proc_dir, "network_graph.rds")
metrics_path <- file.path(proc_dir, "network_metrics.csv")

if (!file.exists(graph_path)) stop("Falta el archivo .rds. Ejecuta primero el script 04.")
if (!file.exists(metrics_path)) stop("Falta network_metrics.csv. Ejecuta el script 04.")

g <- readRDS(graph_path)
metrics <- read.csv(metrics_path, stringsAsFactors = FALSE)

cat("Generando visualizaciones finales en:", fig_dir, "\n")


# ===== 2. Comprobaciones =====
need_cols <- c("gene","degree","betweenness","closeness")
if (!all(need_cols %in% names(metrics))) {
  stop("network_metrics.csv debe contener columnas: ", paste(need_cols, collapse=", "))
}


# ===== 3. Histograma de grado =====
p_hist <- ggplot(metrics, aes(x=degree)) +
  geom_histogram(binwidth = 1, fill="#4682B4", color="black", alpha=0.8) +
  theme_minimal() +
  labs(title = "Distribución de Grado (Conectividad)", 
       subtitle = paste("Nodos:", vcount(g), "| Aristas:", ecount(g)),
       x = "Grado (k)", 
       y = "Frecuencia")

ggsave(file.path(fig_dir, "degree_histogram.png"), p_hist, width = 6, height = 4)


# ===== 4. Red por módulos (Layout FR) =====
# Paleta discreta estable por módulo
set.seed(42)
modules <- if (!is.null(igraph::V(g)$module)) as.integer(igraph::V(g)$module) else NA_integer_
if (all(is.na(modules))) {
  node_colors <- rep("grey60", igraph::vcount(g))
  mod_map <- data.frame(module = integer(), color = character())
} else {
  u_mod <- sort(unique(modules))
  # paleta reproducible: usa hsv espaciado (evita dependencias extra)
  pal <- function(n) hsv(h = seq(0, 1 - 1/n, length.out = n), s = 0.55, v = 0.85)
  cols <- setNames(pal(length(u_mod)), u_mod)
  node_colors <- cols[as.character(modules)]
  mod_map <- data.frame(module = u_mod, color = cols, stringsAsFactors = FALSE)
  write.csv(mod_map, file.path(fig_dir, "module_colors.csv"), row.names = FALSE)
}

# layout y plot base igraph
lay <- igraph::layout_with_fr(g)
png(file.path(fig_dir, "network_modules_plot.png"), width = 1800, height = 1800, res = 220)
par(mar = c(1,1,2,1))
plot(
  g,
  layout       = lay,
  vertex.size  = 5,
  vertex.label = NA,
  vertex.color = node_colors,
  edge.color   = scales::alpha("grey20", 0.35),
  main         = "Red de Interacción Proteína–Proteína (Louvain)"
)
dev.off()

pdf(file.path(fig_dir, "network_modules_plot.pdf"), width = 10, height = 10)
par(mar = c(1,1,2,1))
plot(
  g,
  layout       = lay,
  vertex.size  = 5,
  vertex.label = NA,
  vertex.color = node_colors,
  edge.color   = scales::alpha("grey20", 0.35),
  main         = "Red de Interacción Proteína–Proteína (Louvain)"
)
dev.off()


# ===== 5. Heatmap de Hubs =====
top_hubs <- metrics %>%
  arrange(desc(betweenness)) %>%
  slice(1:10) %>%
  select(gene, degree, betweenness, closeness) %>%
  pivot_longer(cols = -gene, names_to = "metric", values_to = "value") %>%
  group_by(metric) %>%
  mutate(val_norm = (value - min(value)) / (max(value) - min(value))) # Normalización Min-Max

p_heat <- ggplot(top_hubs, aes(x = metric, y = reorder(gene, val_norm), fill = val_norm)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  theme_minimal(base_size = 11) +
  labs(title = "Top 10 Hubs: Perfil de Centralidad (normalizado por métrica)", 
       x = "", y = "", fill = "Score Norm.") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "hubs_heatmap.png"), p_heat, width = 5, height = 6)


# ===== 6. Dispersión grado vs betweenness =====
# Útil para visualizar hubs outliers
p_scatter <- ggplot(metrics, aes(x = degree, y = betweenness)) +
  geom_point(alpha = 0.6) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Hubs en la red: grado vs betweenness",
    x = "Grado",
    y = "Betweenness (normalizada)"
  )
ggsave(file.path(fig_dir, "degree_vs_betweenness.png"), p_scatter, width = 6, height = 4.5, dpi = 300)

cat("Proceso finalizado. Figuras generadas correctamente.\n")