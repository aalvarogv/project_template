library(igraph)
library(ggplot2)
library(tidyverse) 

# Configuración de rutas
base_dir <- dirname(getwd())
proc_dir <- file.path(base_dir, "RESULTS", "processed")
fig_dir  <- file.path(base_dir, "RESULTS", "figures")

# Comprobación rápida
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Carga de datos procesados
graph_path   <- file.path(proc_dir, "network_graph.rds")
metrics_path <- file.path(proc_dir, "network_metrics.csv")

if (!file.exists(graph_path)) stop("Falta el archivo .rds. Ejecuta primero el script 04.")

g <- readRDS(graph_path)
metrics <- read.csv(metrics_path)

cat("Generando visualizaciones finales en:", fig_dir, "\n")

# FIGURA 1: Histograma de Conectividad
p_hist <- ggplot(metrics, aes(x=degree)) +
  geom_histogram(binwidth = 1, fill="#4682B4", color="black", alpha=0.8) +
  theme_minimal() +
  labs(title = "Distribución de Grado (Conectividad)", 
       subtitle = paste("Nodos:", vcount(g), "| Aristas:", ecount(g)),
       x = "Grado (k)", 
       y = "Frecuencia")

ggsave(file.path(fig_dir, "degree_histogram.png"), p_hist, width = 6, height = 4)

# FIGURA 2: Red Topológica (Módulos)
node_colors <- if(!is.null(V(g)$module)) V(g)$module else "grey50"

pdf(file.path(fig_dir, "network_modules_plot.pdf"), width=10, height=10)
set.seed(42)
plot(g, 
     vertex.size = 5, 
     vertex.label = NA, 
     vertex.color = node_colors, 
     edge.color = alpha("grey", 0.4),
     layout = layout_with_fr(g),
     main = "Red de Interacción Proteína-Proteína (Clusterizada)")
dev.off()


# FIGURA 3: Heatmap de Hubs
top_hubs <- metrics %>%
  arrange(desc(betweenness)) %>%
  head(10) %>%
  select(gene, degree, betweenness, closeness) %>%
  pivot_longer(cols = -gene, names_to = "metric", values_to = "value") %>%
  group_by(metric) %>%
  mutate(val_norm = (value - min(value)) / (max(value) - min(value))) # Normalización Min-Max

p_heat <- ggplot(top_hubs, aes(x = metric, y = reorder(gene, val_norm), fill = val_norm)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  theme_minimal() +
  labs(title = "Top 10 Hubs: Perfil de Centralidad", 
       x = "", y = "", fill = "Score Norm.") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "hubs_heatmap.png"), p_heat, width = 5, height = 6)

cat("Proceso finalizado. Figuras generadas correctamente.\n")