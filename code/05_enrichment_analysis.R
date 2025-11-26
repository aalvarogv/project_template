library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

# Configuración de Rutas
base_dir <- dirname(getwd()) 
proc_dir <- file.path(base_dir, "RESULTS", "processed")
fig_dir  <- file.path(base_dir, "RESULTS", "figures")

# Crear carpeta de figuras si no existe
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Carga y Preparación de Datos
metrics_file <- file.path(proc_dir, "network_metrics.csv")

if(!file.exists(metrics_file)) stop("Error: No se encuentra network_metrics.csv")

metrics <- read.csv(metrics_file)

# Seleccionamos los 'Hubs' (top 20 por betweenness)
top_genes <- metrics %>%
  arrange(desc(betweenness)) %>%
  slice(1:20) %>%
  pull(gene)

cat("Analizando enriquecimiento para", length(top_genes), "genes principales...\n")

# Mapeo a ENTREZ
entrez <- mapIds(org.Hs.eg.db, keys = top_genes,
                 keytype = "SYMBOL",
                 column = "ENTREZID",
                 multiVals = "first")

# Análisis GO (Procesos Biológicos)
ego <- enrichGO(gene = entrez,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)

# Guardar tabla
write.csv(as.data.frame(ego), file.path(proc_dir, "enrichment_GO.csv"))

# Visualización
if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  p <- dotplot(ego, showCategory=15) + 
       ggtitle("Enriquecimiento Funcional (GO:BP)")
  
  # Guardamos el plot
  out_plot <- file.path(fig_dir, "enrichment_dotplot.png")
  ggsave(filename = out_plot, plot = p, width = 8, height = 6)
         
  cat("Plot guardado en:", out_plot, "\n")
} else {
  cat("No se encontraron términos significativos (p.adj < 0.05).\n")
}