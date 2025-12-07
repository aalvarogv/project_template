#!/usr/bin/env Rscript
library(tidyverse)
library(gprofiler2)
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

# ===== 1. Rutas =====
base_dir <- dirname(getwd()) 
proc_dir <- file.path(base_dir, "results", "processed")
fig_dir  <- file.path(base_dir, "results", "figures")

# Crear carpeta de figuras si no existe
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Carga y Preparación de Datos
metrics_file <- file.path(proc_dir, "network_metrics.csv")

if(!file.exists(metrics_file)) stop("Error: No se encuentra ", metrics_file)

metrics <- read.csv(metrics_file, stringsAsFactors = FALSE)


# ===== 2. Conjunto de prueba: top hubs =====
top_genes <- metrics %>%
  arrange(desc(betweenness)) %>%
  slice(1:20) %>%
  pull(gene) %>%
  unique()

bg_genes <- unique(metrics$gene)

if (length(top_genes) < 3) stop("Muy pocos genes en el top para hacer enriquecimiento.")
if (length(bg_genes)  < 10) stop("Muy pocos genes de fondo; revisa metrics_file.")

cat("[info] Hubs (n) :", length(top_genes), "\n")
cat("[info] Fondo (n):", length(bg_genes),  "\n")

# ===== 3. Enriquecimiento con g:Profiler (GO:BP) =====
# - organism: hsapiens
# - correction_method: "fdr" (Benjamini–Hochberg)
# - sources: "GO:BP"
# - custom_bg: universo = genes de toda tu red
gp <- gost(
  query = list(hubs = top_genes),
  organism = "hsapiens",
  significant = FALSE,
  correction_method = "fdr",
  sources = c("GO:BP"),
  custom_bg = bg_genes,
  domain_scope = "custom"
)

# Guardar tabla completa y una versión ordenada por p_adj
out_full <- file.path(proc_dir, "enrichment_GO_full.csv")
out_top  <- file.path(proc_dir, "enrichment_GO.csv")

if (!is.null(gp) && !is.null(gp$result) && nrow(gp$result) > 0) {
  res <- gp$result %>%
    # p_value ya viene ajustado (FDR) con correction_method = "fdr"
    dplyr::mutate(p_adj = p_value, .before = term_name) %>%
    # aplastar columnas de lista a texto para escribir CSV sin errores
    dplyr::mutate(dplyr::across(where(is.list), ~sapply(., function(x) paste(x, collapse = ",")))) %>%
    dplyr::arrange(p_adj)
  
  # filtrado razonable dado el fondo pequeño: p_adj <= 0.1 (ajústalo si quieres 0.05)
  res_sig <- res %>% dplyr::filter(p_adj <= 0.1)

  utils::write.csv(res,     out_full, row.names = FALSE)
  utils::write.csv(res_sig, out_top,  row.names = FALSE)

  if (nrow(res_sig) > 0) {
    cat("[ok] Tablas guardadas en:\n  -", out_full, "\n  -", out_top, " (p_adj <= 0.1)\n")
  } else {
    cat("[warn] Sin términos con p_adj <= 0.1; revisa", out_full, "para ver el ranking completo.\n")
  }

  # Dotplot
  if (nrow(res) > 0) {
    topN <- res %>% dplyr::slice_min(order_by = p_adj, n = 15) %>%
      dplyr::mutate(
        term_name     = factor(term_name, levels = rev(unique(term_name))),
        minus_log10_p = -log10(p_adj)
      )

    p <- ggplot2::ggplot(topN, ggplot2::aes(x = minus_log10_p, y = term_name, size = intersection_size)) +
      ggplot2::geom_point() +
      ggplot2::labs(title = "Enriquecimiento GO:BP (g:Profiler — fondo personalizado)",
                    x = expression(-log[10](FDR)),
                    y = NULL, size = "Genes\nintersección") +
      ggplot2::theme_minimal(base_size = 11)

    if (has_ggrepel) {
      p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = term_id),
                                        size = 3, min.segment.length = 0.1, max.overlaps = Inf)
    }

    out_plot <- file.path(fig_dir, "enrichment_dotplot.png")
    ggplot2::ggsave(out_plot, p, width = 8, height = 6)
    cat("[ok] Figura de dotplot guardada en:", out_plot, "\n")
  }
} else {
  cat("[warn] g:Profiler no devolvió filas (con fondo personalizado y FDR). Exporta la lista de hubs y considera relajar umbral.\n")
}


# ===== 4. Enriquecimiento por módulos Louvain =====         
if ("module" %in% names(metrics)) {
  mod_dir <- file.path(proc_dir, "module_enrichment")
  if (!dir.exists(mod_dir)) dir.create(mod_dir, recursive = TRUE, showWarnings = FALSE)

  for (m in sort(unique(metrics$module))) {
    genes_m <- metrics %>% filter(module == m) %>% pull(gene) %>% unique()
    if (length(genes_m) < 5) next

    gp_m <- gost(
      query            = list(set = genes_m),
      organism         = "hsapiens",
      significant      = FALSE,
      correction_method= "fdr",
      sources          = c("GO:BP"),
      custom_bg        = bg_genes,
      domain_scope      = "custom"
    )

    if (!is.null(gp_m) && !is.null(gp_m$result) && nrow(gp_m$result) > 0) {
      res_m <- gp_m$result %>%
        dplyr::mutate(p_adj = p_value) %>%
        dplyr::mutate(dplyr::across(where(is.list), ~sapply(., function(x) paste(x, collapse = ",")))) %>%
        dplyr::arrange(p_adj)
      
      out_mod_full <- file.path(mod_dir, paste0("module_", m, "_GO_full.csv"))
      out_mod_sig  <- file.path(mod_dir, paste0("module_", m, "_GO.csv"))
      
      utils::write.csv(res_m, out_mod_full, row.names = FALSE)
      utils::write.csv(res_m %>% dplyr::filter(p_adj <= 0.1), out_mod_sig, row.names = FALSE)
      
      cat("[ok] Módulo", m, "→", out_mod_sig, "(p_adj <= 0.1); full en", out_mod_full, "\n")
    } else {
      cat("[info] Módulo", m, ": sin términos GO:BP (FDR) con fondo pequeño.\n")
    }
  }
}

# ===== 5. Metadata del análisis =====
meta <- list(
  date_time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  organism = "Homo sapiens",
  db = "GO:BP",
  method = "FDR (Benjamini–Hochberg)",
  query_size = length(top_genes),
  background_n = length(bg_genes),
  tool = paste0("gprofiler2_", as.character(utils::packageVersion("gprofiler2")))

)
jsonlite::write_json(meta, file.path(proc_dir, "enrichment_meta.json"), pretty = TRUE, auto_unbox = TRUE)

cat("Enriquecimiento finalizado. Tablas en results/processed y figuras en results/figures.\n")