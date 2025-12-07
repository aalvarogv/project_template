#!/usr/bin/env Rscript
library(tidyverse)
library(org.Hs.eg.db)

# ===== 1. Rutas =====
base_dir <- dirname(getwd())
raw_path <- file.path(base_dir, "results", "raw", "hpo_genes.csv")
out_dir  <- file.path(base_dir, "results", "processed")
out_path <- file.path(out_dir, "clean_genes.csv")
out_map  <- file.path(out_dir,  "id_mapping_long.csv")
meta_js  <- file.path(out_dir,  "clean_genes.meta.json")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# ===== 2. Cargar datos =====
raw <- read.csv(raw_path, stringsAsFactors = FALSE)

# Quitar duplicados por símbolo
clean <- raw %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  distinct(gene_symbol, .keep_all = TRUE)


# ===== 3. Anotar IDs usando org.Hs.eg.db =====
gene_vec <- unique(clean$gene_symbol)

# mapIds devuelve NA si no encuentra el gen
clean$entrez_id <- mapIds(
  org.Hs.eg.db,
  keys      = clean$gene_symbol,
  column    = "ENTREZID",
  keytype   = "SYMBOL",
  multiVals = "first"
)[clean$gene_symbol]

clean$ensembl_id <- mapIds(
  org.Hs.eg.db,
  keys      = clean$gene_symbol,
  column    = "ENSEMBL",
  keytype   = "SYMBOL",
  multiVals = "first"
)[clean$gene_symbol]

# ===== 4. Registrar genes no encontrados =====
not_found <- clean %>%
  filter(is.na(entrez_id) & is.na(ensembl_id)) %>%
  pull(gene_symbol)

if (length(not_found) > 0) {
  writeLines(not_found, file.path(out_dir, "genes_no_reconocidos.txt"))
  cat("Genes sin anotación:", length(not_found), "\n")
}

# ===== 5. Guardar tabla limpia =====
write.csv(clean, out_path, row.names = FALSE)

mapping_long <- bind_rows(
  tibble(gene_symbol = clean$gene_symbol, id_type = "ENTREZID", id_value = clean$entrez_id),
  tibble(gene_symbol = clean$gene_symbol, id_type = "ENSEMBL",  id_value = clean$ensembl_id)
)
write.csv(mapping_long, out_map, row.names = FALSE)

meta <- list(
  n_input = nrow(raw),
  n_unique_symbols = nrow(clean),
  n_with_entrez = sum(!is.na(clean$entrez_id)),
  n_with_ensembl = sum(!is.na(clean$ensembl_id)),
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
write_json(meta, meta_js, pretty = TRUE, auto_unbox = TRUE)


message("[ok] Guardado: ", out_csv)
message("     Mapeos: ", out_map)
message("     Meta:   ", meta_js)
