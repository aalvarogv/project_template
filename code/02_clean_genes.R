#!/usr/bin/env Rscript
library(tidyverse)
library(org.Hs.eg.db)

# ===== 1. Rutas =====
base_dir <- dirname(getwd())
raw_path <- file.path(base_dir, "RESULTS", "raw", "hpo_genes.csv")
out_dir  <- file.path(base_dir, "RESULTS", "processed")
out_path <- file.path(out_dir, "clean_genes.csv")

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
# mapIds devuelve NA si no encuentra el gen
clean$entrez_id <- mapIds(
  org.Hs.eg.db,
  keys      = clean$gene_symbol,
  column    = "ENTREZID",
  keytype   = "SYMBOL",
  multiVals = "first"
)

clean$ensembl_id <- mapIds(
  org.Hs.eg.db,
  keys      = clean$gene_symbol,
  column    = "ENSEMBL",
  keytype   = "SYMBOL",
  multiVals = "first"
)

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

cat("Guardado en:", out_path, "\n")
cat("Genes totales:", nrow(clean), "\n")
cat("Genes con ENTREZ:", sum(!is.na(clean$entrez_id)), "\n")

