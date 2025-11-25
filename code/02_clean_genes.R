library(tidyverse)

base_dir <- dirname(getwd())

raw_path <- file.path(base_dir, "RESULTS", "raw", "hpo_genes.csv")
out_dir  <- file.path(base_dir, "RESULTS", "processed")
out_path <- file.path(out_dir, "clean_genes.csv")

raw <- read.csv(raw_path, stringsAsFactors = FALSE)

clean <- raw %>%
  distinct(gene_symbol, .keep_all = TRUE)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

write.csv(clean, out_path, row.names = FALSE)

cat("Guardado en RESULTS/processed/clean_genes.csv\n")
cat("Genes totales:", nrow(clean), "\n")
