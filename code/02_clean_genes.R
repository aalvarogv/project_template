library(tidyverse)
library(org.Hs.eg.db)

raw <- read.csv("../../RESULTS/raw/hpo_genes.csv")

clean <- raw %>%
  distinct(gene_symbol, .keep_all = TRUE)

# añadir ENTREZ y ENSEMBL
ids <- mapIds(org.Hs.eg.db,
              keys = clean$gene_symbol,
              keytype = "SYMBOL",
              column = "ENTREZID",
              multiVals = "first")

clean$entrez <- ids

write.csv(clean, "../../RESULTS/processed/clean_genes.csv", row.names = FALSE)

print("Guardado en RESULTS/processed/clean_genes.csv")
