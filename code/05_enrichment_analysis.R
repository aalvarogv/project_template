library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

metrics <- read.csv("../../RESULTS/processed/network_metrics.csv")

# ordenar por centralidad
top_genes <- metrics %>%
  arrange(desc(betweenness)) %>%
  slice(1:20) %>%
  pull(gene)

# convertir a ENTREZ
entrez <- mapIds(org.Hs.eg.db, keys = top_genes,
                 keytype = "SYMBOL",
                 column = "ENTREZID",
                 multiVals = "first")

ego <- enrichGO(gene = entrez,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)

write.csv(as.data.frame(ego), "../../RESULTS/processed/enrichment_GO.csv")
