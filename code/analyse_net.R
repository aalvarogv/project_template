#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(igraph)
})

option_list <- list(
  make_option(c("-i","--input_file"), type="character", help="Fichero de aristas (STRING filtrado)"),
  make_option(c("-g","--graph_file"), type="character", help="Salida PDF de la red"),
  make_option(c("-d","--dict_file"),  type="character", help="Diccionario protein_id <-> HGNC_symbol"),
  make_option(c("-G","--gene_subset"),type="character", help="Lista de genes (símbolos HGNC) del fenotipo HPO")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(file.exists(opt$input_file), file.exists(opt$dict_file), file.exists(opt$gene_subset))

# --- Carga de datos ---
net_table <- read.table(opt$input_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(net_table) <- c("protein1","protein2","combined_score")

gene_list <- read.table(opt$gene_subset, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# Si el archivo tiene una sola columna sin cabecera, fuerza nombre:
if (ncol(gene_list) == 1 && is.null(colnames(gene_list))) {
  colnames(gene_list) <- "name"
}

dict_gene <- read.table(opt$dict_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(dict_gene) <- c("protein_id","hgnc_symbol")

# --- Mapea símbolos HPO -> protein_id STRING ---
mv <- match(gene_list$name, dict_gene$hgnc_symbol)
prot_ids <- dict_gene$protein_id[mv]
prot_ids <- unique(na.omit(prot_ids))
if (length(prot_ids) < 2) {
  stop("Muy pocos genes mapeados a STRING. Revisa el diccionario o el HPO_ID.")
}

# --- Construye grafo y reduce a subgrafo inducido por tus genes ---
g_all <- graph_from_data_frame(d = net_table[, c("protein1","protein2")], directed = FALSE)
# etiqueta los vértices con protein_id
V(g_all)$name <- V(g_all)$name

nodes_in_graph <- intersect(V(g_all)$name, prot_ids)
g <- induced_subgraph(g_all, vids = nodes_in_graph)

# Simplifica (combina paralelas, elimina loops)
g <- igraph::simplify(
  g, remove.multiple = TRUE, remove.loops = TRUE,
  edge.attr.comb = list(weight = "max", .default = "first")
)

cat("Grafo inducido (HPO):\n")
cat("  Nodos:", vcount(g), "\n")
cat("  Aristas:", ecount(g), "\n")

# Métrica global
trans_glob <- transitivity(g, type = "global", isolates = "zero")
cat("  Transitivity (global):", sprintf("%.4f", trans_glob), "\n")

# --- Dibuja a PDF ---
pdf(file = opt$graph_file, width = 10, height = 10)
set.seed(123)
lay <- layout_with_fr(g)
plot(
  g,
  layout       = lay,
  vertex.color = "orange",
  vertex.size  = 5,
  vertex.label = NA,
  edge.width   = 0.5,
  main         = "STRING subgraph (HPO gene set) — Force-directed"
)
dev.off()

cat("[ok] PDF guardado en:", opt$graph_file, "\n")
