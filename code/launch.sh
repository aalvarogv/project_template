#!/usr/bin/env bash
set -euo pipefail

# Variables
export R_LIBS_USER=./R_LIBS
BASE_DIR="$(cd "$(dirname "$0")" && cd .. && pwd)"
DATA_DIR="${BASE_DIR}/results/data"
RES_DIR="${BASE_DIR}/results"
mkdir -p "${DATA_DIR}" "${RES_DIR}"

# === 1) Descargas STRING v12.0 ===
echo "[info] Descargando STRING (edges, aliases) ..."
wget -q "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz" \
     -O "${DATA_DIR}/string.links.txt.gz"
gunzip -f "${DATA_DIR}/string.links.txt.gz"

# Filtra por combined_score > 900 (alta confianza)
awk -F'\t' 'NR==1 {next} $3>900 {print $1 "\t" $2 "\t" $3}' \
  "${DATA_DIR}/string.links.txt" > "${DATA_DIR}/string_network.txt"

wget -q "https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz" \
     -O "${DATA_DIR}/aliases.txt.gz"
gunzip -f "${DATA_DIR}/aliases.txt.gz"

# Diccionario: protein_id <tab> HGNC_symbol (fuente = Ensembl_HGNC_symbol)
grep -P '\tEnsembl_HGNC_symbol$' "${DATA_DIR}/aliases.txt" | cut -f1,2 > "${DATA_DIR}/string_dict.txt"

# === 2) Genes HPO del fenotipo ===
# Usa tu fenotipo (por defecto HP:0030082 – Abnormal drinking behaviour)
HPO_ID="${HPO_ID:-HP:0030082}"
echo "[info] Descargando genes HPO para ${HPO_ID} ..."
wget -q "https://ontology.jax.org/api/network/annotation/${HPO_ID}/download/gene" \
     -O "${DATA_DIR}/hpo_query.txt"

# Mantén solo la segunda columna (gene symbol), y limpia espacios
cut -f2 "${DATA_DIR}/hpo_query.txt" | tr -d ' ' > "${DATA_DIR}/hpo_gene.txt"

# === 3) Ejecuta el análisis y exporta PDF ===
echo "[info] Ejecutando analyse_net.R ..."
Rscript "$(dirname "$0")/analyse_net.R" \
  -i "${DATA_DIR}/string_network.txt" \
  -g "${RES_DIR}/net.pdf" \
  -G "${DATA_DIR}/hpo_gene.txt" \
  -d "${DATA_DIR}/string_dict.txt"

echo "[ok] Hecho. Revisa:"
echo " - Red (PDF): ${RES_DIR}/net.pdf"
echo " - Log en consola (nº nodos/aristas, transitivity)"
