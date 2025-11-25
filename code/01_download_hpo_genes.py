import sys
from pathlib import Path

import pandas as pd

# Fenotipo HPO de interés
HPO_ID = "HP:0030082"

# Fichero HPO: PHENOTYPE -> GENES (anotaciones propagadas, "inferred")
PHENO_URL = (
    "https://github.com/obophenotype/human-phenotype-ontology/"
    "releases/latest/download/phenotype_to_genes.txt"
)

print(f"Descargando anotaciones phenotype_to_genes desde:\n  {PHENO_URL}\n")

try:
    all_annots = pd.read_csv(
        PHENO_URL,
        sep="\t",
        comment="#",
        dtype=str
    )
except Exception as e:
    print("Error al descargar o leer phenotype_to_genes.txt")
    print(e)
    sys.exit(1)

print("Columnas detectadas en phenotype_to_genes.txt:")
print(list(all_annots.columns))

# Columna con el ID HPO
hpo_cols = [c for c in all_annots.columns if "hpo" in c.lower() and "id" in c.lower()]
if not hpo_cols:
    print("\nNo se ha encontrado una columna de ID HPO.")
    sys.exit(1)

hpo_col = hpo_cols[0]
print(f"\nUsando la columna '{hpo_col}' como identificador HPO.")

# Filtramos solo nuestro fenotipo HP:0030082
sub = all_annots[all_annots[hpo_col] == HPO_ID].copy()

if sub.empty:
    print(f"\nNo se han encontrado anotaciones para {HPO_ID} en phenotype_to_genes.txt.")
    sys.exit(1)

# Símbolo de gen
symbol_cols = [c for c in sub.columns if "symbol" in c.lower()]
if symbol_cols:
    gene_symbol_col = symbol_cols[0]
    sub["gene_symbol"] = sub[gene_symbol_col]
else:
    sub["gene_symbol"] = pd.NA
    print("\nAviso: no se ha encontrado columna de símbolo de gen, 'gene_symbol' se deja vacío.")

# ID Entrez (si existe alguna columna con 'entrez' o 'ncbi')
entrez_cols = [c for c in sub.columns if "entrez" in c.lower() or "ncbi" in c.lower()]
if entrez_cols:
    entrez_col = entrez_cols[0]
    sub["entrez_id"] = sub[entrez_col]
else:
    sub["entrez_id"] = pd.NA  # se podría mapear después en R con org.Hs.eg.db

# Construimos tabla de salida base
cols_out = ["gene_symbol", "entrez_id"]
if "hpo_id" in sub.columns:
    cols_out.append("hpo_id")
if "hpo_name" in sub.columns:
    cols_out.append("hpo_name")

df = (
    sub[cols_out]
    .drop_duplicates()
    .reset_index(drop=True)
)
df["source"] = "HPO phenotype_to_genes (inferred)"

print(f"\nSe han obtenido {len(df)} genes para {HPO_ID} a partir de phenotype_to_genes.txt.")


base_dir = Path(__file__).resolve().parents[1]
output_path = base_dir / "RESULTS" / "raw" / "hpo_genes.csv"
output_path.parent.mkdir(parents=True, exist_ok=True)

df.to_csv(output_path, index=False)

print(f"\nArchivo final guardado en: {output_path}")


