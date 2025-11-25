import requests
import pandas as pd

HPO_ID = "HP:0030082"
URL = f"https://hpo.jax.org/api/hpo/term/{HPO_ID}"

print(f"Descargando información para {HPO_ID}...")

data = requests.get(URL).json()

genes = []

for g in data["genes"]:
    genes.append({
        "gene_symbol": g["geneSymbol"],
        "entrez_id": g.get("entrezGeneId", None),
        "source": g.get("source", None)
    })

df = pd.DataFrame(genes)
df.to_csv("../../RESULTS/raw/hpo_genes.csv", index=False)

print("Guardado en RESULTS/raw/hpo_genes.csv")
