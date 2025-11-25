import pandas as pd
import requests
from tqdm import tqdm

df = pd.read_csv("../../RESULTS/processed/clean_genes.csv")
genes = df["gene_symbol"].tolist()

results = []

for g in tqdm(genes):
    url = f"https://string-db.org/api/json/network?identifiers={g}&species=9606"
    res = requests.get(url).json()

    for edge in res:
        if edge["score"] > 0.7:
            results.append({
                "gene1": edge["preferredName_A"],
                "gene2": edge["preferredName_B"],
                "score": edge["score"]
            })

pd.DataFrame(results).to_csv("../../RESULTS/raw/string_edges.csv", index=False)
print("Interacciones guardadas en RESULTS/raw/string_edges.csv")
