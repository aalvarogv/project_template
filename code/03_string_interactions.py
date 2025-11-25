import pandas as pd
import requests
from tqdm import tqdm
from pathlib import Path


base_dir = Path(__file__).resolve().parents[1]

genes_path = base_dir/"RESULTS"/"processed"/"clean_genes.csv"
edges_path = base_dir/"RESULTS"/"raw"/"string_edges.csv"

print(f"Leyendo genes desde: {genes_path}")
df = pd.read_csv(genes_path)

genes = (
    df["gene_symbol"]
    .dropna()
    .astype(str)
    .drop_duplicates()
    .tolist()
)

print(f"Total de genes a consultar en String: {len(genes)}")


CALLER_ID = "UMA_BioSis_AbnormalDrinking"
results = []

for g in tqdm(genes,desc="Consultando STRING"):
    url = (
        "https://string-db.org/api/json/network"
        f"?identifiers={g}"
        "&species=9606"
        "&caller_identity=" + CALLER_ID
    )
    
    try:
        r = requests.get(url, timeout=20)
        if r.status_code != 200:
            print(f"\n[AVISO] STRING devolvió status {r.status_code} para el gen {g}")
            continue

        res = r.json()
    except Exception as e:
        print(f"\n[ERROR] Problema al consultar STRING para {g}: {e}")
        continue

    for edge in res:
        score = edge.get("score",0)
        if score > 0.7:
            results.append({
                "gene1": edge.get("preferredName_A"),
                "gene2": edge.get("preferredName_B"),
                "score": score,
                "nscore":edge.get("neighborhood",None),
                "fscore": edge.get("fusion", None),
                "pscore": edge.get("phylogenetic", None),
                "ascore": edge.get("coexpression", None),
                "escore": edge.get("experimental", None),
                "dscore": edge.get("database", None),
                "tscore": edge.get("textmining", None),
                "combined_score": edge.get("combined_score", None),
            })

edges_df= pd.DataFrame(results)

if edges_df.empty:
    print("No se han encontrado interacciones con score > 0.7")
else:
    edges_df = edges_df.drop_duplicates(subset=["gene1","gene2","score"])
    edges_path.parent.mkdir(parents=True,exist_ok=True)
    edges_df.to_csv(edges_path,index=False)

print(f"Interacciones guardadas en: {edges_path}")
print(f"Total de aristas guardadas: {len(edges_df)}")


# ChEMBL: consultar si el gen es una diana farmacológica 
import time

chembl_targets = []

for g in genes:
    url = f"https://www.ebi.ac.uk/chembl/api/data/target/search.json?q={g}"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            js = r.json()
            if js["page_meta"]["total_count"] > 0:
                chembl_id = js["targets"][0]["target_chembl_id"]
                chembl_targets.append({"gene": g, "chembl_id": chembl_id})
            else:
                chembl_targets.append({"gene": g, "chembl_id": None})
    except:
        chembl_targets.append({"gene": g, "chembl_id": None})

    time.sleep(0.2)   # para no saturar ChEMBL


pd.DataFrame(chembl_targets).to_csv(
    base_dir / "RESULTS" / "raw" / "chembl_targets.csv",
    index=False
)

print("Targets de ChEMBL guardados en RESULTS/raw/chembl_targets.csv")

