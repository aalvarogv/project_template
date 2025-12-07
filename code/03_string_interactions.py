#!/usr/bin/env python3
import sys, time, json
import pandas as pd
import requests
from tqdm import tqdm
from pathlib import Path
from urllib.parse import quote

CALLER_ID = "UMA_BioSis_AbnormalDrinking"
SPECIES   = 9606
REQUIRED_SCORE = 0.7   # umbral 0–1 para JSON
INDUCED_SUBGRAPH_ONLY = True  # True = solo pares dentro del set de genes HPO


# ===== 1. Rutas =====
base_dir = Path(__file__).resolve().parents[1]
genes_path = base_dir/"results"/"processed"/"clean_genes.csv"
edges_path = base_dir/"results"/"raw"/"string_edges.csv"
meta_path  = edges_path.with_suffix(".meta.json")

print(f"Leyendo genes desde: {genes_path}")
df = pd.read_csv(genes_path, dtype=str)
genes = (
    df["gene_symbol"]
    .dropna()
    .astype(str)
    .drop_duplicates()
    .tolist()
)
gene_set = set(genes)
print(f"[info] Total de genes a consultar en STRING: {len(genes)}")


# ===== 2. Sesión HTTP con reintentos =====
session = requests.Session()
adapter = requests.adapters.HTTPAdapter(max_retries=3)
session.mount("https://", adapter)
session.headers.update({"User-Agent": f"biosys-pipeline/1.0; {CALLER_ID}"})


# ===== 3. Llamada batch a STRING =====
# identifiers separados por %0d (carriage return URL-encoded)
identifiers = "%0d".join(quote(g) for g in genes)
url = (
    "https://string-db.org/api/json/network"
    f"?identifiers={identifiers}"
    f"&species={SPECIES}"
    f"&required_score={int(REQUIRED_SCORE*1000)}"   # STRING acepta 0-1000 en required_score
    f"&caller_identity={CALLER_ID}"
)

print("[info] Consultando STRING...")
try:
    r = session.get(url, timeout=60)
    r.raise_for_status()
    res = r.json()
except Exception as e:
    sys.exit(f"[error] STRING falló: {e}")


# ===== 4. Parsear y normalizar =====
results = []
for edge in res:
    gene1 = edge.get("preferredName_A")
    gene2 = edge.get("preferredName_B")
    score = edge.get("score", 0.0)
    if not gene1 or not gene2:
        continue
    # par no dirigido canónico
    n1, n2 = sorted((gene1, gene2))
    results.append({
        "gene1": n1,
        "gene2": n2,
        "score": float(score),
        "nscore": edge.get("neighborhood"),
        "fscore": edge.get("fusion"),
        "pscore": edge.get("phylogenetic"),
        "ascore": edge.get("coexpression"),
        "escore": edge.get("experimental"),
        "dscore": edge.get("database"),
        "tscore": edge.get("textmining"),
    })

edges_df= pd.DataFrame(results)

if edges_df.empty:
    print("[warn] No se han encontrado interacciones por encima del umbral.")
else:
    # nos quedamos con el máximo score por par
    edges_df = (
        edges_df.groupby(["gene1","gene2"], as_index=False)
                .agg({
                    "score":"max",
                    "nscore":"max","fscore":"max","pscore":"max","ascore":"max",
                    "escore":"max","dscore":"max","tscore":"max"
                })
                .sort_values("score", ascending=False)
    )
    edges_path.parent.mkdir(parents=True,exist_ok=True)
    edges_df.to_csv(edges_path,index=False)

print(f"[ok] Interacciones guardadas en: {edges_path}")
print(f"[ok] Total de aristas guardadas: {len(edges_df)}")



# ===== 5. ChEMBL: dianas farmacológicas =====
print("[info] Consultando ChEMBL...")
chembl_targets = []
for g in tqdm(genes, desc="ChEMBL"):
    url = f"https://www.ebi.ac.uk/chembl/api/data/target/search.json?q={g}"
    try:
        r = session.get(url, timeout=15)
        if r.status_code != 200:
            chembl_targets.append({"gene": g, "chembl_id": None, "pref_name": None})
            continue
        js = r.json()
        items = js.get("targets", [])
        # escoger el primero; si hay 'organism', preferir Homo sapiens
        best = None
        for it in items:
            if it.get("organism","").lower().startswith("homo sapiens"):
                best = it; break
        if best is None and items:
            best = items[0]
        chembl_targets.append({
            "gene": g,
            "chembl_id": best.get("target_chembl_id") if best else None,
            "pref_name": best.get("pref_name") if best else None
        })
        time.sleep(0.2)
    except Exception:
        chembl_targets.append({"gene": g, "chembl_id": None, "pref_name": None})


pd.DataFrame(chembl_targets).to_csv(
    base_dir / "results" / "raw" / "chembl_targets.csv",
    index=False
)

chembl_df = pd.DataFrame(chembl_targets)
chembl_path = base_dir / "results" / "raw" / "chembl_targets.csv"
chembl_df.to_csv(chembl_path, index=False)
print(f"[ok] Targets de ChEMBL guardados en: {chembl_path}")


# ===== 6. Metadatos =====
meta = {
    "string_api": "https://string-db.org/api/json/network",
    "species": SPECIES,
    "required_score": REQUIRED_SCORE,
    "induced_subgraph_only": INDUCED_SUBGRAPH_ONLY,
    "n_input_genes": len(genes),
    "n_edges": int(len(edges_df)),
    "chembl_api": "https://www.ebi.ac.uk/chembl/api/data/target/search.json",
    "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
}
meta_path.write_text(json.dumps(meta, indent=2))
print(f"[ok] Metadatos: {meta_path}")