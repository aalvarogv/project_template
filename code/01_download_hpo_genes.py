#!/usr/bin/env python3
import sys, json, time
from pathlib import Path
import argparse
import pandas as pd
import requests

HPO_ID = "HP:0030082" # Fenotipo HPO de interés

# Fichero HPO: PHENOTYPE -> GENES (anotaciones propagadas, "inferred")
PHENO_URL = (
    "https://github.com/obophenotype/human-phenotype-ontology/"
    "releases/latest/download/phenotype_to_genes.txt"
)


def download_table(url: str, max_retries: int = 3, sleep_s: int = 2) -> pd.DataFrame:
    # Descarga robusta con requests y luego lee con pandas desde el buffer
    last_err = None
    headers = {"User-Agent": "biosys-pipeline/1.0"}
    for attempt in range(1, max_retries + 1):
        try:
            r = requests.get(url, headers=headers, timeout=30)
            r.raise_for_status()
            # comprobación rápida de tamaño/contenido
            if len(r.content) < 1024:
                raise RuntimeError("Descarga demasiado pequeña; puede estar corrupta.")
            from io import StringIO
            buf = StringIO(r.text)
            df = pd.read_csv(buf, sep="\t", comment="#", dtype=str, engine="python")
            return df
        except Exception as e:
            last_err = e
            print(f"[warn] intento {attempt}/{max_retries} fallido: {e}")
            time.sleep(sleep_s)
    raise SystemExit(f"[error] no se pudo descargar/leér {url}: {last_err}")


def main():
    p = argparse.ArgumentParser(description="Descarga genes HPO para un término dado")
    p.add_argument("--hpo-id", default=HPO_ID, help="ID HPO (p.ej. HP:0030082)")
    p.add_argument("--url", default=PHENO_URL, help="URL TSV phenotype_to_genes.txt")
    p.add_argument("--out-csv", default=None, help="Ruta de salida CSV (por defecto results/raw/hpo_genes.csv)")
    p.add_argument("--metadata-json", default=None, help="Ruta para guardar metadatos JSON")
    args = p.parse_args()

    print(f"[info] Descargando anotaciones desde:\n  {args.url}\n")
    all_annots = download_table(args.url)

    # detectar columnas HPO
    hpo_cols = [c for c in all_annots.columns if "hpo" in c.lower() and "id" in c.lower()]
    if not hpo_cols:
        raise SystemExit("[error] no se encontró ninguna columna tipo 'HPO ID' en el TSV")
    hpo_col = hpo_cols[0]

    # filtrar por HPO
    sub = all_annots[all_annots[hpo_col] == args.hpo_id].copy()
    if sub.empty:
        raise SystemExit(f"[error] sin anotaciones para {args.hpo_id} en el TSV")

    # columnas símbolo y entrez
    symbol_cols = [c for c in sub.columns if "symbol" in c.lower()]
    entrez_cols = [c for c in sub.columns if ("entrez" in c.lower()) or ("ncbi" in c.lower())]

    sub["gene_symbol"] = sub[symbol_cols[0]] if symbol_cols else pd.NA
    sub["entrez_id"]   = sub[entrez_cols[0]] if entrez_cols else pd.NA

    cols_out = ["gene_symbol", "entrez_id"]
    # conserva nombre e id HPO si existen
    for k in ["hpo_id", "HPO-ID", "HPO ID"]:
        if k in sub.columns:
            sub.rename(columns={k: "hpo_id"}, inplace=True)
            break
    for k in ["hpo_name", "HPO-Name", "HPO Name"]:
        if k in sub.columns:
            sub.rename(columns={k: "hpo_name"}, inplace=True)
            break
    for k in ["hpo_id", "hpo_name"]:
        if k in sub.columns:
            cols_out.append(k)

    df = sub[cols_out].drop_duplicates().reset_index(drop=True)
    if df["gene_symbol"].isna().all():
        raise SystemExit("[error] no hay columna de símbolo de gen utilizable")

    df["source"] = "HPO phenotype_to_genes (inferred)"

    # ---- rutas de salida ----
    project_root = Path(__file__).resolve().parents[1]
    default_csv = project_root / "results" / "raw" / "hpo_genes.csv"
    out_csv = Path(args.out_csv) if args.out_csv else default_csv
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"[ok] {len(df)} genes para {args.hpo_id}. CSV: {out_csv}")

    # ---- metadatos ----
    md_path = Path(args.metadata_json) if args.metadata_json else (out_csv.with_suffix(".meta.json"))
    meta = {
        "hpo_id": args.hpo_id,
        "url": args.url,
        "n_genes": int(df["gene_symbol"].notna().sum()),
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "notes": "Anotaciones inferred desde phenotype_to_genes.txt (HPO release latest)."
    }
    md_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2))
    print(f"[ok] metadata JSON: {md_path}")


if __name__ == "__main__":
    main()