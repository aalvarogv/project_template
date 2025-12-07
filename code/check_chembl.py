import pandas as pd
from pathlib import Path

base = Path(__file__).resolve().parents[1]

chembl = pd.read_csv(base / "results" / "raw" / "chembl_targets.csv")

print(chembl.head())
print("Genes con diana en ChEMBL:", chembl["chembl_id"].notna().sum())
