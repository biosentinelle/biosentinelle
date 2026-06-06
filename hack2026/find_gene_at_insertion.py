import json
import pandas as pd
from pathlib import Path

JSON_PATH = Path(__file__).parent / "insertion_site.json"
CSV_PATH  = Path(__file__).parent / "ncbi_dataset_R64-1-1.csv"

# --- Load insertion site ---
if JSON_PATH.exists():
    with open(JSON_PATH) as f:
        site = json.load(f)
    INSERTION_CHROM = site["chromosome"]
    INSERTION_POS   = site["estimated_position"]
    print(f"Loaded from {JSON_PATH.name}: {INSERTION_CHROM}  pos={INSERTION_POS}")
else:
    INSERTION_CHROM = "ref|NC_001136|"
    INSERTION_POS   = 1342334
    print("No insertion_site.json found — using hardcoded fallback values.")

# --- Load annotation ---
df = pd.read_csv(CSV_PATH, sep="\t")

# The CSV accession looks like "NC_001136.9"; the SAM uses "ref|NC_001136|".
chrom_id = INSERTION_CHROM.replace("ref|", "").replace("|", "")
df["chrom_id"] = df["Accession"].str.replace(r"\.\d+$", "", regex=True)

chrom_df = df[df["chrom_id"] == chrom_id].copy()
print(f"Genes on {chrom_id}: {len(chrom_df)}")

# --- Look for overlap ---
overlapping = chrom_df[(chrom_df["Begin"] <= INSERTION_POS) & (chrom_df["End"] >= INSERTION_POS)]

if not overlapping.empty:
    print(f"\nInsertion at {INSERTION_POS} falls INSIDE {len(overlapping)} gene(s):")
    print(overlapping[["Symbol", "Name", "Begin", "End", "Orientation", "Locus tag", "Gene Type"]].to_string(index=False))
else:
    print(f"\nNo gene directly overlaps position {INSERTION_POS}. Finding nearest neighbors...")

    upstream   = chrom_df[chrom_df["End"]   < INSERTION_POS].sort_values("End",   ascending=False).head(1)
    downstream = chrom_df[chrom_df["Begin"] > INSERTION_POS].sort_values("Begin", ascending=True).head(1)

    if not upstream.empty:
        dist = INSERTION_POS - upstream.iloc[0]["End"]
        print(f"\nNearest gene UPSTREAM (ends {dist} bp before insertion):")
        print(upstream[["Symbol", "Name", "Begin", "End", "Orientation", "Locus tag"]].to_string(index=False))

    if not downstream.empty:
        dist = downstream.iloc[0]["Begin"] - INSERTION_POS
        print(f"\nNearest gene DOWNSTREAM (starts {dist} bp after insertion):")
        print(downstream[["Symbol", "Name", "Begin", "End", "Orientation", "Locus tag"]].to_string(index=False))
