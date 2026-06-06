import pandas as pd
from pathlib import Path

# --- Insertion site from find_kan_insertion.py output ---
# Paste the chromosome and estimated insertion position here
INSERTION_CHROM = "ref|NC_001136|"   # top chromosome from the SAM analysis
INSERTION_POS   = 1344617            # estimated insertion position (1-based, median of densest bin)

CSV_PATH = Path(__file__).parent / "ncbi_dataset_R64-1-1.csv"

# --- Load annotation ---
df = pd.read_csv(CSV_PATH, sep="\t")

# The CSV accession looks like "NC_001136.9"; the SAM uses "ref|NC_001136|".
# Strip the version suffix and "ref|...|" wrapper so both sides match.
chrom_id = INSERTION_CHROM.replace("ref|", "").replace("|", "")  # → NC_001136
df["chrom_id"] = df["Accession"].str.replace(r"\.\d+$", "", regex=True)  # NC_001136.9 → NC_001136

# --- Filter to the right chromosome ---
chrom_df = df[df["chrom_id"] == chrom_id].copy()
print(f"Genes on {chrom_id}: {len(chrom_df)}")

# --- Look for overlap: insertion point falls inside a gene's Begin-End range ---
overlapping = chrom_df[(chrom_df["Begin"] <= INSERTION_POS) & (chrom_df["End"] >= INSERTION_POS)]

if not overlapping.empty:
    print(f"\nInsertion at {INSERTION_POS} falls INSIDE {len(overlapping)} gene(s):")
    print(overlapping[["Symbol", "Name", "Begin", "End", "Orientation", "Locus tag", "Gene Type"]].to_string(index=False))
else:
    # --- Fallback: find the nearest gene upstream and downstream ---
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
