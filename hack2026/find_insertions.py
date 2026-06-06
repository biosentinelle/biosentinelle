import json
import pysam
import pandas as pd
from collections import defaultdict
from pathlib import Path

ALIGNMENT_DIR  = Path(__file__).parent / "alignment2"
SAM1           = ALIGNMENT_DIR / "SRR26617638.lite.1_1.sam"
SAM2           = ALIGNMENT_DIR / "SRR26617638.lite.1_2.sam"
CSV_PATH       = Path(__file__).parent / "ncbi_dataset_R64-1-1.csv"
OUTPUT_PATH    = Path(__file__).parent / "output" / "insertions.json"

MIN_JUNCTION_READS = 10  # discard (element, chromosome) pairs with fewer reads


def is_yeast(contig: str) -> bool:
    return contig.startswith("ref|NC_")


# --- Annotation (optional: flanking genes) ---
annotation = None
if CSV_PATH.exists():
    df = pd.read_csv(CSV_PATH, sep="\t")
    df["chrom_id"] = df["Accession"].str.replace(r"\.\d+$", "", regex=True)
    annotation = df

def flanking_genes(chrom: str, pos: int):
    if annotation is None:
        return None, None
    chrom_id = chrom.replace("ref|", "").replace("|", "")
    sub = annotation[annotation["chrom_id"] == chrom_id]
    up   = sub[sub["End"]   < pos].sort_values("End",   ascending=False).head(1)
    down = sub[sub["Begin"] > pos].sort_values("Begin", ascending=True ).head(1)
    left  = {"symbol": up.iloc[0]["Symbol"],   "locus": up.iloc[0]["Locus tag"],   "end":   int(up.iloc[0]["End"])}   if not up.empty   else None
    right = {"symbol": down.iloc[0]["Symbol"], "locus": down.iloc[0]["Locus tag"], "begin": int(down.iloc[0]["Begin"])} if not down.empty else None
    return left, right


# --- Pass 1: read names that mapped to any synthetic contig ---
print("Pass 1: collecting reads mapped to synthetic contigs...")
# read_name -> set of synthetic contig names
synthetic_reads: dict[str, set[str]] = defaultdict(set)

for path in (SAM1, SAM2):
    with pysam.AlignmentFile(str(path), "r") as sam:
        for read in sam:
            if not read.is_unmapped and not is_yeast(read.reference_name):
                synthetic_reads[read.query_name].add(read.reference_name)

print(f"  {len(synthetic_reads)} read names mapped to synthetic contigs")


# --- Pass 2: find their mates on yeast chromosomes ---
print("Pass 2: finding junction reads on yeast chromosomes...")
# junction_data[element][chromosome] = [0-based positions]
junction_data: dict[str, dict[str, list[int]]] = defaultdict(lambda: defaultdict(list))

for path in (SAM1, SAM2):
    with pysam.AlignmentFile(str(path), "r") as sam:
        for read in sam:
            if (read.query_name in synthetic_reads
                    and not read.is_unmapped
                    and is_yeast(read.reference_name)):
                for element in synthetic_reads[read.query_name]:
                    junction_data[element][read.reference_name].append(read.reference_start)


# --- Build results ---
print("Building results...")
results = []

for element in sorted(junction_data):
    for chrom in sorted(junction_data[element]):
        positions = sorted(junction_data[element][chrom])
        if len(positions) < MIN_JUNCTION_READS:
            continue

        # 1-based median as best insertion estimate
        insertion_pos = positions[len(positions) // 2] + 1
        left_gene, right_gene = flanking_genes(chrom, insertion_pos)

        results.append({
            "element":            element,
            "chromosome":         chrom,
            "insertion_position": insertion_pos,
            "junction_positions": [p + 1 for p in positions],  # 1-based
            "left_gene":          left_gene,
            "right_gene":         right_gene,
        })
        print(f"  {element:20s}  {chrom}  {len(positions):4d} reads  pos={insertion_pos}")

with open(OUTPUT_PATH, "w") as f:
    json.dump(results, f, indent=2)

print(f"\n{len(results)} insertion(s) written to {OUTPUT_PATH}")
