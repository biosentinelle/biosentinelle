import json
import pysam
from collections import Counter
from pathlib import Path

ALIGNMENT_DIR = Path(__file__).parent / "alignment"
SAM1 = ALIGNMENT_DIR / "SRR26617638_1.sam"
SAM2 = ALIGNMENT_DIR / "SRR26617638_2.sam"

# --- Step 1: collect read names that mapped to KanMX in each file ---
print("Step 1: scanning for KanMX-mapped reads...")

kan_names_1 = set()
with pysam.AlignmentFile(str(SAM1), "r") as sam:
    for read in sam:
        if not read.is_unmapped and read.reference_name == "KanMX":
            kan_names_1.add(read.query_name)

kan_names_2 = set()
with pysam.AlignmentFile(str(SAM2), "r") as sam:
    for read in sam:
        if not read.is_unmapped and read.reference_name == "KanMX":
            kan_names_2.add(read.query_name)

print(f"  KanMX reads in _1.sam: {len(kan_names_1)}")
print(f"  KanMX reads in _2.sam: {len(kan_names_2)}")

# --- Step 2: find where the MATES of those reads landed on yeast chromosomes ---
# The two SAM files were aligned independently, so mate info is absent.
# We cross-reference by read name: a read whose pair mapped to KanMX but
# which itself maps to a yeast chromosome is a "junction read" — it sits
# right at the boundary of the insertion site.
print("\nStep 2: finding junction reads on yeast chromosomes...")

junction_positions = []  # (chromosome, position)

# mates of _1.sam KanMX reads are in _2.sam
with pysam.AlignmentFile(str(SAM2), "r") as sam:
    for read in sam:
        if (read.query_name in kan_names_1
                and not read.is_unmapped
                and read.reference_name != "KanMX"):
            junction_positions.append((read.reference_name, read.reference_start))

# mates of _2.sam KanMX reads are in _1.sam
with pysam.AlignmentFile(str(SAM1), "r") as sam:
    for read in sam:
        if (read.query_name in kan_names_2
                and not read.is_unmapped
                and read.reference_name != "KanMX"):
            junction_positions.append((read.reference_name, read.reference_start))

print(f"  Junction reads found: {len(junction_positions)}")

# --- Step 3: find where the positions cluster ---
print("\nStep 3: clustering positions to locate insertion site...")

chrom_counts = Counter(chrom for chrom, _ in junction_positions)
print("\nTop chromosomes near KanMX insertion:")
for chrom, count in chrom_counts.most_common(5):
    print(f"  {chrom}: {count} reads")

top_chrom = chrom_counts.most_common(1)[0][0]
all_positions = sorted(pos for chrom, pos in junction_positions if chrom == top_chrom)

# Bin positions into small windows and find the densest one.
WINDOW = 200
bin_counts = Counter(pos // WINDOW for pos in all_positions)
best_bin = bin_counts.most_common(1)[0][0]
best_bin_count = bin_counts.most_common(1)[0][1]

cluster_positions = [p for p in all_positions if p // WINDOW == best_bin]

# pysam uses 0-based positions; add 1 to convert to standard 1-based genomic coordinates
estimated_pos = cluster_positions[len(cluster_positions) // 2] + 1

print(f"\n=== Result ===")
print(f"Chromosome : {top_chrom}")
print(f"Densest {WINDOW} bp window: {best_bin * WINDOW + 1} – {(best_bin + 1) * WINDOW} bp  ({best_bin_count} reads)")
print(f"Cluster range (1-based): {cluster_positions[0] + 1} – {cluster_positions[-1] + 1}")
print(f"Estimated insertion position (1-based): {estimated_pos}")

print(f"\nTop 10 densest {WINDOW} bp windows on {top_chrom}:")
print(f"  {'Window start':>12}  {'Window end':>10}  {'Reads':>6}")
for b, cnt in bin_counts.most_common(10):
    print(f"  {b * WINDOW + 1:>12}  {(b + 1) * WINDOW:>10}  {cnt:>6}")

# --- Save result for use by find_gene_at_insertion.py ---
out_path = Path(__file__).parent / "output" / "insertion_site.json"
with open(out_path, "w") as f:
    json.dump({"chromosome": top_chrom, "estimated_position": estimated_pos}, f, indent=2)
print(f"\nResult saved to {out_path}")
