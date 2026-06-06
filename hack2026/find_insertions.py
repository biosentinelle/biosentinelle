import argparse
import json
import pysam
import pandas as pd
from collections import defaultdict
from pathlib import Path

CSV_PATH   = Path(__file__).parent / "ncbi_dataset_R64-1-1.csv"
OUTPUT_DIR = Path(__file__).parent / "output"

MIN_JUNCTION_READS = 10   # discard (element, chromosome) pairs below this
ZOOM_WINDOW        = 10_000  # bp — sliding window to find the enriched zone
MIN_GAP            = 500     # bp — minimum gap to be considered the cassette boundary


def is_yeast(contig: str) -> bool:
    return contig.startswith("ref|NC_")


# --- Annotation (optional) ---
annotation = None
if CSV_PATH.exists():
    df = pd.read_csv(CSV_PATH, sep="\t")
    df["chrom_id"] = df["Accession"].str.replace(r"\.\d+$", "", regex=True)
    annotation = df


def get_genes(chrom: str, left_bp: int, right_bp: int):
    """Return (candidate_deleted_gene, left_gene, right_gene) — all 1-based."""
    if annotation is None:
        return None, None, None

    chrom_id = chrom.replace("ref|", "").replace("|", "")
    sub = annotation[annotation["chrom_id"] == chrom_id]

    # Candidate: gene overlapping the breakpoint range; pick the one with most overlap
    overlapping = sub[(sub["Begin"] <= right_bp) & (sub["End"] >= left_bp)].copy()
    candidate = None
    if not overlapping.empty:
        overlapping["overlap"] = (
            overlapping["End"].clip(upper=right_bp)
            - overlapping["Begin"].clip(lower=left_bp)
        )
        row = overlapping.sort_values("overlap", ascending=False).iloc[0]
        candidate = {
            "symbol": row["Symbol"],
            "locus":  row["Locus tag"],
            "begin":  int(row["Begin"]),
            "end":    int(row["End"]),
        }

    # Left neighbor: nearest gene ending before the left breakpoint
    up = sub[sub["End"] < left_bp].sort_values("End", ascending=False).head(1)
    left_gene = (
        {"symbol": up.iloc[0]["Symbol"], "locus": up.iloc[0]["Locus tag"], "end": int(up.iloc[0]["End"])}
        if not up.empty else None
    )

    # Right neighbor: nearest gene starting after the right breakpoint
    down = sub[sub["Begin"] > right_bp].sort_values("Begin", ascending=True).head(1)
    right_gene = (
        {"symbol": down.iloc[0]["Symbol"], "locus": down.iloc[0]["Locus tag"], "begin": int(down.iloc[0]["Begin"])}
        if not down.empty else None
    )

    return candidate, left_gene, right_gene


def estimate_breakpoints(positions: list[int]) -> tuple[int, int, int]:
    """
    Find the enriched zone (10 kb sliding window) then the largest gap inside
    it. Returns (left_breakpoint, right_breakpoint, insertion_position) in
    1-based coordinates.
    Falls back to (min, max, median) of all positions if no gap is found.
    """
    # Step 1: enriched zone
    best_count, best_start, left = 0, positions[0], 0
    for right_idx, p in enumerate(positions):
        while positions[left] < p - ZOOM_WINDOW:
            left += 1
        if right_idx - left + 1 > best_count:
            best_count = right_idx - left + 1
            best_start = positions[left]

    zone = [p for p in positions if best_start <= p <= best_start + ZOOM_WINDOW]

    # Step 2: largest gap inside zone
    best_gap, left_bp, right_bp = 0, zone[0], zone[-1]
    for i in range(1, len(zone)):
        gap = zone[i] - zone[i - 1]
        if gap > best_gap and gap >= MIN_GAP:
            best_gap = gap
            left_bp, right_bp = zone[i - 1], zone[i]

    insertion_pos = (left_bp + right_bp) // 2 + 1
    return left_bp + 1, right_bp + 1, insertion_pos  # convert to 1-based


def main():
    parser = argparse.ArgumentParser(
        description="Detect engineered gene insertions from paired-end SAM files."
    )
    parser.add_argument("alignment_dir", help="Directory containing SAM files")
    parser.add_argument("--output", default=None, help="Output JSON path (default: output/insertions.json)")
    args = parser.parse_args()

    alignment_dir = Path(args.alignment_dir)
    sam_files = sorted(alignment_dir.glob("*.sam"))
    if not sam_files:
        print(f"No SAM files found in {alignment_dir}")
        return
    print(f"Found SAM files: {[f.name for f in sam_files]}")

    output_path = Path(args.output) if args.output else OUTPUT_DIR / "insertions.json"
    output_path.parent.mkdir(exist_ok=True)

    # --- Pass 1: read names mapped to any synthetic contig ---
    print("Pass 1: collecting reads mapped to synthetic contigs...")
    synthetic_reads: dict[str, set[str]] = defaultdict(set)
    for path in sam_files:
        with pysam.AlignmentFile(str(path), "r") as sam:
            for read in sam:
                if not read.is_unmapped and not is_yeast(read.reference_name):
                    synthetic_reads[read.query_name].add(read.reference_name)
    print(f"  {len(synthetic_reads)} read names mapped to synthetic contigs")

    # --- Pass 2: junction reads on yeast chromosomes ---
    print("Pass 2: finding junction reads on yeast chromosomes...")
    junction_data: dict[str, dict[str, list[int]]] = defaultdict(lambda: defaultdict(list))
    for path in sam_files:
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

            left_bp, right_bp, insertion_pos = estimate_breakpoints(positions)
            candidate, left_gene, right_gene = get_genes(chrom, left_bp, right_bp)

            results.append({
                "element":                element,
                "chromosome":             chrom,
                "left_breakpoint":        left_bp,
                "right_breakpoint":       right_bp,
                "insertion_position":     insertion_pos,
                "candidate_deleted_gene": candidate,
                "left_gene":              left_gene,
                "right_gene":             right_gene,
                "junction_positions":     [p + 1 for p in positions],
            })
            print(f"  {element:20s}  {chrom}  {len(positions):4d} reads  [{left_bp} – {right_bp}]")

    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n{len(results)} insertion(s) written to {output_path}")


if __name__ == "__main__":
    main()
