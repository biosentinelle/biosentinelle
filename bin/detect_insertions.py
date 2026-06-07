#!/usr/bin/env python3

import argparse
import csv
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

MIN_JUNCTION_READS = 10
ZOOM_WINDOW = 10_000
MIN_GAP = 500


def run(cmd):
    return subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        errors="replace",
    )


def read_manifest(path):
    alignments = []
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 2:
                raise ValueError(f"Invalid manifest line in {path}: {row}")
            alignments.append((row[0], row[1]))
    return alignments


def sam_records(path):
    path = str(path)
    if path.endswith(".bam"):
        proc = run(["samtools", "view", "-h", path])
        stream = proc.stdout
    else:
        proc = None
        stream = open(path, errors="replace")

    try:
        for line in stream:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            flag = int(fields[1])
            if flag & 4:
                continue
            yield fields[0], fields[2], int(fields[3])
    finally:
        stream.close()
        if proc is not None:
            stderr = proc.stderr.read()
            code = proc.wait()
            if code != 0:
                raise RuntimeError(f"samtools view failed for {path}:\n{stderr}")


def load_annotation(path):
    if not path or not Path(path).exists():
        return None

    by_chrom = defaultdict(list)
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            accession = (row.get("Accession") or "").split(".")[0]
            try:
                begin = int(row.get("Begin") or 0)
                end = int(row.get("End") or 0)
            except ValueError:
                continue
            by_chrom[accession].append({
                "symbol": row.get("Symbol") or "",
                "locus": row.get("Locus tag") or "",
                "begin": begin,
                "end": end,
            })
    return by_chrom


def chrom_key(chrom):
    return chrom.replace("ref|", "").replace("|", "")


def get_genes(annotation, chrom, left_bp, right_bp):
    if not annotation:
        return None, None, None

    genes = annotation.get(chrom_key(chrom), [])
    overlapping = [
        gene for gene in genes
        if gene["begin"] <= right_bp and gene["end"] >= left_bp
    ]
    candidate = None
    if overlapping:
        candidate = max(
            overlapping,
            key=lambda gene: min(gene["end"], right_bp) - max(gene["begin"], left_bp),
        )

    left_candidates = [gene for gene in genes if gene["end"] < left_bp]
    right_candidates = [gene for gene in genes if gene["begin"] > right_bp]
    left_gene = max(left_candidates, key=lambda gene: gene["end"]) if left_candidates else None
    right_gene = min(right_candidates, key=lambda gene: gene["begin"]) if right_candidates else None

    if left_gene:
        left_gene = {
            "symbol": left_gene["symbol"],
            "locus": left_gene["locus"],
            "end": left_gene["end"],
        }
    if right_gene:
        right_gene = {
            "symbol": right_gene["symbol"],
            "locus": right_gene["locus"],
            "begin": right_gene["begin"],
        }

    return candidate, left_gene, right_gene


def estimate_breakpoints(positions):
    best_count, best_start, left = 0, positions[0], 0
    for right_idx, pos in enumerate(positions):
        while positions[left] < pos - ZOOM_WINDOW:
            left += 1
        count = right_idx - left + 1
        if count > best_count:
            best_count = count
            best_start = positions[left]

    zone = [pos for pos in positions if best_start <= pos <= best_start + ZOOM_WINDOW]
    best_gap, left_bp, right_bp = 0, zone[0], zone[-1]
    for idx in range(1, len(zone)):
        gap = zone[idx] - zone[idx - 1]
        if gap > best_gap and gap >= MIN_GAP:
            best_gap = gap
            left_bp, right_bp = zone[idx - 1], zone[idx]

    return left_bp, right_bp, (left_bp + right_bp) // 2


def main():
    parser = argparse.ArgumentParser(
        description="Detect insertion junctions by crossing feature and reference alignments."
    )
    parser.add_argument("--refs", required=True, help="TSV manifest: reference_name<TAB>alignment.bam_or.sam")
    parser.add_argument("--features", required=True, help="TSV manifest: feature_name<TAB>alignment.bam_or.sam")
    parser.add_argument("--annotation", default="", help="Optional NCBI gene TSV annotation")
    parser.add_argument("--output", default="insertions.json")
    parser.add_argument("--min-junction-reads", type=int, default=MIN_JUNCTION_READS)
    args = parser.parse_args()

    ref_alignments = read_manifest(args.refs)
    feature_alignments = read_manifest(args.features)
    annotation = load_annotation(args.annotation)

    print(f"Reference alignments: {len(ref_alignments)}", file=sys.stderr)
    print(f"Feature alignments: {len(feature_alignments)}", file=sys.stderr)

    synthetic_reads = defaultdict(set)
    for feature_name, alignment in feature_alignments:
        count = 0
        for read_name, _contig, _pos in sam_records(alignment):
            synthetic_reads[read_name].add(feature_name)
            count += 1
        print(f"Feature {feature_name}: {count} mapped reads", file=sys.stderr)

    junction_data = defaultdict(lambda: defaultdict(list))
    for ref_name, alignment in ref_alignments:
        count = 0
        for read_name, chrom, pos in sam_records(alignment):
            if read_name not in synthetic_reads:
                continue
            for feature_name in synthetic_reads[read_name]:
                junction_data[feature_name][chrom].append(pos)
                count += 1
        print(f"Reference {ref_name}: {count} junction reads", file=sys.stderr)

    results = []
    for feature_name in sorted(junction_data):
        for chrom in sorted(junction_data[feature_name]):
            positions = sorted(junction_data[feature_name][chrom])
            if len(positions) < args.min_junction_reads:
                continue
            left_bp, right_bp, insertion_pos = estimate_breakpoints(positions)
            candidate, left_gene, right_gene = get_genes(annotation, chrom, left_bp, right_bp)
            results.append({
                "element": feature_name,
                "chromosome": chrom,
                "left_breakpoint": left_bp,
                "right_breakpoint": right_bp,
                "insertion_position": insertion_pos,
                "candidate_deleted_gene": candidate,
                "left_gene": left_gene,
                "right_gene": right_gene,
                "junction_positions": positions,
            })

    with open(args.output, "w") as handle:
        json.dump(results, handle, indent=2)
        handle.write("\n")

    print(f"Detected insertion signals: {len(results)}", file=sys.stderr)


if __name__ == "__main__":
    main()
