#!/usr/bin/env python3
"""
Turn read alignments against the tool element database into per-element scores.

Two independent measurements, both needed:

  breadth    fraction of the element's bases covered at or above --min-depth.
             Answers "is the whole element there". Partial breadth is the
             signature of a truncated cassette, a shared promoter or terminator
             matching only its conserved part, or a chance match.

  identity   fraction of aligned bases that match the element exactly, taken
             from the NM tag over the aligned length. Answers "is it the same
             sequence". An element can be covered end to end by reads from a
             diverged homologue: yeast's own TRP1 will pull reads onto the TRP1
             in a cassette. Identity is what separates the two.

Input is a `samtools depth -a` file and a `samtools view` stream reduced to
reference name, CIGAR and NM tag. Both are produced by the ToolScan module.
"""

import argparse
import csv
import re
import sys
from collections import defaultdict

CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")


def aligned_length(cigar):
    """Bases of the read that take part in the alignment, excluding clipping."""
    if cigar == "*":
        return 0
    return sum(int(length) for length, operation in CIGAR_PATTERN.findall(cigar)
               if operation in "MID=X")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--depth", required=True,
                        help="samtools depth -a output: reference, position, depth")
    parser.add_argument("--alignments", required=True,
                        help="reference, CIGAR, NM per aligned read, tab separated")
    parser.add_argument("--manifest", required=True, help="tools.tsv")
    parser.add_argument("--out", required=True)
    parser.add_argument("--min-depth", type=int, default=1,
                        help="depth at which a base counts as covered (default 1)")
    args = parser.parse_args()

    lengths, order = {}, []
    with open(args.manifest, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            lengths[row["element_id"]] = int(row["length_bp"])
            order.append(row["element_id"])

    covered = defaultdict(int)
    depth_sum = defaultdict(int)
    with open(args.depth) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 3:
                continue
            element, depth = fields[0], int(fields[2])
            depth_sum[element] += depth
            if depth >= args.min_depth:
                covered[element] += 1

    matched = defaultdict(int)
    aligned = defaultdict(int)
    reads = defaultdict(int)
    with open(args.alignments) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            element, cigar, nm = fields[0], fields[1], fields[2]
            if element == "*":
                continue
            span = aligned_length(cigar)
            if not span:
                continue
            try:
                mismatches = int(nm)
            except ValueError:
                # No NM tag on this record. Counting it as a perfect match would
                # inflate identity, so the read is left out of the identity
                # calculation entirely and only contributes to breadth.
                continue
            aligned[element] += span
            matched[element] += max(span - mismatches, 0)
            reads[element] += 1

    rows = []
    for element in order:
        length = lengths[element]
        breadth = covered.get(element, 0) / float(length) if length else 0.0
        identity = (matched.get(element, 0) / float(aligned[element])
                    if aligned.get(element) else 0.0)
        mean_depth = depth_sum.get(element, 0) / float(length) if length else 0.0
        rows.append({
            "element_id": element,
            "length_bp": length,
            "covered_bases": covered.get(element, 0),
            "breadth": round(breadth, 4),
            "identity": round(identity, 4),
            "mean_depth": round(mean_depth, 2),
            "reads": reads.get(element, 0),
        })

    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t",
            fieldnames=["element_id", "length_bp", "covered_bases", "breadth",
                        "identity", "mean_depth", "reads"])
        writer.writeheader()
        writer.writerows(rows)

    detected = sum(1 for row in rows if row["breadth"] >= 0.99 and row["identity"] >= 0.99)
    print("scored %d elements, %d fully covered at high identity" % (len(rows), detected))
    return 0


if __name__ == "__main__":
    sys.exit(main())
