#!/usr/bin/env python3
"""
Build the k-mers that read the mating type of a sample.

Every S. cerevisiae genome carries both mating-type cassettes: the silent
HMLalpha and HMRa copies, plus the expressed MAT locus, which is either a or
alpha. So the presence of a or alpha sequence says nothing; their ratio does.
Counting one copy per cassette:

    MATa haploid      Ya at MAT and HMR, Yalpha at HML      a : alpha = 2 : 1
    MATalpha haploid  Ya at HMR, Yalpha at MAT and HML      a : alpha = 1 : 2
    MATa/MATalpha     the same summed over two homologs     a : alpha = 3 : 3

That is what separates FY4 (MATa) from S288C (MATalpha), and BY4741 (MATa)
from BY4742 (MATalpha), which no SNP does.

The k-mers come from the a1 and a2 genes at HMR and the alpha1 and alpha2 genes
at MAT, with coordinates fetched from Ensembl by systematic name. The X region
is shared by all three cassettes and carries no information, so a k-mer is
kept only when it is specific: an a k-mer must occur once in the reference
(at HMR, the reference being MATalpha) and an alpha k-mer exactly twice (MAT
and HML), and neither may occur in the other set.
"""

import argparse
import csv
import os
import sys

from build_marker_loci import ENSEMBL, chromosome_names, fetch_json
from marker_effects import revcomp

GENES = {"YCR097W": "a", "YCR096C": "a",            # HMRA1, HMRA2
         "YCR040W": "alpha", "YCR039C": "alpha"}    # MATALPHA1, MATALPHA2
EXPECTED_COPIES = {"a": 1, "alpha": 2}


def canonical(kmer):
    reverse = revcomp(kmer)
    return kmer if kmer <= reverse else reverse


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--species", default="saccharomyces_cerevisiae")
    parser.add_argument("--assembly", default="R64-1-1")
    parser.add_argument("-k", type=int, default=21)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    names, sequences = chromosome_names(args.reference)
    candidates = {"a": set(), "alpha": set()}
    sources = []
    for gene, cassette in GENES.items():
        info = fetch_json("%s/lookup/id/%s?species=%s" % (ENSEMBL, gene, args.species))
        if info.get("assembly_name") != args.assembly:
            sys.exit("Ensembl reports %s for %s, the reference is %s"
                     % (info.get("assembly_name"), gene, args.assembly))
        record = names.get(info["seq_region_name"])
        sequence = sequences[record][info["start"] - 1:info["end"]]
        for i in range(len(sequence) - args.k + 1):
            candidates[cassette].add(canonical(sequence[i:i + args.k]))
        sources.append("%s %s:%d-%d" % (gene, info["seq_region_name"], info["start"],
                                        info["end"]))

    wanted = candidates["a"] | candidates["alpha"]
    counts = dict.fromkeys(wanted, 0)
    for sequence in sequences.values():
        for i in range(len(sequence) - args.k + 1):
            key = canonical(sequence[i:i + args.k])
            if key in counts:
                counts[key] += 1

    rows = []
    for cassette, kmers in candidates.items():
        other = candidates["alpha" if cassette == "a" else "a"]
        for kmer in sorted(kmers):
            if kmer in other or counts[kmer] != EXPECTED_COPIES[cassette]:
                continue
            rows.append((kmer, "MAT" + cassette))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["kmer", "class"])
        writer.writerows(rows)
    print("mating-type k-mers (k=%d) from %s" % (args.k, "; ".join(sources)))
    for cassette in ("a", "alpha"):
        print("  MAT%-6s %d specific k-mers" % (
            cassette, sum(1 for _, c in rows if c == "MAT" + cassette)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
