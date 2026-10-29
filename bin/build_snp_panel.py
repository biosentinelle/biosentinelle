#!/usr/bin/env python3
"""
Turn a strain-versus-reference VCF into a k-mer SNP panel.

For every biallelic SNV that passes the filters, two k-mers are emitted: the one
carrying the reference allele and the one carrying the strain allele, each
centred on the variant position. Detection then becomes exact k-mer counting in
the raw reads, with no alignment and no variant calling on the query sample.

Three filters decide whether a site is usable, and all three matter:

  flank purity   The k-mer is only a valid probe if the k-1 bases around the
                 SNP are identical between the reference and the strain. If a
                 second variant sits 4 bp away, the strain's real reads carry a
                 k-mer that matches neither probe, and the site silently reads
                 as absent. Sites with another variant inside the window are
                 dropped.

  uniqueness     A k-mer that occurs more than once in the reference genome is
                 not evidence about one locus. Yeast has Ty elements, tandem
                 rDNA and subtelomeric repeats, so this is not a rare case.
                 Both probes must be unique in the reference.

  composition    Homopolymers and very low-complexity k-mers accumulate
                 sequencing errors and cross-match; they are dropped.
"""

import argparse
import csv
import gzip
import os
import sys

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def canonical(kmer):
    """Strand-independent form, so a read from either strand gives the same key."""
    reverse = revcomp(kmer)
    return kmer if kmer <= reverse else reverse


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta_dict(path):
    sequences, name, chunks = {}, None, []
    with opener(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.upper())
    if name is not None:
        sequences[name] = "".join(chunks)
    return sequences


def parse_vcf(path, min_qual):
    """Yield (chrom, pos_1based, ref, alt) for biallelic SNVs passing quality."""
    with opener(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            chrom, pos, _, ref, alt, qual = fields[:6]
            if "," in alt:                      # multiallelic, ambiguous expectation
                continue
            if len(ref) != 1 or len(alt) != 1:  # indel, breaks a fixed-length k-mer
                continue
            if ref not in "ACGT" or alt not in "ACGT":
                continue
            try:
                if qual != "." and float(qual) < min_qual:
                    continue
            except ValueError:
                pass
            yield chrom, int(pos), ref, alt


def low_complexity(kmer, max_homopolymer, min_distinct):
    if len(set(kmer)) < min_distinct:
        return True
    run, longest, previous = 1, 1, ""
    for base in kmer:
        run = run + 1 if base == previous else 1
        longest = max(longest, run)
        previous = base
    return longest > max_homopolymer


def count_in_genome(sequences, wanted):
    """Count occurrences of each wanted canonical k-mer across both strands."""
    if not wanted:
        return {}
    k = len(next(iter(wanted)))
    counts = dict.fromkeys(wanted, 0)
    for sequence in sequences.values():
        for i in range(len(sequence) - k + 1):
            key = canonical(sequence[i:i + k])
            if key in counts:
                counts[key] += 1
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vcf", required=True, help="strain vs reference VCF")
    parser.add_argument("--reference", required=True, help="reference genome FASTA")
    parser.add_argument("--strain", required=True, help="strain identifier")
    parser.add_argument("--out", required=True, help="output panel TSV")
    parser.add_argument("--out-positions",
                        help="companion file listing every variant position of this "
                             "strain, kept or rejected. merge_snp_panels.py needs it "
                             "to reject sites whose k-mer window is broken in a "
                             "different panel strain.")
    parser.add_argument("-k", type=int, default=31, help="k-mer length, odd (default 31)")
    parser.add_argument("--min-qual", type=float, default=30.0)
    parser.add_argument("--max-homopolymer", type=int, default=6)
    parser.add_argument("--min-distinct-bases", type=int, default=3)
    parser.add_argument("--max-sites", type=int, default=0,
                        help="cap the panel size, 0 for no cap; sites are thinned evenly")
    args = parser.parse_args()

    if args.k % 2 == 0:
        parser.error("k must be odd so the variant sits in the middle of the k-mer")
    flank = args.k // 2

    reference = read_fasta_dict(args.reference)
    variants = list(parse_vcf(args.vcf, args.min_qual))
    if not variants:
        sys.exit("no usable SNV in %s" % args.vcf)

    # Positions of every variant, used to reject sites whose window is not clean.
    occupied = {}
    for chrom, pos, _, _ in variants:
        occupied.setdefault(chrom, set()).add(pos)

    candidates, rejected = [], {"no_contig": 0, "edge": 0, "ref_mismatch": 0,
                                "dirty_flank": 0, "low_complexity": 0, "repeat": 0}
    for chrom, pos, ref, alt in variants:
        sequence = reference.get(chrom)
        if sequence is None:
            rejected["no_contig"] += 1
            continue
        start, end = pos - 1 - flank, pos + flank
        if start < 0 or end > len(sequence):
            rejected["edge"] += 1
            continue
        window = sequence[start:end]
        if window[flank] != ref:
            # The VCF and the FASTA disagree; usually a contig naming or build
            # mismatch. Dropping is right: we cannot trust the coordinate.
            rejected["ref_mismatch"] += 1
            continue
        neighbours = occupied.get(chrom, ())
        if any(p in neighbours for p in range(pos - flank, pos + flank + 1) if p != pos):
            rejected["dirty_flank"] += 1
            continue
        if "N" in window:
            rejected["low_complexity"] += 1
            continue
        ref_kmer = window
        alt_kmer = window[:flank] + alt + window[flank + 1:]
        if (low_complexity(ref_kmer, args.max_homopolymer, args.min_distinct_bases) or
                low_complexity(alt_kmer, args.max_homopolymer, args.min_distinct_bases)):
            rejected["low_complexity"] += 1
            continue
        candidates.append((chrom, pos, ref, alt, canonical(ref_kmer), canonical(alt_kmer)))

    # Uniqueness is checked once, in a single pass over the genome, for the
    # reference probes of all surviving candidates. The alternative probe does
    # not exist in the reference by construction, but it can still collide with
    # an unrelated locus, so it is checked in the same pass.
    wanted = set()
    for _, _, _, _, ref_kmer, alt_kmer in candidates:
        wanted.add(ref_kmer)
        wanted.add(alt_kmer)
    counts = count_in_genome(reference, wanted)

    sites = []
    for chrom, pos, ref, alt, ref_kmer, alt_kmer in candidates:
        if counts.get(ref_kmer, 0) != 1 or counts.get(alt_kmer, 0) != 0:
            rejected["repeat"] += 1
            continue
        sites.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                      "ref_kmer": ref_kmer, "alt_kmer": alt_kmer, "strain": args.strain})

    if args.max_sites and len(sites) > args.max_sites:
        step = len(sites) / float(args.max_sites)
        sites = [sites[int(i * step)] for i in range(args.max_sites)]

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    columns = ["strain", "chrom", "pos", "ref", "alt", "ref_kmer", "alt_kmer"]
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for site in sites:
            writer.writerow(site)

    # Every variant position, not only the kept ones. A site rejected here can
    # still be another strain's panel site, and that strain's probe is only
    # valid if this strain's window is clean too.
    if args.out_positions:
        with open(args.out_positions, "w", newline="") as handle:
            handle.write("chrom\tpos\n")
            for chrom, pos, _, _ in variants:
                handle.write("%s\t%d\n" % (chrom, pos))

    print("%s: %d variants in, %d sites kept (k=%d)" %
          (args.strain, len(variants), len(sites), args.k))
    for reason in sorted(rejected):
        if rejected[reason]:
            print("  rejected %-16s %d" % (reason, rejected[reason]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
