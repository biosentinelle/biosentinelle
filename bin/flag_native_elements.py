#!/usr/bin/env python3
"""
Mark which tool elements are native to the host genome.

A cassette is assembled from parts, and some of those parts are the host's own
sequence: the GAL1 promoter, the ADH1 and CYC1 terminators and the TRP1 gene
are all S. cerevisiae sequence. Reads from a completely unmodified strain align
to them end to end at full identity, because the strain carries the native
locus. Others are foreign by construction, such as kanMX from a bacterial
transposon, NrsR, HygR, ble, GFP, or the Kluyveromyces lactis URA3 and TRP1
used precisely so that they do not recombine with the host copy.

Presence means two different things for the two groups:

  foreign   Covered at high identity means the sequence is in the sample and
            got there by engineering. This is real evidence.

  native    Covered at high identity means only that the strain has its own
            genome. It is not evidence of a cassette. Deciding that needs
            either reads spanning the junction between the element and its
            insertion site, or a copy number above one.

Without the distinction the report says "TRP1 present" for a wild-type strain
and anyone reading it concludes there is a TRP1 cassette. Measured on
SRR26617475, the only elements scoring above 90 were TRP1, the TRP1 and GAL1
promoters and the ADH1 and CYC1 terminators: every one of them native, and the
strain carries no cassette at all.

The split is decided empirically rather than by hand, by asking what fraction
of each element's k-mers occur in the reference genome, so it stays correct
when elements are added and transfers to any organism with a reference.
"""

import argparse
import csv
import os
import sys

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def canonical(kmer):
    reverse = revcomp(kmer)
    return kmer if kmer <= reverse else reverse


def read_fasta(path):
    name, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.upper())
    if name is not None:
        yield name, "".join(chunks)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tools", required=True, help="tools.fasta")
    parser.add_argument("--manifest", required=True, help="tools.tsv to annotate in place")
    parser.add_argument("--reference", required=True, help="host reference genome FASTA")
    parser.add_argument("-k", type=int, default=21)
    parser.add_argument("--native-threshold", type=float, default=0.50,
                        help="fraction of an element's k-mers that must occur in the "
                             "host genome for it to count as native (default 0.50). "
                             "The split is strongly bimodal, so the threshold is not "
                             "delicate: on the yeast set at k=21 the native elements "
                             "score 0.71 to 1.00 and every foreign one scores 0.00. "
                             "k matters more than the threshold. Exact k-mers are "
                             "unforgiving of divergence, and a native element that "
                             "differs from the host by 1.5 per cent keeps only "
                             "0.985^k of its k-mers, so a long k reads a diverged "
                             "native element as foreign. Below about k=15 chance "
                             "matches start appearing instead.")
    args = parser.parse_args()

    elements = list(read_fasta(args.tools))
    if not elements:
        sys.exit("no sequence in %s" % args.tools)

    # All element k-mers at once, then one pass over the genome.
    wanted = set()
    element_kmers = {}
    for name, sequence in elements:
        kmers = {canonical(sequence[i:i + args.k])
                 for i in range(len(sequence) - args.k + 1)}
        element_kmers[name] = kmers
        wanted |= kmers
    if not wanted:
        sys.exit("every element is shorter than k=%d" % args.k)

    found = set()
    for _, sequence in read_fasta(args.reference):
        for i in range(len(sequence) - args.k + 1):
            key = canonical(sequence[i:i + args.k])
            if key in wanted:
                found.add(key)

    fractions = {}
    for name, kmers in element_kmers.items():
        fractions[name] = (len(kmers & found) / float(len(kmers))) if kmers else 0.0

    rows = []
    with open(args.manifest, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = list(reader.fieldnames)
        rows = list(reader)
    for column in ("host_kmer_fraction", "origin"):
        if column not in columns:
            columns.append(column)

    native = foreign = unknown = 0
    for row in rows:
        element = row["element_id"]
        if element not in fractions:
            # Shorter than k, so this test cannot be applied to it.
            row["host_kmer_fraction"] = ""
            row["origin"] = "unknown"
            unknown += 1
            continue
        fraction = fractions[element]
        row["host_kmer_fraction"] = "%.4f" % fraction
        if fraction >= args.native_threshold:
            row["origin"] = "native"
            native += 1
        else:
            row["origin"] = "foreign"
            foreign += 1

    with open(args.manifest, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print("annotated %s: %d native, %d foreign, %d too short to test (k=%d)"
          % (args.manifest, native, foreign, unknown, args.k))
    for row in sorted(rows, key=lambda r: r.get("host_kmer_fraction") or "", reverse=True):
        if row["origin"] == "native":
            print("  native  %-24s %s" % (row["name"], row["host_kmer_fraction"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
