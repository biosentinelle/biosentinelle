#!/usr/bin/env python3
"""
Count panel k-mers in raw reads, alignment-free.

Reads FASTQ or FASTA (optionally gzipped), counts how many times each panel
k-mer occurs, and writes per-site reference and alternative allele counts.
Reads are canonicalised, so strand does not matter.

Throughput note
---------------
This is pure Python: expect roughly 1-3 Mbp of reads per second. That is fine
for a subsampled run and for the yeast proof of concept, and too slow for a
deep whole-genome sample scanned end to end. Two ways out, in order of
preference:

  * --max-bases  Subsample to a number of bases. Coverage, not read count, is
    what decides whether a site gets called, and read length varies by an order
    of magnitude between platforms: 400000 reads is 60 Mbp at 150 bp and 14 Mbp
    at 35 bp, which is 5x and 1.2x on a 12 Mb genome. Capping by reads silently
    produces an unusable scan on short-read data, so the cap is in bases.
    --max-reads still exists and still works, but prefer --max-bases.

  * --jellyfish  Hand the counting to jellyfish when it is installed. The panel
    k-mers are then queried against its database instead of being counted here.

A site is called only when its total depth reaches --min-depth; below that the
site is reported as not covered rather than as reference, because absence of a
k-mer is also what low coverage looks like.
"""

import argparse
import csv
import gzip
import os
import subprocess
import sys

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def canonical(kmer):
    reverse = revcomp(kmer)
    return kmer if kmer <= reverse else reverse


def opener(path):
    if path == "-":
        return sys.stdin
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def iter_sequences(paths, max_reads, max_bases=0):
    """Yield read sequences from FASTQ or FASTA inputs, stopping at either cap."""
    emitted = 0
    emitted_bases = 0
    for path in paths:
        with opener(path) as handle:
            first = handle.readline()
            if not first:
                continue
            is_fastq = first.startswith("@")
            handle.seek(0) if path != "-" else None
            if is_fastq:
                while True:
                    header = handle.readline()
                    if not header:
                        break
                    sequence = handle.readline().strip().upper()
                    handle.readline()
                    handle.readline()
                    if not sequence:
                        break
                    yield sequence
                    emitted += 1
                    emitted_bases += len(sequence)
                    if max_reads and emitted >= max_reads:
                        return
                    if max_bases and emitted_bases >= max_bases:
                        return
            else:
                chunks = []
                for line in handle:
                    line = line.strip()
                    if line.startswith(">"):
                        if chunks:
                            sequence = "".join(chunks).upper()
                            yield sequence
                            emitted += 1
                            emitted_bases += len(sequence)
                            if max_reads and emitted >= max_reads:
                                return
                            if max_bases and emitted_bases >= max_bases:
                                return
                        chunks = []
                    elif line:
                        chunks.append(line)
                if chunks:
                    yield "".join(chunks).upper()


def count_python(paths, wanted, k, max_reads, progress_every, max_bases=0):
    counts = dict.fromkeys(wanted, 0)
    bases = 0
    for index, sequence in enumerate(iter_sequences(paths, max_reads, max_bases), start=1):
        bases += len(sequence)
        # Reads with ambiguous bases are split: an N can never be part of a
        # panel k-mer, and skipping the whole read would bias against sites in
        # lower-quality regions.
        for piece in sequence.split("N"):
            for i in range(len(piece) - k + 1):
                key = canonical(piece[i:i + k])
                if key in counts:
                    counts[key] += 1
        if progress_every and index % progress_every == 0:
            sys.stderr.write("  %d reads, %.1f Mbp\n" % (index, bases / 1e6))
            sys.stderr.flush()
    return counts, bases


def count_jellyfish(paths, wanted, k, threads, hash_size, workdir):
    """Count with jellyfish, then query only the panel k-mers."""
    database = os.path.join(workdir, "reads.jf")
    command = ["jellyfish", "count", "-C", "-m", str(k), "-s", hash_size,
               "-t", str(threads), "-o", database] + list(paths)
    sys.stderr.write("running: %s\n" % " ".join(command))
    subprocess.run(command, check=True)

    query_input = os.path.join(workdir, "panel.kmers")
    with open(query_input, "w") as handle:
        for kmer in wanted:
            handle.write(kmer + "\n")
    result = subprocess.run(["jellyfish", "query", "-s", query_input, database],
                            check=True, capture_output=True, text=True)
    counts = dict.fromkeys(wanted, 0)
    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) == 2:
            counts[canonical(parts[0])] = int(parts[1])
    return counts, 0


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--kmers", required=True, help="panel k-mer lookup TSV")
    parser.add_argument("--reads", required=True, nargs="+", help="FASTQ or FASTA inputs")
    parser.add_argument("--out", required=True, help="output per-site allele counts TSV")
    parser.add_argument("--min-depth", type=int, default=3,
                        help="total allele count needed to call a site (default 3)")
    parser.add_argument("--min-allele-fraction", type=float, default=0.8,
                        help="fraction of reads needed to call a homozygous allele")
    parser.add_argument("--max-bases", type=int, default=0,
                        help="stop after this many bases, 0 to use all (the default). "
                             "Prefer this over --max-reads: it is coverage that "
                             "decides whether a site is called, and bases are the "
                             "unit of coverage. For a 12 Mb genome, 300000000 is "
                             "about 25x, which is ample.")
    parser.add_argument("--max-reads", type=int, default=0,
                        help="stop after this many reads, 0 to use all (the default). "
                             "Kept for convenience, but a read cap means different "
                             "coverage on different platforms")
    parser.add_argument("--jellyfish", action="store_true",
                        help="count with jellyfish instead of in Python")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--hash-size", default="100M")
    parser.add_argument("--workdir", default=".")
    parser.add_argument("--progress-every", type=int, default=100000)
    args = parser.parse_args()

    sites = {}
    wanted = set()
    with open(args.kmers, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            kmer = row["kmer"]
            wanted.add(kmer)
            sites.setdefault(row["site_id"], {})[row["allele"]] = kmer
    if not wanted:
        sys.exit("empty k-mer panel: %s" % args.kmers)
    k = len(next(iter(wanted)))
    if any(len(kmer) != k for kmer in wanted):
        sys.exit("panel mixes k-mer lengths, rebuild it with one k")

    sys.stderr.write("panel: %d sites, %d k-mers, k=%d\n" % (len(sites), len(wanted), k))
    if args.jellyfish:
        counts, bases = count_jellyfish(args.reads, wanted, k, args.threads,
                                        args.hash_size, args.workdir)
    else:
        counts, bases = count_python(args.reads, wanted, k, args.max_reads,
                                     args.progress_every, args.max_bases)
        sys.stderr.write("scanned %.1f Mbp\n" % (bases / 1e6))

    rows = []
    for site_id, alleles in sites.items():
        ref_count = counts.get(alleles.get("R", ""), 0)
        alt_count = counts.get(alleles.get("A", ""), 0)
        depth = ref_count + alt_count
        if depth < args.min_depth:
            call, fraction = "?", 0.0
        else:
            fraction = alt_count / float(depth)
            if fraction >= args.min_allele_fraction:
                call = "A"
            elif fraction <= 1.0 - args.min_allele_fraction:
                call = "R"
            else:
                # Both alleles well represented: a heterozygous diploid, a mixed
                # culture, or a cross. Reported as its own state, never forced
                # into one of the two homozygous calls.
                call = "H"
        rows.append({"site_id": site_id, "ref_count": ref_count, "alt_count": alt_count,
                     "depth": depth, "alt_fraction": round(fraction, 4), "call": call})

    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["site_id", "ref_count", "alt_count", "depth",
                                "alt_fraction", "call"], delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: row["site_id"]))

    tally = {}
    for row in rows:
        tally[row["call"]] = tally.get(row["call"], 0) + 1
    sys.stderr.write("calls: %s\n" % ", ".join(
        "%s=%d" % (call, tally[call]) for call in sorted(tally)))
    if tally.get("?", 0) > len(rows) * 0.5:
        sys.stderr.write(
            "warning: over half the panel is below --min-depth after %.1f Mbp. "
            "Either the scan was capped too low (raise --max-bases, 0 scans "
            "everything), the run is too shallow, or k is too close to the read "
            "length, which leaves too few k-mers per read\n" % (bases / 1e6))
    return 0


if __name__ == "__main__":
    sys.exit(main())
