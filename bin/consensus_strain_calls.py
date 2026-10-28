#!/usr/bin/env python3
"""
Combine the calls of several independent read runs of one strain into the
strain's own variants and callable regions.

Two stocks of the same strain sequenced in two laboratories share the strain's
variants, and each carries its own: mutations picked up since the stocks were
separated, and sequencing or mapping artefacts. Only the first kind identifies
the strain. A variant is therefore kept only when every run is callable there
and every run calls the same genotype. Where the runs disagree, or one of them
cannot see the position, the strain is recorded as unknown there, never as
reference: claiming the reference allele would be claiming knowledge the runs
do not give.

For a haploid strain, a heterozygous call is not biology. It is reads of two
paralogous copies landing on one, or a duplication, so the position is made
unknown too.

With a single run the run's calls are taken as they are, and the report says
the strain rests on one stock.
"""

import argparse
import gzip
import sys


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_calls(path, min_qual, min_depth):
    calls, het = {}, set()
    with opener(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, _, ref, alt, qual = fields[:6]
            keys = fields[8].split(":")
            values = fields[9].split(":")
            sample = dict(zip(keys, values))
            gt = sample.get("GT", "./.").replace("|", "/")
            if "." in gt or gt in ("0/0", "0"):
                continue
            try:
                if float(qual) < min_qual:
                    continue
            except ValueError:
                pass
            depth = sample.get("DP", "")
            if depth.isdigit() and int(depth) < min_depth:
                continue
            alleles = gt.split("/")
            key = (chrom, int(pos), ref, alt)
            if len(set(alleles)) > 1:
                het.add((chrom, int(pos)))
                calls[key] = "0/1"
            else:
                calls[key] = "1/1"
    return calls, het


def load_bed(path):
    regions = {}
    with opener(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 3 or line.startswith(("#", "track")):
                continue
            regions.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    for chrom in regions:
        regions[chrom].sort()
    return regions


def intersect(a, b):
    out = {}
    for chrom in set(a) & set(b):
        x, y, i, j, merged = a[chrom], b[chrom], 0, 0, []
        while i < len(x) and j < len(y):
            start, end = max(x[i][0], y[j][0]), min(x[i][1], y[j][1])
            if start < end:
                merged.append((start, end))
            if x[i][1] < y[j][1]:
                i += 1
            else:
                j += 1
        out[chrom] = merged
    return out


def subtract_points(regions, points):
    """Remove single 1-based positions from the intervals."""
    by_chrom = {}
    for chrom, pos in points:
        by_chrom.setdefault(chrom, []).append(pos - 1)
    for chrom, holes in by_chrom.items():
        holes.sort()
        out = []
        for start, end in regions.get(chrom, []):
            cursor = start
            for hole in holes:
                if hole < cursor or hole >= end:
                    continue
                if hole > cursor:
                    out.append((cursor, hole))
                cursor = hole + 1
            if cursor < end:
                out.append((cursor, end))
        regions[chrom] = out
    return regions


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--strain", required=True)
    parser.add_argument("--run", action="append", required=True, metavar="RUN=VCF=BED")
    parser.add_argument("--ploidy", type=int, default=1, choices=(1, 2))
    parser.add_argument("--min-qual", type=float, default=30.0)
    parser.add_argument("--min-depth", type=int, default=10)
    parser.add_argument("--out-vcf", required=True)
    parser.add_argument("--out-bed", required=True)
    args = parser.parse_args()

    runs = []
    for item in args.run:
        name, vcf, bed = item.split("=")
        calls, het = load_calls(vcf, args.min_qual, args.min_depth)
        runs.append((name, calls, het, load_bed(bed)))

    callable_regions = runs[0][3]
    for _, _, _, regions in runs[1:]:
        callable_regions = intersect(callable_regions, regions)

    keys = set()
    for _, calls, _, _ in runs:
        keys.update(calls)
    kept, holes, disagreements = {}, set(), 0
    for key in keys:
        genotypes = {calls.get(key) for _, calls, _, _ in runs}
        if len(genotypes) == 1 and None not in genotypes:
            gt = genotypes.pop()
            if gt == "0/1" and args.ploidy == 1:
                holes.add(key[:2])
                continue
            kept[key] = gt
        else:
            disagreements += 1
            for offset in range(max(len(key[2]), 1)):
                holes.add((key[0], key[1] + offset))
    callable_regions = subtract_points(callable_regions, holes)

    with open(args.out_vcf, "w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##source=consensus_strain_calls.py runs=%s\n"
                     % ",".join(name for name, _, _, _ in runs))
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % args.strain)
        for chrom, pos, ref, alt in sorted(kept):
            handle.write("%s\t%d\t.\t%s\t%s\t60\tPASS\t.\tGT\t%s\n"
                         % (chrom, pos, ref, alt, kept[(chrom, pos, ref, alt)]))
    with open(args.out_bed, "w") as handle:
        for chrom in sorted(callable_regions):
            for start, end in callable_regions[chrom]:
                handle.write("%s\t%d\t%d\n" % (chrom, start, end))

    span = sum(end - start for regions in callable_regions.values() for start, end in regions)
    print("%s: %d run(s), %d variants kept, %d dropped for disagreement between runs, "
          "%d positions made unknown, %.2f Mb callable"
          % (args.strain, len(runs), len(kept), disagreements, len(holes), span / 1e6))
    if len(runs) == 1:
        print("warning: %s rests on one run; its private variants include that stock's "
              "own mutations and cannot be told apart from the strain's" % args.strain)
    return 0


if __name__ == "__main__":
    sys.exit(main())
