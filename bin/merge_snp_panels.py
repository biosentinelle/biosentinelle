#!/usr/bin/env python3
"""
Merge per-strain k-mer panels into one discriminant site matrix.

Each input panel lists the sites where one strain differs from the reference.
Merging them gives, for every site, the allele expected in every panel strain,
which is what actually identifies a background: not "this strain has a variant
here" but "this combination of alleles across these sites is this strain".

Expected allele per strain at a site:
  A   the strain carries the alternative allele (the site is in its panel)
  R   the strain carries the reference allele
  ?   unknown, the strain's assembly did not cover the site

Without a callable-region file a strain absent from a site is recorded as R,
and `assumed_ref` is set on that site so the report can say how much of the
evidence rests on that assumption. Alignment gaps in subtelomeres and repeats
are common enough that silently treating them as reference would overstate
confidence.

Sites that carry no information are dropped: if every strain in the panel has
the same expected allele, the site cannot separate any pair of them.
"""

import argparse
import csv
import os
import sys
from collections import OrderedDict


def load_panel(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_callable(path):
    """Read a BED of callable regions into per-contig sorted interval lists."""
    regions = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            regions.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    for contig in regions:
        regions[contig].sort()
    return regions


def is_callable(regions, chrom, pos):
    """pos is 1-based; BED is half-open 0-based."""
    import bisect
    intervals = regions.get(chrom)
    if not intervals:
        return False
    index = bisect.bisect_right(intervals, (pos - 1, float("inf"))) - 1
    if index < 0:
        return False
    start, end = intervals[index]
    return start <= pos - 1 < end


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--panel", action="append", required=True, metavar="STRAIN=FILE",
                        help="per-strain panel, repeatable, e.g. --panel W303=w303.snp.tsv")
    parser.add_argument("--callable", action="append", default=[], metavar="STRAIN=BED",
                        help="optional callable regions per strain, repeatable")
    parser.add_argument("--positions", action="append", default=[], metavar="STRAIN=FILE",
                        help="per-strain variant position files from build_snp_panel.py "
                             "--out-positions, repeatable. Used to drop sites whose "
                             "k-mer window is broken in another panel strain.")
    parser.add_argument("-k", type=int, default=31,
                        help="k-mer length the panels were built with, used for the "
                             "cross-strain window check (default 31)")
    parser.add_argument("--reference-strain", metavar="NAME",
                        help="add the reference strain as a panel member expecting the "
                             "reference allele everywhere. Use this instead of giving "
                             "the reference its own assembly panel: aligning an "
                             "assembly of the reference strain back to the reference "
                             "yields assembly and version differences, not biology, "
                             "and those would be scored as real mismatches.")
    parser.add_argument("--out-sites", required=True, help="output site matrix TSV")
    parser.add_argument("--out-kmers", required=True, help="output k-mer lookup TSV")
    parser.add_argument("--keep-uninformative", action="store_true",
                        help="keep sites where every strain expects the same allele")
    args = parser.parse_args()

    panels, callable_regions = OrderedDict(), {}
    for item in args.panel:
        if "=" not in item:
            parser.error("--panel expects STRAIN=FILE, got %r" % item)
        strain, path = item.split("=", 1)
        panels[strain] = load_panel(path)
    for item in args.callable:
        strain, path = item.split("=", 1)
        callable_regions[strain] = load_callable(path)

    # Variant positions of every strain, pooled. A probe k-mer has to be valid
    # for all of them: the sample is one unknown strain, and whichever it is,
    # its reads must carry either the reference probe or the alternative probe
    # at that site. If a second variant of some strain sits inside the window,
    # that strain's reads match neither, and the site reads as uncovered for
    # everyone. Rejecting those sites costs panel size and buys real power.
    all_positions = {}
    for item in args.positions:
        strain, path = item.split("=", 1)
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                all_positions.setdefault(row["chrom"], set()).add(int(row["pos"]))

    strains = list(panels)
    if args.reference_strain:
        if args.reference_strain in strains:
            parser.error("--reference-strain %s also has its own panel; pick one"
                         % args.reference_strain)
        strains.append(args.reference_strain)

    # Collect every site seen in any panel. Two strains can carry a different
    # alternative allele at the same position; that is informative, so the
    # alternative allele is part of the site key.
    sites = OrderedDict()
    for strain, rows in panels.items():
        for row in rows:
            key = (row["chrom"], int(row["pos"]), row["ref"], row["alt"])
            site = sites.setdefault(key, {
                "chrom": row["chrom"], "pos": int(row["pos"]),
                "ref": row["ref"], "alt": row["alt"],
                "ref_kmer": row["ref_kmer"], "alt_kmer": row["alt_kmer"],
                "carriers": set(),
            })
            site["carriers"].add(strain)
            # Panels built against the same reference must agree on the probes.
            if site["ref_kmer"] != row["ref_kmer"]:
                sys.exit("k-mer conflict at %s:%s between panels, were they built "
                         "against the same reference and the same k?" %
                         (row["chrom"], row["pos"]))

    flank = args.k // 2
    rows_out, kmer_rows = [], []
    dropped_uninformative = dropped_dirty = 0
    for index, (key, site) in enumerate(sites.items()):
        neighbours = all_positions.get(site["chrom"], ())
        if neighbours and any(pos in neighbours
                              for pos in range(site["pos"] - flank, site["pos"] + flank + 1)
                              if pos != site["pos"]):
            dropped_dirty += 1
            continue

        expected = {}
        for strain in strains:
            if strain in site["carriers"]:
                expected[strain] = "A"
            elif strain in callable_regions:
                expected[strain] = "R" if is_callable(
                    callable_regions[strain], site["chrom"], site["pos"]) else "?"
            else:
                expected[strain] = "R"
        known = [value for value in expected.values() if value != "?"]
        if not args.keep_uninformative and len(set(known)) < 2:
            dropped_uninformative += 1
            continue

        site_id = "s%06d" % index
        # Weight by how evenly the site splits the panel: a site carried by one
        # strain out of five separates that strain from four others and is worth
        # more than one carried by every strain but one in a different sense.
        # Both extremes are informative, a 50/50 split is the most informative.
        carriers = sum(1 for value in known if value == "A")
        total = len(known)
        balance = (carriers * (total - carriers)) / float(total * total / 4.0) if total else 0.0
        row = {"site_id": site_id, "chrom": site["chrom"], "pos": site["pos"],
               "ref": site["ref"], "alt": site["alt"],
               "n_callable_strains": total,
               # The reference strain is excluded: it carries the reference
               # allele everywhere by definition, not by assumption, so
               # counting it here would flag every site and make the warning
               # meaningless.
               "assumed_ref": "yes" if any(
                   strain not in site["carriers"]
                   and strain not in callable_regions
                   and strain != args.reference_strain
                   for strain in strains) else "no",
               "weight": round(max(balance, 0.1), 4)}
        for strain in strains:
            row["exp_" + strain] = expected[strain]
        rows_out.append(row)
        kmer_rows.append({"kmer": site["ref_kmer"], "site_id": site_id, "allele": "R"})
        kmer_rows.append({"kmer": site["alt_kmer"], "site_id": site_id, "allele": "A"})

    # A k-mer that ended up attached to two different sites is ambiguous
    # evidence and is removed from both, rather than counted twice.
    seen, duplicates = {}, set()
    for row in kmer_rows:
        if row["kmer"] in seen and seen[row["kmer"]] != row["site_id"]:
            duplicates.add(row["kmer"])
        seen[row["kmer"]] = row["site_id"]
    if duplicates:
        bad_sites = {row["site_id"] for row in kmer_rows if row["kmer"] in duplicates}
        kmer_rows = [row for row in kmer_rows if row["site_id"] not in bad_sites]
        rows_out = [row for row in rows_out if row["site_id"] not in bad_sites]

    for path in (args.out_sites, args.out_kmers):
        directory = os.path.dirname(os.path.abspath(path))
        os.makedirs(directory, exist_ok=True)

    site_columns = (["site_id", "chrom", "pos", "ref", "alt", "weight",
                     "n_callable_strains", "assumed_ref"] +
                    ["exp_" + strain for strain in strains])
    with open(args.out_sites, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=site_columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_out)

    with open(args.out_kmers, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["kmer", "site_id", "allele"],
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(kmer_rows)

    print("strains: %s" % ", ".join(strains))
    print("informative sites: %d  (k-mers: %d)" % (len(rows_out), len(kmer_rows)))
    if dropped_uninformative:
        print("dropped %d sites that no pair of strains disagrees on" % dropped_uninformative)
    if dropped_dirty:
        print("dropped %d sites whose k-mer window carries another strain's variant"
              % dropped_dirty)
    elif not args.positions:
        print("note: no --positions given, so sites broken by a neighbouring variant "
              "in another panel strain were not removed; they will read as uncovered")
    if duplicates:
        print("dropped %d sites whose k-mers collided with another site" % len(bad_sites))
    assumed = sum(1 for row in rows_out if row["assumed_ref"] == "yes")
    if assumed:
        print("warning: %d sites assume reference allele for at least one strain "
              "(no callable regions given)" % assumed)
    return 0


if __name__ == "__main__":
    sys.exit(main())
