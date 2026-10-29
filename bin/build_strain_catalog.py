#!/usr/bin/env python3
"""
Build the strain catalogue: every strain the panel describes, once, with the
source its SNPs are taken from.

Three kinds of source, in order of preference:

  reads     Illumina runs aligned to the reference and called here. Used for
            the S288C family, where the real differences number in the tens
            and the base errors of a short-read assembly would be as many.
            Several independent runs per strain are combined, and a variant
            is kept only when every run agrees on it.
  matrix    The 1011 Yeast Genomes population matrix (Peter et al. 2018),
            genotypes called from reads at every position variable in the
            population. Preferred to an assembly of the same isolate: one
            caller, one reference, and explicit no-calls.
  assembly  A whole-genome assembly, aligned to the reference. ScRAPdb
            long-read assemblies first, then SGD's.

The same strain often appears in several sources under slightly different
names (RM11-1a, RM11-1A, RM, AAA). Names are normalised before comparison,
and aliases.tsv holds the matches normalisation cannot make. A strain that
survives under two names would be scored as two strains that can never be
told apart, which is a wrong answer, not a cautious one.
"""

import argparse
import csv
import os
import re
import sys


def normalise(name):
    """Upper case, alphanumerics only, without the _1b / _b sub-clone suffix.

    A separator between two digits is kept: K1_1 (Kyokai no. 1) and K11 are
    different sake strains, and dropping it would make them one.
    """
    name = re.sub(r"[_-](1[a-d]|[a-d])$", "", name.strip(), flags=re.IGNORECASE).upper()
    name = re.sub(r"(?<=[0-9])[^A-Z0-9]+(?=[0-9])", "-", name)
    return re.sub(r"[^A-Z0-9-]", "", name)


def safe_id(name):
    return re.sub(r"[^A-Za-z0-9.\-]+", "_", name.strip()).strip("_")


def load_aliases(path):
    """alias -> canonical name, both normalised."""
    aliases = {}
    if path:
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                aliases[normalise(row["alias"])] = normalise(row["strain"])
    return aliases


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--matrix-table", required=True,
                        help="Table S1 of Peter et al. 2018 exported as TSV")
    parser.add_argument("--matrix-samples", required=True,
                        help="sample names of the matrix, one per line")
    parser.add_argument("--reads", required=True,
                        help="TSV strain, run, role, study of the read runs")
    parser.add_argument("--reads-dir", required=True,
                        help="directory holding <run>_1.fastq.gz / <run>.fastq.gz")
    parser.add_argument("--scrapdb-dir", required=True)
    parser.add_argument("--sgd-dir", required=True)
    parser.add_argument("--exclude", default="",
                        help="comma separated assemblies to leave out (other species, "
                             "mislabelled). A name prefixed scrapdb: or sgd: is left out "
                             "of that collection only")
    parser.add_argument("--aliases", help="TSV alias, strain")
    parser.add_argument("--path-map", default="", metavar="LOCAL=TARGET",
                        help="rewrite the leading LOCAL of written paths as TARGET, so "
                             "the catalogue points at the paths the build will see")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    aliases = load_aliases(args.aliases)
    excluded = set()
    for item in (x.strip() for x in args.exclude.split(",") if x.strip()):
        collection, _, name = item.rpartition(":")
        for scope in ([collection] if collection else ["scrapdb", "sgd"]):
            excluded.add((scope, normalise(name)))

    def canonical(name):
        key = normalise(name)
        return aliases.get(key, key)

    def written(path, root):
        path = path.replace("\\", "/")
        if args.path_map:
            local, target = args.path_map.split("=", 1)
            if path.startswith(local):
                path = target + path[len(local):]
        return path

    rows, taken = [], {}

    def add(row, names):
        # Two matrix samples sharing a name are kept apart: the 1011 set holds
        # some isolates twice, as the original and as a spore of it (_1b),
        # which differ wherever the original was heterozygous.
        keys = {canonical(n) for n in names if n}
        clash = sorted(k for k in keys if k in taken
                       and not (row["source_type"] == "matrix" == taken[k][1]))
        if clash:
            return taken[clash[0]][0]
        for key in keys:
            taken.setdefault(key, (row["strain"], row["source_type"]))
        rows.append(row)
        return None

    # 1. reads, the S288C family
    runs = {}
    with open(args.reads, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            runs.setdefault(row["strain"], []).append(row)
    for strain, items in runs.items():
        build = [r["run"] for r in items if r["role"] == "build"]
        heldout = [r["run"] for r in items if r["role"] == "heldout"]
        add({"strain": strain, "display_name": strain, "source_type": "reads",
             "source": ",".join(build), "clade": "S288C family",
             "ploidy": items[0].get("ploidy", ""),
             "zygosity": "", "origin": "Laboratory",
             "notes": "runs %s%s" % (", ".join("%s (%s)" % (r["run"], r["study"])
                                               for r in items if r["role"] == "build"),
                                     "; held out for validation: " + ", ".join(heldout)
                                     if heldout else "")},
            [strain])

    # 2. the 1011 matrix
    samples = [line.strip() for line in open(args.matrix_samples) if line.strip()]
    table = {}
    with open(args.matrix_table, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("Standardized name") in samples:
                table[row["Standardized name"]] = row
    used_ids = {r["strain"] for r in rows}
    for code in samples:
        info = table.get(code, {})
        isolate = info.get("Isolate name") or code
        strain = safe_id(isolate)
        if strain in used_ids:
            strain = "%s_%s" % (strain, code)
        duplicate = add({"strain": strain, "display_name": "%s (%s)" % (isolate, code),
                         "source_type": "matrix", "source": code,
                         "clade": info.get("Clades", ""), "ploidy": info.get("Ploidy", ""),
                         "zygosity": info.get("Zygosity", ""),
                         "origin": info.get("Ecological origins", ""),
                         "notes": "1011 Yeast Genomes, Peter et al. 2018, doi:10.1038/s41586-018-0030-5"},
                        [isolate, code])
        if duplicate:
            sys.stderr.write("matrix %s (%s) is %s, kept from its earlier source\n"
                             % (code, isolate, duplicate))
        used_ids.add(strain)

    # 3. assemblies: ScRAPdb (long reads) before SGD
    scrap = {}
    for name in sorted(os.listdir(args.scrapdb_dir)):
        match = re.match(r"^(.+?)\.(asm\d+)\.([^.]+)\.nuclear_genome", name)
        if not match:
            continue
        strain, asm, hap = match.groups()
        rank = (asm, {"HP0": 0, "collapsed": 1, "HP": 2, "HP1": 3, "HP2": 4}.get(hap, 9))
        if strain not in scrap or rank < scrap[strain][0]:
            scrap[strain] = (rank, os.path.join(args.scrapdb_dir, name))
    sgd = {}
    for root, _, files in os.walk(args.sgd_dir):
        if "/." in root.replace("\\", "/") or "archive" in root:
            continue
        for name in files:
            if not re.search(r"\.fsa(\.gz)?$", name):
                continue
            match = re.match(r"^(.+?)_([A-Za-z]+)_(\d{4})_", name)
            if not match:
                continue
            strain, centre, year = match.groups()
            # The SGD 2014-2015 re-sequencing is the most uniform set; within a
            # strain, the most recent assembly is preferred.
            rank = -int(year)
            if strain not in sgd or rank < sgd[strain][0]:
                sgd[strain] = (rank, os.path.join(root, name), centre, year)

    for strain, (_, path) in sorted(scrap.items()):
        if ("scrapdb", normalise(strain)) in excluded:
            continue
        duplicate = add({"strain": safe_id(strain), "display_name": strain,
                         "source_type": "assembly",
                         "source": written(path, args.scrapdb_dir), "clade": "", "ploidy": "",
                         "zygosity": "", "origin": "",
                         "notes": "ScRAPdb 2024-06-10, O'Donnell et al. 2023, doi:10.1038/s41588-023-01459-y"},
                        [strain])
        if duplicate:
            sys.stderr.write("ScRAPdb %s is %s, kept from its earlier source\n"
                             % (strain, duplicate))
    for strain, (_, path, centre, year) in sorted(sgd.items()):
        if ("sgd", normalise(strain)) in excluded:
            continue
        duplicate = add({"strain": safe_id(strain), "display_name": strain,
                         "source_type": "assembly",
                         "source": written(path, args.sgd_dir), "clade": "", "ploidy": "",
                         "zygosity": "", "origin": "",
                         "notes": "SGD strain sequences, %s %s" % (centre, year)},
                        [strain])
        if duplicate:
            sys.stderr.write("SGD %s is %s, kept from its earlier source\n"
                             % (strain, duplicate))

    columns = ["strain", "display_name", "source_type", "source", "clade", "ploidy",
               "zygosity", "origin", "notes"]
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    counts = {}
    for row in rows:
        counts[row["source_type"]] = counts.get(row["source_type"], 0) + 1
    print("%d strains: %s" % (len(rows), ", ".join("%s %d" % kv for kv in sorted(counts.items()))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
