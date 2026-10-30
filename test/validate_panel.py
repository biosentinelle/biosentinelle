#!/usr/bin/env python3
"""
Check the strain panel against genomes it was not built from.

For each case, reads are simulated from a genome assembly and scored against
the panel; the expected strain must come first. The point is to use a genome
that did not go into the panel for that strain: ScRAPdb holds long-read
assemblies of about 130 isolates that the catalogue takes from the 1011
matrix instead (short reads, another laboratory, another caller), so reads
from the assembly against the panel's matrix strain is an independent check
of both the panel and the matrix conversion.

    python3 test/validate_panel.py --cases cases.tsv --panel-dir ../biosentinelle_db/db/yeast/panels \\
        --workdir <scratch> --threads 8

cases.tsv: expected_strain, assembly path. Reads are 100 bp at 20x unless set
otherwise, with 0.3 % substitution errors, from both strands.
"""

import argparse
import csv
import gzip
import json
import os
import random
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def read_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    sequences, chunks = [], []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if chunks:
                    sequences.append("".join(chunks))
                chunks = []
            else:
                chunks.append(line.strip().upper())
    if chunks:
        sequences.append("".join(chunks))
    return sequences


def simulate(assembly, out, depth, length, seed):
    rng = random.Random(seed)
    contigs = [s for s in read_fasta(assembly) if len(s) > length]
    total = sum(len(s) for s in contigs)
    with open(out, "w") as handle:
        n = 0
        for contig in contigs:
            for _ in range(len(contig) * depth // length):
                start = rng.randrange(0, len(contig) - length)
                read = contig[start:start + length]
                if rng.random() < 0.5:
                    read = read.translate(COMPLEMENT)[::-1]
                read = "".join(rng.choice("ACGT") if rng.random() < 0.003 else b for b in read)
                handle.write(">r%d\n%s\n" % (n, read))
                n += 1
    return total


def run_case(case, args):
    name = case["expected_strain"]
    work = os.path.join(args.workdir, name)
    os.makedirs(work, exist_ok=True)
    reads = os.path.join(work, "reads.fa")
    if not os.path.exists(reads):
        simulate(case["assembly"], reads, args.depth, args.length, hash(name) & 0xFFFF)
    calls = os.path.join(work, "calls.tsv")
    if args.rescan or not os.path.exists(calls):
        subprocess.run([sys.executable, os.path.join(BIN, "scan_kmers.py"),
                        "--kmers", os.path.join(args.panel_dir, "panel_kmers.tsv.gz"),
                        "--reads", reads, "--out", calls, "--progress-every", "0"],
                       check=True, capture_output=True)
    out = os.path.join(work, "strain_id.json")
    subprocess.run([sys.executable, os.path.join(BIN, "score_strains.py"),
                    "--sites", os.path.join(args.panel_dir, "panel_sites.tsv.gz"),
                    "--strains", os.path.join(args.panel_dir, "panel_strains.tsv"),
                    "--calls", calls, "--sample", name, "--detail-top", "0",
                    "--out-ranking", os.path.join(work, "ranking.tsv"), "--out", out],
                   check=True, capture_output=True)
    report = json.load(open(out))
    ranking = list(csv.DictReader(open(os.path.join(work, "ranking.tsv")), delimiter="\t"))
    rank = next((int(r["rank"]) for r in ranking if r["strain_id"] == name), None)
    best = report["strains"][0]
    tied = [row["strain"] for row in report["resolution"]["tied"]]
    return {"tied": tied, "in_tied": name in tied,"expected": name, "rank_of_expected": rank, "best": best["strain_id"],
            "concordance": best["concordance"],
            "separating_vs_next": best.get("discriminating_sites_vs_next"),
            "runner_up": best.get("runner_up"),
            "uncovered": report["coverage"]["uncovered_fraction"],
            "unassigned": report["unassigned"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cases", required=True)
    parser.add_argument("--panel-dir", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--depth", type=int, default=20)
    parser.add_argument("--length", type=int, default=100)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--rescan", action="store_true",
                        help="scan the reads again even when a calls file exists "
                             "(needed after the panel changes)")
    args = parser.parse_args()

    cases = list(csv.DictReader(open(args.cases), delimiter="\t"))
    with ThreadPoolExecutor(args.threads) as pool:
        results = list(pool.map(lambda c: run_case(c, args), cases))
    print("expected\trank\tbest\tconcordance\tseparating_vs_next\trunner_up\tuncovered"
          "\tunassigned\ttied_with_best")
    correct = within = 0
    for r in results:
        correct += r["rank_of_expected"] == 1
        within += r["rank_of_expected"] == 1 or r["in_tied"]
        print("%s\t%s\t%s\t%.4f\t%s\t%s\t%.3f\t%s\t%s" % (
            r["expected"], r["rank_of_expected"], r["best"], r["concordance"],
            r["separating_vs_next"], r["runner_up"], r["uncovered"], r["unassigned"],
            ",".join(r["tied"])))
    print("\n%d of %d ranked first; %d of %d first or among the strains reported as not "
          "separable from the first" % (correct, len(results), within, len(results)))
    return 0 if correct == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
