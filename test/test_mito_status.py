#!/usr/bin/env python3
"""
Offline check of the mitochondrial genome (rho) call, hence growth on glycerol.

No aligner: the depth over the mtDNA of the reference and the read counts per
chromosome are written directly, as samtools depth and samtools idxstats would
give them, for the situations the call must tell apart:

    rho+        every respiratory gene covered, at 20 times the nuclear depth
    rho+ intron the optional introns of COX1 and COB missing, exons covered:
                strains differ in introns, this must stay rho+
    rho+ mixed  one region amplified 20 fold over the rest: rho- cells in a
                rho+ culture, rho+ with a note
    rho-        only the COX2 region left, amplified
    rho0        no read on the mtDNA
    uncertain   VAR1 covered at half its exon length only
    no data     nuclear depth under the threshold
    RNA         not assessed

    python test/test_mito_status.py
"""

import csv
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
GENES = os.path.join(os.path.dirname(ROOT), "biosentinelle_db", "db", "yeast", "mito",
                     "mito_genes.tsv")
MITO_LEN = 85779
NUCLEAR = [("ref|NC_001133|", 230218), ("ref|NC_001134|", 813184)]
READ_LEN = 35


def load_genes():
    with open(GENES, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def exons(row):
    return [tuple(map(int, x.split("-"))) for x in row["exons"].split(";")]


def run_case(work, label, depth_of, nuclear_depth=30.0, extra=()):
    genes = load_genes()
    chrom = genes[0]["chrom"]
    depth_path = os.path.join(work, label + ".depth")
    with open(depth_path, "w") as handle:
        for p in range(MITO_LEN):
            handle.write("%s\t%d\t%d\n" % (chrom, p + 1, depth_of(p)))
    mito_reads = int(sum(depth_of(p) for p in range(0, MITO_LEN, 100)) * 100 / READ_LEN)
    idx_path = os.path.join(work, label + ".idxstats")
    with open(idx_path, "w") as handle:
        for name, length in NUCLEAR:
            handle.write("%s\t%d\t%d\t0\n" % (name, length, int(length * nuclear_depth / READ_LEN)))
        handle.write("%s\t%d\t%d\t0\n*\t0\t0\t0\n" % (chrom, MITO_LEN, mito_reads))
    out = os.path.join(work, label + ".json")
    subprocess.run([sys.executable, os.path.join(BIN, "mito_status.py"), "--genes", GENES,
                    "--depth", depth_path, "--idxstats", idx_path,
                    "--read-length", str(READ_LEN), "--out", out] + list(extra),
                   check=True, stderr=subprocess.DEVNULL)
    with open(out) as handle:
        return json.load(handle)


def main():
    genes = {g["gene"]: g for g in load_genes()}
    in_exon = set()
    for g in genes.values():
        for a, b in exons(g):
            in_exon.update(range(a, b))
    introns = set()
    for name in ("COX1", "COB"):
        ex = exons(genes[name])
        for (a, b), (c, d) in zip(ex, ex[1:]):
            introns.update(range(b, c))
    cox2 = genes["COX2"]
    cox2_region = range(int(cox2["start"]) - 2000, int(cox2["end"]) + 2000)
    var1 = exons(genes["VAR1"])[0]
    var1_half = range(var1[0], (var1[0] + var1[1]) // 2)
    amplified = range(int(genes["ATP6"]["start"]) - 500, int(genes["ATP6"]["end"]) + 500)

    cases = [
        ("rho+", lambda p: 600, (), "rho+", "growth expected"),
        ("rho+ intron", lambda p: 0 if p in introns else 600, (), "rho+", "growth expected"),
        ("rho+ mixed", lambda p: 12000 if p in amplified else 600, (), "rho+",
         "growth expected"),
        ("rho-", lambda p: 15000 if p in cox2_region else 0, (), "rho-", "no growth"),
        ("rho0", lambda p: 0, (), "rho0", "no growth"),
        ("uncertain", lambda p: 0 if p in var1_half else 600, (), "uncertain", "unknown"),
        ("no data", lambda p: 600, ("--min-depth", "4"), "no_data", "unknown"),
        ("RNA", lambda p: 600, ("--seq-molecule", "RNA"), "not_assessed", "unknown"),
    ]
    failed = 0
    with tempfile.TemporaryDirectory() as work:
        for label, depth_of, extra, status, glycerol in cases:
            nuclear = 2.0 if label == "no data" else 30.0
            result = run_case(work, label.replace(" ", "_").replace("+", "p"), depth_of,
                              nuclear, extra)
            ok = result["status"] == status and result["glycerol"] == glycerol
            if label == "rho+ mixed":
                ok = ok and bool(result.get("notes"))
            if label == "rho-":
                # Only the COX2 region is left: the reason names what is kept.
                ok = ok and "kept:" in result["reason"] and "COX2" in result["reason"]
            failed += not ok
            print("%-4s %-12s -> %-12s glycerol: %-15s %s" % (
                "ok" if ok else "FAIL", label, result["status"], result["glycerol"],
                result["reason"][:90]))
    print("%d of %d cases failed" % (failed, len(cases)) if failed else "all cases pass")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
