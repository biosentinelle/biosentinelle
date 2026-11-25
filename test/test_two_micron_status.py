#!/usr/bin/env python3
"""
Offline check of the 2-micron plasmid status (cir+, cir0) and copy number.

No aligner: the depth over the plasmid and the read counts per chromosome are
written directly, as samtools depth -aa and samtools idxstats would give them,
for the situations the call must tell apart (nuclear depth 30x):

    cir+           50 copies per nuclear genome over the whole plasmid
    cir+ mito      the same, with a mitochondrial chromosome at 20 times the
                   nuclear depth, left out of the nuclear depth
    cir0           no read on the plasmid
    stray reads    1x over the plasmid, under the floor: cir0
    low copy       2 copies per nuclear genome: cir+ with a note
    YEp fragment   only the region of the origin and STB and the repeats,
                   RAF1 not covered: partial, missing features named
    variant        half of REP2 not covered: partial, every feature there
    no data        nuclear depth under the threshold
    RNA reads      transcripts over the genes: present, no copy number
    RNA none       no transcript: not assessed
    assembly       a contig covers the plasmid: present, no copy number
    assembly none  no 2-micron contig: not assessed
    not aligned    the 2-micron sequence was not in the alignment

    python test/test_two_micron_status.py [db folder, default ../yeast_db/db]
"""

import csv
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
DB = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(ROOT), "yeast_db", "db")
FEATURES = os.path.join(DB, "yeast", "two_micron", "two_micron_features.tsv")
LENGTH = 6318
NUCLEAR = [("ref|NC_001133|", 230218), ("ref|NC_001134|", 813184)]
MITO = ("ref|NC_001224|", 85779)
READ_LEN = 35


def load_features():
    with open(FEATURES, newline="") as handle:
        return {row["feature"]: row for row in csv.DictReader(handle, delimiter="\t")}


def span(row):
    out = set()
    for item in row["segments"].split(";"):
        a, b = map(int, item.split("-"))
        out.update(range(a, b))
    return out


def run_case(work, label, depth_of, nuclear_depth=30.0, extra=(), mito=False, aligned=True):
    chrom = next(iter(load_features().values()))["chrom"]
    depth_path = os.path.join(work, label + ".depth")
    with open(depth_path, "w") as handle:
        if aligned:
            for p in range(LENGTH):
                handle.write("%s\t%d\t%d\n" % (chrom, p + 1, depth_of(p)))
    idx_path = os.path.join(work, label + ".idxstats")
    with open(idx_path, "w") as handle:
        for name, length in NUCLEAR:
            handle.write("%s\t%d\t%d\t0\n" % (name, length, int(length * nuclear_depth / READ_LEN)))
        if mito:
            handle.write("%s\t%d\t%d\t0\n" % (MITO[0], MITO[1],
                                              int(MITO[1] * 20 * nuclear_depth / READ_LEN)))
        handle.write("*\t0\t0\t0\n")
    out = os.path.join(work, label + ".json")
    command = [sys.executable, os.path.join(BIN, "two_micron_status.py"), "--features", FEATURES,
               "--depth", depth_path, "--idxstats", idx_path, "--read-length", str(READ_LEN),
               "--out", out] + list(extra)
    if mito:
        command += ["--mito-chrom", MITO[0]]
    subprocess.run(command, check=True, stderr=subprocess.DEVNULL)
    with open(out) as handle:
        return json.load(handle)


def main():
    features = load_features()
    genes = set()
    for name in ("FLP1", "REP1", "REP2", "RAF1"):
        genes |= span(features[name])
    # A 2-micron vector carries the origin and STB, between RAF1 and IR1, and
    # its reads of the repeat land on both copies of the repeat.
    yep = set(range(int(features["RAF1"]["end"]) + 100, int(features["IR1"]["end"]))) \
        | span(features["IR2"])
    rep2 = sorted(span(features["REP2"]))
    rep2_half = set(rep2[:len(rep2) // 2])

    cases = [
        # label, depth, options, status, cir, check
        ("cir+", lambda p: 1500, {}, "present", "cir+",
         lambda r: abs(r["copies_per_nuclear_genome"] - 50) < 1 and not r["notes"]),
        ("cir+ mito", lambda p: 1500, {"mito": True}, "present", "cir+",
         lambda r: abs(r["copies_per_nuclear_genome"] - 50) < 1),
        ("cir0", lambda p: 0, {}, "absent", "cir0",
         lambda r: "copies_per_nuclear_genome" not in r),
        ("stray reads", lambda p: 1, {}, "absent", "cir0", None),
        ("low copy", lambda p: 60, {}, "present", "cir+",
         lambda r: abs(r["copies_per_nuclear_genome"] - 2) < 0.1 and r["notes"]),
        ("YEp fragment", lambda p: 600 if p in yep else 0, {}, "partial", None,
         lambda r: "RAF1 (0%)" in r["reason"].split("missing:")[1]
         and "IR1 (100%)" in r["reason"].split("missing:")[0]
         and "2-micron vector" in r["reason"]),
        ("variant", lambda p: 0 if p in rep2_half else 1500, {}, "partial", None,
         lambda r: "REP2" in r["reason"] and "missing" not in r["reason"]),
        ("no data", lambda p: 1500, {"nuclear": 2.0}, "no_data", None, None),
        ("RNA reads", lambda p: 200 if p in genes else 0, {"extra": ("--seq-molecule", "RNA")},
         "present", None, lambda r: "copies_per_nuclear_genome" not in r),
        ("RNA none", lambda p: 0, {"extra": ("--seq-molecule", "RNA")}, "not_assessed", None,
         None),
        ("assembly", lambda p: 2, {"extra": ("--seq-kind", "assembly"), "nuclear": 2.0},
         "present", "cir+", lambda r: "copies_per_nuclear_genome" not in r),
        ("assembly none", lambda p: 0, {"extra": ("--seq-kind", "assembly"), "nuclear": 2.0},
         "not_assessed", None, None),
        ("not aligned", lambda p: 0, {"aligned": False}, "not_assessed", None, None),
    ]
    failed = 0
    with tempfile.TemporaryDirectory() as work:
        for label, depth_of, options, status, cir, check in cases:
            result = run_case(work, label.replace(" ", "_").replace("+", "p"), depth_of,
                              options.get("nuclear", 30.0), options.get("extra", ()),
                              options.get("mito", False), options.get("aligned", True))
            ok = result["status"] == status and result.get("cir") == cir
            ok = ok and (check is None or check(result))
            failed += not ok
            copies = result.get("copies_per_nuclear_genome")
            print("%-4s %-13s -> %-12s %-5s %-6s %s" % (
                "ok" if ok else "FAIL", label, result["status"], result.get("cir") or "",
                "" if copies is None else "%.1f" % copies, result["reason"][:80]))
    print("%d of %d cases failed" % (failed, len(cases)) if failed else "all cases pass")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
