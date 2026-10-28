#!/usr/bin/env python3
"""
Offline end-to-end check of the strain identity chain. No aligner, no network,
no containers: a synthetic genome stands in for the reference, so the whole
build-and-score path can be exercised on any machine.

    python test/smoke_test_strain_id.py

Three controls; the last one matters most:

  positive   Reads simulated from a strain that is in the panel must be ranked
             first, with high concordance and no uncovered sites.

  relative   Reads from a strain that differs from a panel strain at four sites
             only must be ranked first: the panel must keep a strain's rare
             variants, which are all that separate clonal relatives.

  negative  Reads simulated from a strain that is NOT in the panel must come
             back unassigned. This is the case that a concordance threshold
             alone does not catch: an unknown strain whose variants miss the
             panel's sites looks exactly like the reference strain, and will be
             reported as such with full confidence unless the uncovered-site
             check is working. If this control regresses, the pipeline has
             started inventing answers.
"""

import json
import os
import random
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
GENOME_LEN = 120000
CONTIG = "chrTest"
READ_LEN = 150
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def run(script, args, work):
    command = [sys.executable, os.path.join(BIN, script)] + args
    result = subprocess.run(command, capture_output=True, text=True, cwd=work)
    if result.returncode != 0:
        print("FAILED: %s\n%s\n%s" % (script, result.stdout, result.stderr))
        sys.exit(1)
    return result.stdout + result.stderr


def write_fasta(path, sequence):
    with open(path, "w") as handle:
        handle.write(">%s\n" % CONTIG)
        for i in range(0, len(sequence), 70):
            handle.write(sequence[i:i + 70] + "\n")


def make_variants(reference, count, seed):
    rng = random.Random(seed)
    out = []
    for pos in sorted(rng.sample(range(100, len(reference) - 100), count)):
        ref_base = reference[pos - 1]
        out.append((pos, ref_base, rng.choice([b for b in "ACGT" if b != ref_base])))
    return out


def write_vcf(path, variants):
    with open(path, "w") as handle:
        handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for pos, ref_base, alt in variants:
            handle.write("%s\t%d\t.\t%s\t%s\t60\tPASS\t.\n" % (CONTIG, pos, ref_base, alt))


def simulate_reads(path, genome, depth, seed):
    rng = random.Random(seed)
    with open(path, "w") as handle:
        for i in range(len(genome) * depth // READ_LEN):
            start = rng.randrange(0, len(genome) - READ_LEN)
            read = genome[start:start + READ_LEN]
            read = "".join(
                rng.choice([b for b in "ACGT" if b != base]) if rng.random() < 0.003 else base
                for base in read)
            if rng.random() < 0.5:
                read = read.translate(COMPLEMENT)[::-1]
            handle.write("@read%d\n%s\n+\n%s\n" % (i, read, "I" * READ_LEN))


def apply_variants(reference, variants):
    genome = list(reference)
    for pos, _, alt in variants:
        genome[pos - 1] = alt
    return "".join(genome)


def main():
    work = tempfile.mkdtemp(prefix="biosentinelle_smoke_")
    path = lambda name: os.path.join(work, name)
    print("work directory: %s\n" % work)

    rng = random.Random(20261023)
    reference = "".join(rng.choice("ACGT") for _ in range(GENOME_LEN))
    write_fasta(path("ref.fa"), reference)

    # Two panel strains sharing some differences, as related lab strains do,
    # and a close relative of StrainX that differs from it at four sites only,
    # as BY4741 differs from S288C. StrainY also carries variants 5 bp from
    # some of StrainX's, so those sites need per-strain k-mers to be read.
    shared = make_variants(reference, 150, 3)
    strains = {}
    strains["StrainX"] = sorted(set(make_variants(reference, 900, 1)) | set(shared))
    near = []
    for pos, _, _ in strains["StrainX"][::20]:
        base = reference[pos + 4]
        near.append((pos + 5, base, "A" if base != "A" else "C"))
    strains["StrainY"] = sorted(set(make_variants(reference, 900, 2)) | set(shared) | set(near))
    taken = {pos for pos, _, _ in strains["StrainX"]}
    private = [v for v in make_variants(reference, 40, 5) if v[0] not in taken][:4]
    strains["StrainX2"] = sorted(set(strains["StrainX"]) | set(private))
    strains["RefStrain"] = []
    with open(path("strains.tsv"), "w") as handle:
        handle.write("strain\tdisplay_name\tsource_type\tsource\tclade\tploidy\t"
                     "zygosity\torigin\tnotes\n")
        for name in strains:
            handle.write("%s\t%s\tassembly\t%s.vcf\t\t1\t\t\t\n" % (name, name, name))
    with open(path("all.bed"), "w") as handle:
        handle.write("%s\t0\t%d\n" % (CONTIG, GENOME_LEN))
    singles = []
    for name, variants in strains.items():
        write_vcf(path("%s.vcf" % name), variants)
        singles += ["--single", "%s=%s=%s" % (name, path("%s.vcf" % name), path("all.bed"))]

    print(run("build_population_panel.py", [
        "--strains", path("strains.tsv"), "--reference", path("ref.fa"), "-k", "31",
        "--rare-max-carriers", "1", "--per-strain-sites", "5000",
        "--common-sites", "5000"] + singles + [
        "--out-sites", path("panel_sites.tsv"), "--out-kmers", path("panel_kmers.tsv"),
        "--out-strains", path("panel_strains.tsv")], work))

    failures = []

    def score(reads, label, out_json):
        run("scan_kmers.py", ["--kmers", path("panel_kmers.tsv"), "--reads", path(reads),
                              "--out", path(out_json + ".calls"), "--max-reads", "0",
                              "--progress-every", "0"], work)
        run("score_strains.py", ["--sites", path("panel_sites.tsv"),
                                 "--strains", path("panel_strains.tsv"),
                                 "--calls", path(out_json + ".calls"),
                                 "--min-discriminating-sites", "3",
                                 "--sample", label, "--out", path(out_json)], work)
        return json.load(open(path(out_json)))

    # --- positive control ---------------------------------------------------
    simulate_reads(path("reads_x.fq"), apply_variants(reference, strains["StrainX"]), 25, 99)
    result = score("reads_x.fq", "POSITIVE", "x.json")
    best = result["strains"][0]
    print("positive control: %s p=%.3f concordance=%.3f uncovered=%.1f%% unassigned=%s" % (
        best["strain_id"], best["probability"], best["concordance"],
        100 * result["coverage"]["uncovered_fraction"], result["unassigned"]))
    if best["strain_id"] != "StrainX":
        failures.append("positive control ranked %s first, expected StrainX" % best["strain_id"])
    if result["unassigned"]:
        failures.append("positive control was reported unassigned")
    if best["concordance"] < 0.95:
        failures.append("positive control concordance %.3f below 0.95" % best["concordance"])

    # --- close relative -----------------------------------------------------
    # Four private sites are all that separate StrainX2 from StrainX. They are
    # rare by construction, so the panel must keep them, and the reads must
    # put StrainX2 first.
    simulate_reads(path("reads_x2.fq"), apply_variants(reference, strains["StrainX2"]), 25, 31)
    result = score("reads_x2.fq", "RELATIVE", "x2.json")
    best = result["strains"][0]
    print("close relative:   %s p=%.3f concordance=%.3f, %d sites from %s" % (
        best["strain_id"], best["probability"], best["concordance"],
        best.get("discriminating_sites_vs_next", -1), best.get("runner_up")))
    if best["strain_id"] != "StrainX2":
        failures.append("close relative ranked %s first, expected StrainX2" % best["strain_id"])

    # --- negative control ---------------------------------------------------
    off_panel = make_variants(reference, 1200, 4242)
    simulate_reads(path("reads_z.fq"), apply_variants(reference, off_panel), 25, 7)
    result = score("reads_z.fq", "NEGATIVE", "z.json")
    best = result["strains"][0]
    print("negative control: %s p=%.3f concordance=%.3f uncovered=%.1f%% unassigned=%s" % (
        best["strain_id"], best["probability"], best["concordance"],
        100 * result["coverage"]["uncovered_fraction"], result["unassigned"]))
    if not result["unassigned"]:
        failures.append("negative control was assigned to %s with concordance %.3f; "
                        "a strain outside the panel must not be matched to a panel "
                        "member" % (best["strain_id"], best["concordance"]))

    # --- report rendering ---------------------------------------------------
    run("render_report.py", ["--json", path("x.json"), "--out", path("report.html")], work)
    size = os.path.getsize(path("report.html"))
    print("report rendered: %d bytes" % size)
    if size < 5000:
        failures.append("rendered report is suspiciously small (%d bytes)" % size)

    print()
    if failures:
        for failure in failures:
            print("FAIL: %s" % failure)
        return 1
    print("PASS: all three controls behaved correctly")
    return 0


if __name__ == "__main__":
    sys.exit(main())
