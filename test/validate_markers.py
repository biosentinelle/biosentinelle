#!/usr/bin/env python3
"""
Check native marker genotyping against strains of known genotype.

Reads are simulated from the strain assemblies around every marker locus, put
through the same alignment and genotype_markers.py as a real sample, and the
call is compared with the documented genotype of the strain. This is the test
that matters: it exercises point mutations (W303, YPH499 ochre and amber
alleles), small frameshifts (his3-11,15, leu2-3,112), internal and complete
deletions (BY4742) and a transposon insertion (ura3-52), all from reads.

It needs the ScRAPdb assemblies and a bowtie2 index of the reference, so it
is not part of the offline smoke test:

    python3 test/validate_markers.py \\
        --assembly-dir dataset/genome_references/ScRAPdb_20240610_nuclear_genome \\
        --index <bowtie2 index prefix> --reference <reference FASTA> \\
        --run 'apptainer exec -B /mnt/c <bowtie2 image>' --workdir <dir>

A real strain is never this clean: the simulated reads carry substitution
errors, but no library bias, no contamination and no copy number variation.
"""

import argparse
import csv
import json
import os
import random
import shlex
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))
from marker_effects import read_fasta, revcomp  # noqa: E402

DB = os.path.join(ROOT, "db", "yeast", "markers")


def load_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def simulate(regions, read_len, depth, error, rng, handle):
    n = 0
    for seq in regions:
        for _ in range(int(len(seq) * depth / read_len)):
            start = rng.randrange(0, max(1, len(seq) - read_len))
            read = list(seq[start:start + read_len])
            for i in range(len(read)):
                if rng.random() < error:
                    read[i] = rng.choice([b for b in "ACGT" if b != read[i]])
            read = "".join(read)
            if rng.random() < 0.5:
                read = revcomp(read)
            handle.write("@r%d\n%s\n+\n%s\n" % (n, read, "I" * len(read)))
            n += 1
    return n


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--assembly-dir", required=True)
    parser.add_argument("--index", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--run", default="", help="prefix for bowtie2 and samtools")
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--strains", nargs="+",
                        default=["BY4742", "W303", "YPH499", "CENPK21D", "S288C"])
    parser.add_argument("--read-len", type=int, default=100)
    parser.add_argument("--depth", type=float, default=25)
    parser.add_argument("--error", type=float, default=0.003)
    parser.add_argument("--flank", type=int, default=4000)
    args = parser.parse_args()

    catalogue = load_tsv(os.path.join(DB, "marker_alleles.tsv"))
    documented = dict((r["strain"], r) for r in load_tsv(os.path.join(DB, "documented_genotypes.tsv")))
    os.makedirs(args.workdir, exist_ok=True)
    prefix = shlex.split(args.run)
    bed = os.path.join(DB, "marker_loci.bed")
    rng = random.Random(1)
    failures = 0

    for strain in args.strains:
        assembly = os.path.join(args.assembly_dir, "%s.asm01.HP0.nuclear_genome.tidy.fa" % strain)
        contigs = dict((n.split()[0], s) for n, s in read_fasta(assembly))
        # The contig and span the catalogue builder chose for each locus, so the
        # reads come from the same place the catalogue allele was read from.
        regions = []
        for row in catalogue:
            if row["assembly"] == strain and row["contig"]:
                seq = contigs[row["contig"]]
                start, end = int(row["contig_start"]), int(row["contig_end"])
                regions.append(seq[max(0, start - args.flank):end + args.flank])

        work = os.path.join(args.workdir, strain)
        os.makedirs(work, exist_ok=True)
        fastq = os.path.join(work, "reads.fastq")
        with open(fastq, "w") as handle:
            simulate(regions, args.read_len, args.depth, args.error, rng, handle)

        bam = os.path.join(work, "markers.bam")
        script = (
            "{p} bowtie2 --sensitive-local -p 4 -x {idx} -U {fq} 2> {w}/bt2.log"
            " | {p} samtools view -b -L {bed} - | {p} samtools sort -o {bam} - &&"
            " {p} samtools index {bam} &&"
            " {p} samtools mpileup -B -q 20 -Q 13 -d 10000 -l {bed} -f {ref} {bam} > {w}/markers.pileup &&"
            " {p} samtools depth -a -b {bed} -Q 0 {bam} > {w}/markers.depth &&"
            " {p} samtools view -q 1 {bam} > {w}/markers.sam"
        ).format(p=" ".join(prefix), idx=args.index, fq=fastq, w=work, bed=bed, bam=bam,
                 ref=args.reference)
        subprocess.run(["bash", "-c", script], check=True,
                       stderr=subprocess.DEVNULL)
        subprocess.run([sys.executable, os.path.join(ROOT, "bin", "genotype_markers.py"),
                        "--loci", os.path.join(DB, "marker_loci.tsv"),
                        "--loci-fasta", os.path.join(DB, "marker_loci.fasta"),
                        "--alleles", os.path.join(DB, "marker_alleles.tsv"),
                        "--pileup", os.path.join(work, "markers.pileup"),
                        "--depth", os.path.join(work, "markers.depth"),
                        "--sam", os.path.join(work, "markers.sam"),
                        "--repeats", os.path.join(DB, "marker_repeats.tsv"),
                        "--tools", os.path.join(ROOT, "db", "yeast", "tools", "tools.fasta"),
                        "--tools-manifest", os.path.join(ROOT, "db", "yeast", "tools", "tools.tsv"),
                        "--out", os.path.join(work, "markers.json")],
                       check=True, stderr=subprocess.DEVNULL)
        with open(os.path.join(work, "markers.json")) as handle:
            result = json.load(handle)

        tokens = documented.get(strain, {}).get("documented_genotype", "").split()
        print("\n%s  documented: %s" % (strain, " ".join(tokens)))
        print("%s  called:     %s" % (" " * len(strain), result["genotype"]))
        for record in result["markers"]:
            gene = record["gene"]
            # Documented genotype when there is one, otherwise what the strain's
            # own assembly carries according to the catalogue.
            if tokens:
                hits = [t for t in tokens if t.lower().startswith(gene.lower())
                        and not t.startswith(gene.upper())]
                expected = "inactive" if hits else "active"
            else:
                own = [r for r in catalogue if r["assembly"] == strain and r["gene"] == gene]
                expected = own[0]["activity"] if own else "?"
                if expected == "uncertain":
                    expected = "active"  # a natural form of an isolate
            ok = record["status"] == expected
            failures += 0 if ok else 1
            allele = record["alleles"][0] if record["alleles"] else {}
            print("  %-4s %-6s expected %-8s called %-9s %-14s %s | %s" % (
                "ok" if ok else "FAIL", gene, expected, record["status"],
                allele.get("effect", ""), allele.get("protein_change", ""),
                allele.get("reason", "")))
    print("\n%d mismatches" % failures)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
