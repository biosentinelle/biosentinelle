#!/usr/bin/env python3
"""
Writes the sequence files of galaxy/test-data used by the tests of the Galaxy tool (galaxy/lab_strain_watch.xml),
from the reference files of the repository, so that they can be rebuilt identically (fixed random seed, gzip without date):

    S288C_chrIII.fasta.gz     chromosome III of the S288C reference genome (references/genome), an assembly
                              carrying the MAT, HML and HMR loci and the LEU2 marker gene (test 1, seq_kind = "assembly")
    chrIII_MAT_region.fasta   chrIII:192001-207000 of S288C, around the MAT locus (test 2, own genome of the variant counts)
    reads_R1.fastq.gz         simulated paired-end reads (2 x 150 bp, fragments of 300 bp on average, 30x, 0.2% substitutions)
    reads_R2.fastq.gz         of chrIII_MAT_region and of the pieces of the pFA6a-3HA-His3MX6 cassette of references/features
                              joined end to end (test 2), so that the features and mating type sequences of the pipeline are covered

The other files of galaxy/test-data (parameter file, expected lines) are written by hand.

Usage, from any folder:
    python3 galaxy/local/make_test_data.py
"""

import glob
import gzip
import os
import random

REPO_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT_DIR = os.path.join(REPO_DIR, "galaxy", "test-data")
GENOME = os.path.join(REPO_DIR, "references", "genome", "S288C_reference_sequence_R64-1-1_20110203.fsa")
FEATURES = sorted(glob.glob(os.path.join(REPO_DIR, "references", "features", "*.fasta")))
CHROMOSOME = "NC_001135"  # chromosome III
REGION = (192000, 207000)  # 0-based, end excluded: chrIII:192001-207000
READ_LENGTH = 150
FRAGMENT_MEAN, FRAGMENT_SD = 300, 30
DEPTH = 30
ERROR_RATE = 0.002
SEED = 1
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def read_fasta(path):
    """Returns [(header, sequence)] of a FASTA file."""
    records, header, chunks = [], None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks).upper()))
                header, chunks = line[1:], []
            elif line:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks).upper()))
    return records


def fasta_text(header, sequence, width=70):
    return ">" + header + "\n" + "".join(sequence[i:i + width] + "\n" for i in range(0, len(sequence), width))


def write_gzip(path, text):
    """gzip without file name nor date in its header: same bytes at each run."""
    with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
        handle.write(text.encode("ascii"))


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def mutate(sequence, rng):
    bases = list(sequence)
    for i, base in enumerate(bases):
        if rng.random() < ERROR_RATE:
            bases[i] = rng.choice([b for b in "ACGT" if b != base])
    return "".join(bases)


def simulate_pairs(sequences, rng):
    """Paired-end reads of the sequences, at DEPTH: returns the R1 and R2 FASTQ texts."""
    total = sum(len(s) for s in sequences)
    weights = [len(s) / total for s in sequences]
    pairs = total * DEPTH // (2 * READ_LENGTH)
    r1, r2 = [], []
    for n in range(1, pairs + 1):
        sequence = rng.choices(sequences, weights)[0]
        length = int(rng.gauss(FRAGMENT_MEAN, FRAGMENT_SD))
        length = max(READ_LENGTH + 20, min(length, FRAGMENT_MEAN + 3 * FRAGMENT_SD, len(sequence)))
        start = rng.randint(0, len(sequence) - length)
        fragment = sequence[start:start + length]
        if rng.random() < 0.5:
            fragment = reverse_complement(fragment)
        read1 = mutate(fragment[:READ_LENGTH], rng)
        read2 = mutate(reverse_complement(fragment)[:READ_LENGTH], rng)
        quality = "I" * READ_LENGTH
        r1.append(f"@sim{n}/1\n{read1}\n+\n{quality}\n")
        r2.append(f"@sim{n}/2\n{read2}\n+\n{quality}\n")
    return "".join(r1), "".join(r2)


def main():
    rng = random.Random(SEED)
    genome = {header.split("|")[1] if header.startswith("ref|") else header.split()[0]: seq for header, seq in read_fasta(GENOME)}
    chromosome = genome[CHROMOSOME]

    path = os.path.join(OUT_DIR, "S288C_chrIII.fasta.gz")
    write_gzip(path, fasta_text(f"{CHROMOSOME} Saccharomyces cerevisiae S288C chromosome III", chromosome))
    print(f"{path}: {len(chromosome)} bp")

    region = chromosome[REGION[0]:REGION[1]]
    path = os.path.join(OUT_DIR, "chrIII_MAT_region.fasta")
    with open(path, "w") as handle:
        handle.write(fasta_text(f"chrIII_{REGION[0] + 1}_{REGION[1]} S288C {CHROMOSOME}:{REGION[0] + 1}-{REGION[1]}", region))
    print(f"{path}: {len(region)} bp")

    cassette = "".join(seq for path in FEATURES for _, seq in read_fasta(path))
    read1, read2 = simulate_pairs([region, cassette], rng)
    for name, text in (("reads_R1.fastq.gz", read1), ("reads_R2.fastq.gz", read2)):
        path = os.path.join(OUT_DIR, name)
        write_gzip(path, text)
        print(f"{path}: {text.count(chr(10)) // 4} reads")


if __name__ == "__main__":
    main()
