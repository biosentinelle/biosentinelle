#!/usr/bin/env python3
"""
Build the native marker gene loci from the reference genome.

The selection markers that matter most in a laboratory yeast strain are not the
foreign cassettes, they are the strain's own URA3, HIS3, LEU2, TRP1, LYS2,
ADE2, MET15 and CAN1, and what has been done to them. Those genes are native,
so the tool element scan cannot say anything about them: a wild-type strain
covers them end to end.

This script writes one locus per gene: the ORF plus a flank on each side, in
the gene's coding orientation, so that build_marker_alleles.py and
genotype_markers.py can judge every allele by the protein it encodes.

Coordinates are not typed in. They are fetched from the Ensembl REST API by
systematic name, and the assembly Ensembl reports must be the one the reference
genome FASTA is (R64-1-1 for yeast). Each ORF is then checked against the
reference sequence: it must start with ATG, end with a stop, be a whole number
of codons and carry no internal stop. A gene that fails is refused rather than
written, because every downstream call is read in this frame.

Single-exon genes only. A spliced gene would need its exons stitched, and none
of the yeast markers has an intron; a gene with more than one exon is refused.
"""

import argparse
import csv
import json
import os
import re
import sys
import time
import urllib.request

from marker_effects import read_fasta, revcomp, translate

ENSEMBL = "https://rest.ensembl.org"


def fetch_json(url):
    request = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return json.loads(response.read().decode())
        except Exception as error:  # network hiccups are retried, then fatal
            if attempt == 3:
                raise
            sys.stderr.write("retrying %s: %s\n" % (url, error))
            time.sleep(2 * (attempt + 1))


def chromosome_names(reference_path):
    """Map chromosome labels (I, II, ..., Mito) to FASTA record names."""
    names, sequences = {}, {}
    for header, sequence in read_fasta(reference_path):
        record = header.split()[0]
        sequences[record] = sequence
        match = re.search(r"\[chromosome=([^\]]+)\]", header)
        if match:
            names[match.group(1)] = record
        elif "mitochondrion" in header:
            names["Mito"] = record
        names.setdefault(record, record)
    return names, sequences


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genes", required=True, help="marker_genes.tsv")
    parser.add_argument("--reference", required=True, help="reference genome FASTA")
    parser.add_argument("--species", default="saccharomyces_cerevisiae")
    parser.add_argument("--assembly", default="R64-1-1",
                        help="assembly name the reference FASTA is; Ensembl must agree")
    parser.add_argument("--flank", type=int, default=1000,
                        help="bases kept on each side of the ORF (default 1000)")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    names, sequences = chromosome_names(args.reference)
    release = fetch_json(ENSEMBL + "/info/data?content-type=application/json")["releases"]

    with open(args.genes, newline="") as handle:
        genes = list(csv.DictReader(handle, delimiter="\t"))

    rows, records, bed, failed = [], [], [], []
    for gene in genes:
        url = "%s/lookup/id/%s?content-type=application/json;expand=1" % (
            ENSEMBL, gene["systematic_id"])
        info = fetch_json(url)
        if info.get("assembly_name") != args.assembly:
            failed.append("%s: Ensembl assembly %s, expected %s"
                          % (gene["gene"], info.get("assembly_name"), args.assembly))
            continue
        exons = [len(t.get("Exon", [])) for t in info.get("Transcript", [])
                 if t.get("biotype") == "protein_coding"]
        if exons != [1]:
            failed.append("%s: %s protein coding exons, only single exon genes are handled"
                          % (gene["gene"], exons))
            continue
        chrom = names.get(info["seq_region_name"])
        if chrom is None:
            failed.append("%s: chromosome %s not in the reference"
                          % (gene["gene"], info["seq_region_name"]))
            continue
        chrom_seq = sequences[chrom]
        orf_start, orf_end = info["start"] - 1, info["end"]  # 0-based half open
        strand = "+" if info["strand"] == 1 else "-"
        locus_start = max(0, orf_start - args.flank)
        locus_end = min(len(chrom_seq), orf_end + args.flank)
        locus = chrom_seq[locus_start:locus_end]
        if strand == "+":
            up, down = orf_start - locus_start, locus_end - orf_end
        else:
            locus = revcomp(locus)
            up, down = locus_end - orf_end, orf_start - locus_start
        orf = locus[up:len(locus) - down]

        protein = translate(orf, to_stop=False)
        problems = []
        if not orf.startswith("ATG"):
            problems.append("no ATG start")
        if len(orf) % 3:
            problems.append("length %d is not a multiple of 3" % len(orf))
        if not protein.endswith("*"):
            problems.append("no stop codon at the end")
        if "*" in protein[:-1]:
            problems.append("internal stop at codon %d" % (protein.index("*") + 1))
        if problems:
            failed.append("%s: %s" % (gene["gene"], ", ".join(problems)))
            continue

        locus_id = "%s__%s" % (gene["gene"], gene["systematic_id"])
        rows.append({
            "locus_id": locus_id,
            "gene": gene["gene"],
            "systematic_id": gene["systematic_id"],
            "ensembl_name": info.get("display_name", ""),
            "chrom": chrom,
            "genome_start": locus_start,
            "genome_end": locus_end,
            "strand": strand,
            "orf_start": up,
            "orf_end": len(locus) - down,
            "orf_length": len(orf),
            "protein_length": len(protein) - 1,
            "loss_of_function_phenotype": gene["loss_of_function_phenotype"],
            "source": "%s/lookup/id/%s (Ensembl release %s, %s)"
                      % (ENSEMBL, gene["systematic_id"], ",".join(map(str, release)),
                         info.get("assembly_name")),
        })
        records.append((locus_id, locus))
        bed.append((chrom, locus_start, locus_end, locus_id))
        sys.stderr.write("%-6s %s %s:%d-%d (%s) ORF %d bp, %d aa\n"
                         % (gene["gene"], gene["systematic_id"], chrom, locus_start + 1,
                            locus_end, strand, len(orf), len(protein) - 1))

    if failed:
        sys.exit("refused:\n  " + "\n  ".join(failed))

    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "marker_loci.tsv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with open(os.path.join(args.outdir, "marker_loci.fasta"), "w", newline="\n") as handle:
        for locus_id, sequence in records:
            handle.write(">%s\n" % locus_id)
            for i in range(0, len(sequence), 80):
                handle.write(sequence[i:i + 80] + "\n")
    with open(os.path.join(args.outdir, "marker_loci.bed"), "w", newline="\n") as handle:
        for chrom, start, end, locus_id in sorted(bed):
            handle.write("%s\t%d\t%d\t%s\n" % (chrom, start, end, locus_id))
    sys.stderr.write("%d loci written to %s\n" % (len(rows), args.outdir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
