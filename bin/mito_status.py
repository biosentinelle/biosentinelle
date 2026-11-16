#!/usr/bin/env python3
"""
Mitochondrial genome status of a sample, and what it says about growth on
glycerol.

Growth on glycerol, or any non-fermentable carbon source, needs respiration,
hence a mitochondrial genome (mtDNA) that still encodes the respiratory chain.
Laboratory strains lose it easily, and a petite strain grows on glucose like
any other, so the loss goes unnoticed until the strain is put on glycerol:

    rho+   every respiratory gene covered          glycerol: growth expected
    rho-   only part of the mtDNA left, often      glycerol: no growth
           amplified many times
    rho0   no mtDNA at all                         glycerol: no growth
    uncertain  every gene there, some only in part glycerol: unknown

The reads are aligned to the whole reference genome, mtDNA included, and the
depth is read over the exons of the genes needed for respiration
(biosentinelle_db bin/build_mito_genes.py). Exons only: the group I and II
introns of COX1, COB and the 21S rRNA are optional and strains differ in which
they carry, so a missing intron is not a rho- deletion.

A gene is present when most of its exon bases are covered, the floor being a
fraction of the nuclear depth: mtDNA is in tens of copies per cell, so a rho+
strain covers its genes far above the nuclear depth.

Not seen from depth: a point mutation in a mitochondrial gene (mit-), and
nuclear mutations that stop respiration (pet). A culture where some cells are
rho- shows as rho+, with the genes the rho- cells kept over-represented; that
is reported as a note.
"""

import argparse
import csv
import json
import statistics
import sys


def load_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_idxstats(path):
    """contig -> (length, mapped reads)."""
    out = {}
    with open(path) as handle:
        for line in handle:
            fields = line.split("\t")
            if len(fields) >= 3 and fields[0] != "*":
                out[fields[0]] = (int(fields[1]), int(fields[2]))
    return out


def parse_depth(path, chrom):
    out = {}
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 3 and fields[0] == chrom:
                out[int(fields[1]) - 1] = int(fields[2])
    return out


def exon_positions(row):
    out = []
    for item in row["exons"].split(";"):
        a, b = item.split("-")
        out.extend(range(int(a), int(b)))
    return out


def assess(genes, depth, nuclear_depth, args):
    """Per-gene coverage and the rho call. Returns the result dict."""
    if args.seq_kind == "assembly":
        # One copy of everything: a base is covered or it is not.
        floor = 1
    else:
        floor = max(args.min_depth, args.floor_fraction * nuclear_depth)
    per_gene = []
    for row in genes:
        values = [depth.get(p, 0) for p in exon_positions(row)]
        covered = sum(1 for v in values if v >= floor)
        per_gene.append({
            "gene": row["gene"], "role": row["role"],
            "exon_bp": len(values),
            "covered_fraction": round(covered / float(len(values)), 3) if values else 0.0,
            "median_depth": statistics.median(values) if values else 0,
        })
    total_bp = sum(g["exon_bp"] for g in per_gene)
    total_cov = sum(g["covered_fraction"] * g["exon_bp"] for g in per_gene) / float(total_bp)
    for g in per_gene:
        if g["covered_fraction"] >= args.present:
            g["state"] = "present"
        elif g["covered_fraction"] < args.absent:
            g["state"] = "absent"
        else:
            g["state"] = "partial"
    present = [g for g in per_gene if g["state"] == "present"]
    lost = [g for g in per_gene if g["state"] != "present"]
    medians = [g["median_depth"] for g in present]
    mito_depth = statistics.median(medians) if medians else 0.0

    notes = []
    if total_cov < args.absent:
        status = "rho0"
        reason = ("no read covers the mitochondrial genes (%.0f%% of %d exon bp at %.1fx or "
                  "more): no mtDNA" % (100 * total_cov, total_bp, floor))
    elif not lost:
        status = "rho+"
        reason = ("all %d respiratory genes covered (%.0f%% of %d exon bp at %.1fx or more)"
                  % (len(per_gene), 100 * total_cov, total_bp, floor))
    elif any(g["state"] == "absent" for g in lost):
        status = "rho-"
        def listed(group, limit=8):
            text = ", ".join("%s (%.0f%% covered)" % (g["gene"], 100 * g["covered_fraction"])
                             for g in group[:limit])
            return text + (", and %d more" % (len(group) - limit) if len(group) > limit else "")
        reason = ("%d of %d respiratory genes are missing or partly missing, the mtDNA is "
                  "reduced to a fragment; %s"
                  % (len(lost), len(per_gene),
                     "kept: " + listed(present) if present and len(present) < len(lost)
                     else "lost: " + listed(lost)))
    else:
        # Every gene is there, some only in part. A rho- deletion removes whole
        # genes; a gene part covered is as likely a stretch where the strain's
        # mtDNA diverges from the reference (VAR1 varies in length between
        # strains) or where these reads cannot map.
        status = "uncertain"
        reason = ("every respiratory gene is covered, but %s only in part: a rho- breakpoint "
                  "inside the gene, or sequence of this strain the reads do not map to"
                  % ", ".join("%s (%.0f%% covered)" % (g["gene"], 100 * g["covered_fraction"])
                              for g in lost))
    if status == "rho+" and len(medians) >= 3 and mito_depth > 0:
        high = [g for g in present if g["median_depth"] > args.amplified * mito_depth]
        low = [g for g in present if g["median_depth"] < mito_depth / args.amplified]
        if high or low:
            notes.append(
                "uneven depth over the mtDNA (median %.0fx): %s. A fragment amplified in part "
                "of the cells is what rho- cells in a rho+ culture look like"
                % (mito_depth, "; ".join(
                    ["%s %.0fx" % (g["gene"], g["median_depth"]) for g in high + low])))
    glycerol = {"rho+": "growth expected", "rho-": "no growth",
                "rho0": "no growth", "uncertain": "unknown"}[status]
    return {
        "status": status, "glycerol": glycerol, "reason": reason,
        "floor_depth": round(floor, 2),
        "mito_median_depth": mito_depth,
        "covered_fraction": round(total_cov, 3),
        "genes": per_gene, "notes": notes,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genes", required=True, help="mito_genes.tsv of the database")
    parser.add_argument("--depth", required=True,
                        help="samtools depth -a over the mitochondrial chromosome")
    parser.add_argument("--idxstats", required=True,
                        help="samtools idxstats of the whole genome alignment")
    parser.add_argument("--read-length", type=float, required=True,
                        help="mean aligned read length, to turn read counts into depth")
    parser.add_argument("--seq-kind", default="reads", choices=["reads", "assembly"])
    parser.add_argument("--seq-molecule", default="DNA")
    parser.add_argument("--min-depth", type=float, default=4,
                        help="nuclear depth needed to say anything, and the lowest floor")
    parser.add_argument("--floor-fraction", type=float, default=0.1,
                        help="a gene base is covered at this fraction of the nuclear depth "
                             "or more (default 0.1)")
    parser.add_argument("--present", type=float, default=0.9,
                        help="covered fraction of the exons for a gene present (default 0.9)")
    parser.add_argument("--absent", type=float, default=0.05,
                        help="covered fraction under which a gene is absent (default 0.05)")
    parser.add_argument("--amplified", type=float, default=5.0,
                        help="fold from the median gene depth reported as uneven (default 5)")
    parser.add_argument("--out", required=True, help="JSON output")
    args = parser.parse_args()

    genes = load_tsv(args.genes)
    chrom = genes[0]["chrom"]
    stats = parse_idxstats(args.idxstats)
    nuclear = [(length, reads) for contig, (length, reads) in stats.items() if contig != chrom]
    nuclear_len = sum(length for length, _ in nuclear)
    nuclear_depth = (sum(reads for _, reads in nuclear) * args.read_length / nuclear_len
                     if nuclear_len else 0.0)
    mito_reads = stats.get(chrom, (0, 0))[1]
    base = {"chrom": chrom, "nuclear_depth": round(nuclear_depth, 2),
            "mito_reads": mito_reads, "read_length": args.read_length}

    if chrom not in stats:
        result = {"status": "not_assessed", "glycerol": "unknown",
                  "reason": "the reference genome has no mitochondrial chromosome %s" % chrom}
    elif args.seq_molecule.upper() == "RNA":
        result = {"status": "not_assessed", "glycerol": "unknown",
                  "reason": "RNA reads: mitochondrial depth follows expression, which "
                            "depends on the carbon source, not on the mtDNA"}
    elif args.seq_kind == "reads" and nuclear_depth < args.min_depth:
        result = {"status": "no_data", "glycerol": "unknown",
                  "reason": "nuclear depth %.1fx, below the %.0fx needed"
                            % (nuclear_depth, args.min_depth)}
    else:
        depth = parse_depth(args.depth, chrom)
        result = assess(genes, depth, nuclear_depth, args)
        if args.seq_kind == "assembly" and result["status"] == "rho0":
            # Assemblies are often published without their mtDNA.
            result.update({"status": "not_assessed", "glycerol": "unknown",
                           "reason": "no mitochondrial contig in the assembly; assemblies "
                                     "are often deposited without their mtDNA"})
        elif args.seq_kind == "reads" and result.get("mito_median_depth") and nuclear_depth:
            result["mtdna_copies_per_nuclear_genome"] = round(
                result["mito_median_depth"] / nuclear_depth, 1)
    result.update(base)
    result["limits"] = [
        "a point mutation in a mitochondrial gene (mit-) and a nuclear mutation that stops "
        "respiration (pet) are not seen from depth: rho+ means the mtDNA is complete, "
        "not that the strain respires"]
    with open(args.out, "w") as handle:
        json.dump(result, handle, indent=1)
    sys.stderr.write("mtDNA: %s, glycerol: %s; %s\n"
                     % (result["status"], result["glycerol"], result["reason"]))
    for note in result.get("notes", []):
        sys.stderr.write("  note: %s\n" % note)
    return 0


if __name__ == "__main__":
    sys.exit(main())
