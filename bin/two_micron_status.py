#!/usr/bin/env python3
"""
2-micron plasmid status of a sample: present (cir+) or absent (cir0), and its
copy number.

The 2-micron plasmid (6,318 bp) lives in the nucleus of most Saccharomyces
cerevisiae strains, at tens of copies per haploid genome (cir+). Its four genes
keep it there: FLP1 amplifies the copy number by recombining the two inverted
repeats (IR1, IR2), REP1 and REP2 partition it between mother and daughter
cells with the STB locus, and RAF1 relieves the repression of FLP1.
The cell does not need it, and some laboratory strains have lost it (cir0). It
matters for the 2-micron vectors (YEp): those that carry only the origin and
STB of the plasmid rely on the REP1 and REP2 proteins of the native plasmid to
stay at high copy.

The plasmid is not part of the reference genome. Its sequence (database
two_micron/two_micron.fasta) is added to the genome the reads are aligned to,
and the depth over it is set against the nuclear depth:

    present   every feature covered (FLP1, REP1, REP2, RAF1, IR1, IR2)    cir+
    absent    no read over the plasmid                                    cir0
    partial   part of the plasmid only: a 2-micron vector without the
              plasmid itself, or a plasmid deleted in part or diverging
              from the reference

A base is covered at a fraction of the nuclear depth or more, as for the
mitochondrial genome: a plasmid in tens of copies covers its bases far above
the nuclear depth, so its absence is clear even at a low nuclear depth. The
copy number is the median depth over the covered bases divided by the nuclear
depth (reads of DNA only).

An assembly has one copy of each sequence: the plasmid is there when a contig
covers it; a missing plasmid is not assessed, assemblies being often deposited
without it. RNA reads show the plasmid when its transcripts are there; without
them nothing is said, since depth follows expression.

Not seen: a 2-micron vector (YEp) in a cir+ strain, whose reads only add to the
depth over the part of the plasmid it carries.
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
    """0-based position -> depth over chrom."""
    out = {}
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 3 and fields[0] == chrom:
                out[int(fields[1]) - 1] = int(fields[2])
    return out


def positions(row):
    out = []
    for item in row["segments"].split(";"):
        a, b = item.split("-")
        out.extend(range(int(a), int(b)))
    return out


def listed(group):
    return ", ".join("%s (%.0f%%)" % (g["feature"], 100 * g["covered_fraction"]) for g in group)


def assess(features, depth, length, nuclear_depth, args):
    """Per-feature coverage and the status. Returns the result dict."""
    if args.seq_kind == "assembly":
        floor = 1   # one copy of everything: a base is covered or it is not
    else:
        floor = max(args.min_depth, args.floor_fraction * nuclear_depth)
    per_feature = []
    for row in features:
        values = [depth.get(p, 0) for p in positions(row)]
        covered = sum(1 for v in values if v >= floor)
        fraction = covered / float(len(values)) if values else 0.0
        per_feature.append({
            "feature": row["feature"], "role": row["role"], "length_bp": len(values),
            "covered_fraction": round(fraction, 3),
            "median_depth": statistics.median(values) if values else 0,
            "state": "present" if fraction >= args.present else
                     "absent" if fraction < args.absent else "partial",
        })
    covered_depths = [depth.get(p, 0) for p in range(length) if depth.get(p, 0) >= floor]
    total_cov = len(covered_depths) / float(length)
    median_depth = statistics.median(covered_depths) if covered_depths else 0

    present = [g for g in per_feature if g["state"] == "present"]
    absent = [g for g in per_feature if g["state"] == "absent"]
    partial = [g for g in per_feature if g["state"] == "partial"]
    if total_cov < args.absent:
        status = "absent"
        reason = ("no read covers the plasmid (%.0f%% of its %d bp at %.1fx or more)"
                  % (100 * total_cov, length, floor))
        if args.seq_kind == "reads" and nuclear_depth:
            reason += (", where a cir+ strain gives tens of times the nuclear depth of %.1fx"
                       % nuclear_depth)
    elif not absent and not partial:
        status = "present"
        reason = ("every feature of the plasmid is covered (%s; %.0f%% of its %d bp at %.1fx "
                  "or more)" % (", ".join(g["feature"] for g in per_feature), 100 * total_cov,
                                length, floor))
    elif absent:
        status = "partial"
        reason = ("only part of the plasmid is covered (%.0f%% of its %d bp at %.1fx or more); "
                  "covered: %s; missing: %s. A 2-micron vector (YEp), which carries the "
                  "origin and STB region of the plasmid but not all its genes, gives this in a "
                  "strain without the plasmid itself (cir0); a plasmid deleted in part, or "
                  "diverging from the reference over the parts not covered, looks the same"
                  % (100 * total_cov, length, floor,
                     listed([g for g in per_feature if g["state"] != "absent"]) or "none",
                     listed(absent)))
    else:
        # Every feature is there, some only in part: the plasmid of this strain
        # diverges from the reference over a stretch, or a deletion ends there.
        status = "partial"
        reason = ("every feature of the plasmid is covered, but %s only in part: a plasmid "
                  "diverging from the reference over these stretches, as the plasmid of "
                  "another strain can, or a deletion ending inside the feature"
                  % listed(partial))
    return {
        "status": status, "cir": {"present": "cir+", "absent": "cir0"}.get(status),
        "reason": reason, "floor_depth": round(floor, 2),
        "covered_fraction": round(total_cov, 3), "median_depth": median_depth,
        "features": per_feature, "notes": [],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--features", required=True,
                        help="two_micron_features.tsv of the database")
    parser.add_argument("--depth", required=True,
                        help="samtools depth -a over the 2-micron sequence of the alignment")
    parser.add_argument("--idxstats", required=True,
                        help="samtools idxstats of the reference genome sequences of the "
                             "alignment, for the nuclear depth")
    parser.add_argument("--mito-chrom", default="",
                        help="mitochondrial chromosome of the reference, left out of the "
                             "nuclear depth")
    parser.add_argument("--read-length", type=float, required=True,
                        help="mean aligned read length, to turn read counts into depth")
    parser.add_argument("--seq-kind", default="reads", choices=["reads", "assembly"])
    parser.add_argument("--seq-molecule", default="DNA")
    parser.add_argument("--min-depth", type=float, default=4,
                        help="nuclear depth needed to say anything, and the lowest floor")
    parser.add_argument("--floor-fraction", type=float, default=0.1,
                        help="a base is covered at this fraction of the nuclear depth or more "
                             "(default 0.1)")
    parser.add_argument("--present", type=float, default=0.9,
                        help="covered fraction of a feature for it to be present (default 0.9)")
    parser.add_argument("--absent", type=float, default=0.05,
                        help="covered fraction under which a feature, or the plasmid, is "
                             "absent (default 0.05)")
    parser.add_argument("--low-copies", type=float, default=5.0,
                        help="copies per nuclear genome under which a plasmid present is "
                             "noted as rare in the culture (default 5)")
    parser.add_argument("--out", required=True, help="JSON output")
    args = parser.parse_args()

    features = load_tsv(args.features)
    chrom = features[0]["chrom"]
    stats = parse_idxstats(args.idxstats)
    nuclear = [(length, reads) for contig, (length, reads) in stats.items()
               if contig not in (chrom, args.mito_chrom)]
    nuclear_len = sum(length for length, _ in nuclear)
    nuclear_depth = (sum(reads for _, reads in nuclear) * args.read_length / nuclear_len
                     if nuclear_len else 0.0)
    # samtools depth -aa: every base of the plasmid, read or not
    depth = parse_depth(args.depth, chrom)
    length = max(depth) + 1 if depth else 0
    base = {"chrom": chrom, "length_bp": length, "nuclear_depth": round(nuclear_depth, 2),
            "read_length": args.read_length}
    rna = args.seq_molecule.upper() == "RNA"

    if not length:
        result = {"status": "not_assessed", "cir": None,
                  "reason": "the 2-micron sequence %s was not part of the alignment" % chrom}
    elif args.seq_kind == "reads" and not rna and nuclear_depth < args.min_depth:
        result = {"status": "no_data", "cir": None,
                  "reason": "nuclear depth %.1fx, below the %.0fx needed"
                            % (nuclear_depth, args.min_depth)}
    else:
        result = assess(features, depth, length, nuclear_depth, args)
        if rna:
            if result["status"] == "absent":
                result.update({"status": "not_assessed", "cir": None,
                               "reason": "RNA reads: no 2-micron transcript seen, which does "
                                         "not show the plasmid absent, since depth follows "
                                         "expression"})
            else:
                result.update({"status": "present", "cir": None,
                               "reason": "2-micron transcripts in the RNA reads cover %.0f%% "
                                         "of the plasmid: it is there, or a 2-micron vector; "
                                         "depth follows expression, so the copy number is not "
                                         "measured" % (100 * result["covered_fraction"])})
        elif args.seq_kind == "assembly":
            if result["status"] == "absent":
                result.update({"status": "not_assessed", "cir": None,
                               "reason": "no 2-micron contig in the assembly; assemblies are "
                                         "often deposited without the plasmid"})
        elif result["median_depth"] and nuclear_depth:
            copies = round(result["median_depth"] / nuclear_depth, 1)
            result["copies_per_nuclear_genome"] = copies
            if result["status"] == "present" and copies < args.low_copies:
                result["notes"].append(
                    "%.1f copies per nuclear genome, far below the tens of copies of a cir+ strain: "
                    "the plasmid may be in part of the cells only (a culture losing it, or a "
                    "mix of strains), or its reads may come from another sample" % copies)
    result.update(base)
    result["limits"] = [
        "a 2-micron vector (YEp) carries part of the plasmid (its origin and STB region): in a "
        "cir+ strain its reads only add to the depth over that part, so it is not told apart "
        "from the plasmid; in a cir0 strain the plasmid shows as partial",
        "the plasmid of a strain can diverge from the reference sequence (A form of strain "
        "A364A, RefSeq NC_001398.1): a stretch that diverges is covered in part or not at all"]
    with open(args.out, "w") as handle:
        json.dump(result, handle, indent=1)
    sys.stderr.write("2-micron plasmid: %s%s; %s\n" % (
        result["status"], " (%s)" % result["cir"] if result.get("cir") else "",
        result["reason"]))
    if result.get("copies_per_nuclear_genome") is not None:
        sys.stderr.write("  %.1f copies per nuclear genome\n" % result["copies_per_nuclear_genome"])
    for note in result.get("notes", []):
        sys.stderr.write("  note: %s\n" % note)
    return 0


if __name__ == "__main__":
    sys.exit(main())
