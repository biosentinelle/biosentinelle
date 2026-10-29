#!/usr/bin/env python3
"""
Build the strain k-mer panel from a population of strains, each strain with its
own k-mers.

Inputs are the strain catalogue and, for its strains, either the population
matrix (one multi-sample genotype table, the 1011 Yeast Genomes) or one VCF
and callable BED per strain (from an assembly, or from reads). Output is a
site table, a k-mer table and a strain table.

What changed from merge_snp_panels.py, and why
----------------------------------------------
A probe k-mer is only valid for a strain if the k-1 bases around the SNP are
what that strain carries. With five panel strains the old rule, dropping any
site with another strain's variant in its window, cost a few per cent of the
sites. With a thousand strains, carrying together a variant every few bases,
it would drop almost all of them. So each strain gets its own k-mers: the
window is rebuilt from the strain's own alleles, and the site keeps one k-mer
per distinct haplotype, each labelled with the allele it carries at the
centre. A strain whose window cannot be rebuilt (an indel in it, a no-call,
more than two heterozygous neighbours) is recorded as unknown at that site,
never as reference.

Indels are sites too. The reference and alternative haplotypes are written
out with the allele between flank bases on each side, and each probe is the
k-mer of its haplotype centred on the allele: the junction for a deletion,
the inserted bases for an insertion. k stays the same, so short reads still
carry the probes. An indel that only shifts a tandem repeat leaves the
reference probe readable in the alternative allele and is dropped. Neighbour
indels are still not rebuilt into a strain's window: a strain with one there
is unknown at the site.

The site table is stored by exception. A strain not listed at a site carries
the reference allele there; the listed ones carry the alternative allele, are
heterozygous, or are unknown. A column per strain would be a thousand columns
by hundreds of thousands of rows.

Which sites are kept
--------------------
Two strains that are clonal relatives differ only at variants that are rare in
the population, often private to one of them. Those sites are what separates
them, so for every strain the rarest sites it carries are kept first, up to
--per-strain-sites. Common variants place a sample in its clade; a thinned,
evenly spread set of them is added on top (--common-sites).

Positions absent from the matrix are taken as reference in every matrix
strain: the matrix lists every position variable in the population, so a
position missing from it is one no isolate was called non-reference at. That
is an assumption about coverage there, stated in the report.
"""

import argparse
import bisect
import csv
import gzip
import itertools
import os
import sys
from collections import OrderedDict

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
SHORT = 50  # records with a reference allele longer than this are kept apart


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def canonical(kmer):
    reverse = revcomp(kmer)
    return kmer if kmer <= reverse else reverse


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path):
    names, sequences, name, chunks = [], {}, None, []
    with opener(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
                names.append(name)
            else:
                chunks.append(line.upper())
    if name is not None:
        sequences[name] = "".join(chunks)
    return names, sequences


def low_complexity(kmer, max_homopolymer, min_distinct):
    if len(set(kmer)) < min_distinct:
        return True
    run, previous = 1, ""
    for base in kmer:
        run = run + 1 if base == previous else 1
        if run > max_homopolymer:
            return True
        previous = base
    return False


def centre_variant(ref, alts, max_indel):
    """True for a record that can be the centre of a site: biallelic, and a
    SNV or an indel of at most max_indel bases (0 for SNVs only)."""
    if len(alts) != 1:
        return False
    alt = alts[0]
    if not ref or not alt or ref == alt or set(ref + alt) - set("ACGT"):
        return False
    if len(ref) == len(alt):
        return len(ref) == 1
    return abs(len(ref) - len(alt)) <= max_indel and max(len(ref), len(alt)) <= SHORT


def left_align(sequence, pos, ref, alt):
    """Left-aligned, parsimonious form of a variant (pos 1-based), so that one
    indel has one key whichever position its source wrote it at. A SNV comes
    back unchanged."""
    if len(ref) == len(alt):
        return pos, ref, alt
    while True:
        changed = False
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt, changed = ref[:-1], alt[:-1], True
        if (not ref or not alt) and pos > 1:
            pos -= 1
            base = sequence[pos - 1]
            ref, alt, changed = base + ref, base + alt, True
        if not changed:
            break
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def centred_kmer(haplotype, flank, allele_length, k):
    """The k-mer of a haplotype centred on its allele, which starts at offset
    flank. For a SNV it is the whole haplotype; for an indel allele longer
    than one base it is centred on the allele's middle base."""
    start = allele_length // 2
    return haplotype[start:start + k]


def parse_gt(gt):
    """'0/1' -> (alleles tuple, missing flag)."""
    if "." in gt:
        return (), True
    return tuple(int(a) for a in gt.replace("|", "/").split("/") if a != ""), False


class Records:
    """Variant records of one contig for one source, with window lookup."""

    def __init__(self):
        self.short_pos, self.short, self.long = [], [], []

    def add(self, record):
        # record: (pos, ref, alts, genotypes) ; genotypes is either a raw
        # 'SAMPLE:GT,...' string (matrix) or a dict {strain_index: gt}
        if len(record[1]) > SHORT:
            self.long.append(record)
        else:
            self.short_pos.append(record[0])
            self.short.append(record)

    def overlapping(self, start, end):
        """Records whose reference span meets [start, end] (1-based, inclusive)."""
        out = []
        low = bisect.bisect_left(self.short_pos, start - SHORT)
        high = bisect.bisect_right(self.short_pos, end)
        for index in range(low, high):
            record = self.short[index]
            if record[0] + len(record[1]) - 1 >= start:
                out.append((index, record))
        for record in self.long:
            if record[0] <= end and record[0] + len(record[1]) - 1 >= start:
                out.append((-1, record))
        return out


def load_bed(path):
    regions = {}
    with opener(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            regions.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    for contig in regions:
        regions[contig].sort()
        merged = []
        for start, end in regions[contig]:
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        regions[contig] = (merged, [s for s, _ in merged])
    return regions


def covered(regions, chrom, start, end):
    """True when the 1-based inclusive span [start, end] lies in one interval."""
    entry = regions.get(chrom)
    if not entry:
        return False
    intervals, starts = entry
    index = bisect.bisect_right(starts, start - 1) - 1
    return index >= 0 and intervals[index][0] <= start - 1 and intervals[index][1] >= end


def load_single_vcf(path, min_qual):
    """Per contig list of (pos, ref, alt, gt) from a one-sample VCF."""
    out = {}
    with opener(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, _, ref, alts = fields[:5]
            qual = fields[5] if len(fields) > 5 else "."
            try:
                if qual != "." and float(qual) < min_qual:
                    continue
            except ValueError:
                pass
            gt = "1/1"
            if len(fields) > 9:
                keys = fields[8].split(":")
                if "GT" in keys:
                    gt = fields[9].split(":")[keys.index("GT")]
            alleles, missing = parse_gt(gt)
            if missing or not any(alleles):
                continue
            out.setdefault(chrom, []).append((int(pos), ref.upper(), alts.upper().split(","), gt))
    for chrom in out:
        out[chrom].sort()
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--strains", required=True, help="strain catalogue TSV")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--matrix", help="sparse genotype table: CHROM POS REF ALT "
                                         "SAMPLE:GT,... listing every sample not 0/0")
    parser.add_argument("--matrix-contig-prefix", default="chromosome",
                        help="matrix contig N is the N-th reference record")
    parser.add_argument("--single", action="append", default=[],
                        metavar="STRAIN=VCF=BED",
                        help="one-sample variants and callable regions, repeatable")
    parser.add_argument("--min-nocall-ref-reads", type=int, default=2,
                        help="a matrix no-call is read as reference when at least this "
                             "many reads carry the reference allele and none another")
    parser.add_argument("--min-qual", type=float, default=30.0)
    parser.add_argument("-k", type=int, default=21)
    parser.add_argument("--rare-max-carriers", type=int, default=20,
                        help="a site carried by at most this many strains is rare")
    parser.add_argument("--per-strain-sites", type=int, default=300,
                        help="rare sites kept per strain, rarest first")
    parser.add_argument("--common-sites", type=int, default=30000,
                        help="common sites kept, evenly spread")
    parser.add_argument("--max-unknown-fraction", type=float, default=0.10,
                        help="drop sites unknown in more than this fraction of strains")
    parser.add_argument("--max-indel", type=int, default=20,
                        help="longest indel, in bases gained or lost, kept as a site; "
                             "0 for SNVs only")
    parser.add_argument("--trace-strain", metavar="STRAIN",
                        help="write, for every variant this strain carries, the step "
                             "that kept or dropped it, to --trace-out")
    parser.add_argument("--trace-out", default="panel_trace.tsv")
    parser.add_argument("--max-het-neighbours", type=int, default=2)
    parser.add_argument("--max-homopolymer", type=int, default=6)
    parser.add_argument("--min-distinct-bases", type=int, default=3)
    parser.add_argument("--out-sites", required=True)
    parser.add_argument("--out-kmers", required=True)
    parser.add_argument("--out-strains", required=True)
    args = parser.parse_args()

    if args.k % 2 == 0:
        parser.error("k must be odd")
    flank = args.k // 2

    with open(args.strains, newline="", encoding="utf-8") as handle:
        catalogue = list(csv.DictReader(handle, delimiter="\t"))
    singles_given = {}
    for item in args.single:
        parts = item.split("=")
        if len(parts) != 3:
            parser.error("--single expects STRAIN=VCF=BED, got %r" % item)
        singles_given[parts[0]] = (parts[1], parts[2])

    strains, code_index, single_index = [], {}, {}
    for row in catalogue:
        if row["source_type"] == "matrix":
            if not args.matrix:
                continue
            code_index[row["source"]] = len(strains)
        elif row["strain"] in singles_given:
            single_index[row["strain"]] = len(strains)
        else:
            sys.stderr.write("no variants given for %s (%s), left out\n"
                             % (row["strain"], row["source_type"]))
            continue
        strains.append(row)
    n_strains = len(strains)
    if not n_strains:
        sys.exit("no strain has variants")
    max_unknown = int(args.max_unknown_fraction * n_strains)
    trace, trace_out = None, []
    if args.trace_strain:
        trace = next((i for i, row in enumerate(strains)
                      if row["strain"] == args.trace_strain), None)
        if trace is None:
            parser.error("--trace-strain %s is not in the catalogue" % args.trace_strain)
    sys.stderr.write("%d strains: %d from the matrix, %d from their own VCF\n"
                     % (n_strains, len(code_index), len(single_index)))

    names, reference = read_fasta(args.reference)
    contig_of = {"%s%d" % (args.matrix_contig_prefix, i + 1): name
                 for i, name in enumerate(names)}

    singles = {}
    for strain, index in single_index.items():
        vcf, bed = singles_given[strain]
        singles[index] = (load_single_vcf(vcf, args.min_qual), load_bed(bed))

    # ------------------------------------------------------------------ matrix
    def matrix_contigs():
        """Yield (reference contig, Records) one contig at a time."""
        if not args.matrix:
            return
        current, records = None, None
        with opener(args.matrix) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                contig = contig_of.get(fields[0], fields[0])
                if contig != current:
                    if records is not None:
                        yield current, records
                    current, records = contig, Records()
                records.add((int(fields[1]), fields[2].upper(),
                             fields[3].upper().split(","), fields[4]))
        if records is not None:
            yield current, records

    parsed_cache = OrderedDict()

    def genotypes(contig_key, index, record):
        """Matrix genotypes as {strain_index: gt}, parsed once per record."""
        if isinstance(record[3], dict):
            return record[3]
        key = (contig_key, index, record[0])
        if key in parsed_cache:
            return parsed_cache[key]
        # Items are SAMPLE:GT:AD. The matrix holds many no-calls over reads
        # that all carry the reference base ('./.:14,0'), left by its
        # joint-genotyping thresholds. Read as unknown they would blank out
        # a tenth of every window, so a no-call with enough reference reads
        # and no alternative read is taken as reference.
        # The AD values are comma separated like the items themselves, so a
        # token without ':' continues the previous item's AD.
        items = []
        for token in record[3].split(","):
            if ":" in token:
                fields = token.split(":")
                items.append([fields[0], fields[1], fields[2:3]])
            elif token and items:
                items[-1][2].append(token)
        out = {}
        for sample, gt, depths in items:
            strain = code_index.get(sample)
            if strain is None:
                continue
            if "." in gt:
                depths = [int(x) for x in depths if x.isdigit()]
                if (depths and depths[0] >= args.min_nocall_ref_reads
                        and not any(depths[1:])):
                    continue
            out[strain] = gt
        parsed_cache[key] = out
        if len(parsed_cache) > 4000:
            parsed_cache.popitem(last=False)
        return out

    stats = dict.fromkeys(["candidates", "edge", "ref_mismatch", "low_complexity",
                           "indel_in_repeat", "ambiguous", "too_unknown",
                           "uninformative", "thinned_common"], 0)
    sites = []  # kept candidates, before uniqueness and selection

    # Common sites are thinned by distance while they are read: keeping
    # every window of a million common candidates in memory to thin them
    # afterwards would cost gigabytes. About three times the requested number
    # are kept, for the later filters to choose from.
    genome_length = sum(len(reference[name]) for name in names)
    common_spacing = max(1, genome_length // (3 * max(args.common_sites, 1)))

    def contig_stream():
        done = set()
        for contig, records in matrix_contigs():
            done.add(contig)
            yield contig, records
        for contig in names:
            if contig not in done and any(contig in v for v, _ in singles.values()):
                yield contig, Records()

    for contig, records in contig_stream():
        sequence = reference.get(contig)
        if sequence is None:
            sys.stderr.write("contig %s is not in the reference, skipped\n" % contig)
            continue
        last_common = -common_spacing

        # Per single strain: its variants on this contig as Records with a
        # genotype dict, so matrix and single strains go through one path.
        single_records = {}
        for index, (variants, _) in singles.items():
            recs = Records()
            for pos, ref, alts, gt in variants.get(contig, []):
                recs.add((pos, ref, alts, {index: gt}))
            single_records[index] = recs

        # Candidate centres: biallelic SNVs and short indels from the matrix,
        # plus single-strain ones at positions the matrix does not list.
        # An indel is keyed by its left-aligned form: the matrix and an
        # assembly VCF can write one event at two positions, which would make
        # two sites with the same probes, both then lost as a collision.
        centres = {}
        for index, record in enumerate(records.short):
            if centre_variant(record[1], record[2], args.max_indel):
                norm = left_align(sequence, record[0], record[1], record[2][0])
                centres.setdefault(norm[0], []).append(("m", index, record, norm))
        for strain, recs in single_records.items():
            for record in recs.short:
                if centre_variant(record[1], record[2], args.max_indel):
                    norm = left_align(sequence, record[0], record[1], record[2][0])
                    centres.setdefault(norm[0], []).append(("s", strain, record, norm))

        for pos in sorted(centres):
            entries = centres[pos]
            alleles_here = {(entry[3][1], entry[3][2]) for entry in entries}
            stats["candidates"] += 1
            traced = trace is not None and any(
                trace in (genotypes(contig, owner, record) if kind == "m" else record[3])
                for kind, owner, record, _ in entries)

            def reject(reason):
                stats[reason] += 1
                if traced:
                    trace_out.append((contig, pos, reason))

            if len(alleles_here) != 1:
                if traced:
                    trace_out.append((contig, pos, "multiallelic"))
                continue  # two different alternative alleles: not a biallelic site
            ref, alt = alleles_here.pop()

            def same_event(record):
                if len(record[2]) != 1:
                    return False
                if len(record[1]) == len(record[2][0]):
                    return record[0] == pos and record[1] == ref and record[2][0] == alt
                return left_align(sequence, record[0], record[1], record[2][0]) == (pos, ref, alt)
            indel = len(ref) != len(alt)
            # The reference segment: the centre allele and flank bases on each
            # side. For a SNV it is the k-mer window itself; for an indel it is
            # longer, and the probes are cut from it once the allele is in.
            last = pos + len(ref) - 1
            start, end = pos - 1 - flank, last + flank
            if start < 0 or end > len(sequence):
                reject("edge")
                continue
            window = sequence[start:end]
            if window[flank:flank + len(ref)] != ref:
                reject("ref_mismatch")
                continue
            ref_probe = centred_kmer(window, flank, len(ref), args.k)
            alt_haplotype = window[:flank] + alt + window[flank + len(ref):]
            alt_probe = centred_kmer(alt_haplotype, flank, len(alt), args.k)
            if "N" in window or any(low_complexity(probe, args.max_homopolymer,
                                                   args.min_distinct_bases)
                                    for probe in (ref_probe, alt_probe)):
                reject("low_complexity")
                continue
            # An indel in a tandem repeat shifts the repeat rather than
            # changing its sequence: the reference probe can then still be
            # read in the alternative allele, and the site cannot tell them
            # apart. The alternative probe is checked against the whole
            # reference later, which covers the other direction.
            if indel and (ref_probe in alt_haplotype or revcomp(ref_probe) in alt_haplotype):
                reject("indel_in_repeat")
                continue

            # Centre states from the centre records.
            state = {}
            for kind, owner, record, _ in entries:
                gts = genotypes(contig, owner, record) if kind == "m" else record[3]
                for strain, gt in gts.items():
                    alleles, missing = parse_gt(gt)
                    if missing:
                        state[strain] = "U"
                    elif alleles and all(a == 1 for a in alleles):
                        state[strain] = "A"
                    elif 1 in alleles:
                        state[strain] = "H"
            carriers = sum(1 for s in state.values() if s in "AH")
            rare = carriers <= args.rare_max_carriers
            if not rare:
                if pos - last_common < common_spacing:
                    reject("thinned_common")
                    continue
                last_common = pos

            # Single strains not carrying the variant: reference if callable
            # over the whole window, unknown otherwise.
            for strain, (variants, callable_regions) in singles.items():
                if strain in state:
                    continue
                if not covered(callable_regions, contig, pos - flank, last + flank):
                    state[strain] = "U"

            # Neighbours: every record meeting the window other than the
            # centre variant itself, from the matrix and from each single strain.
            neighbours = []
            for index, record in records.overlapping(pos - flank, last + flank):
                if same_event(record):
                    continue
                neighbours.append((record, genotypes(contig, index, record)))
            for strain, recs in single_records.items():
                for _, record in recs.overlapping(pos - flank, last + flank):
                    if same_event(record):
                        continue
                    neighbours.append((record, record[3]))

            # Per strain substitutions in the window.
            subs = {}  # strain -> list of (offset, ref_seq, [allele seqs], het)
            for record, gts in neighbours:
                rpos, rref, ralts = record[0], record[1], record[2]
                for strain, gt in gts.items():
                    if state.get(strain) == "U":
                        continue
                    alleles, missing = parse_gt(gt)
                    if missing:
                        state[strain] = "U"
                        continue
                    nonref = sorted({a for a in alleles if a != 0})
                    if not nonref:
                        continue
                    seqs = []
                    for a in nonref:
                        allele = ralts[a - 1] if a - 1 < len(ralts) else "*"
                        if allele == "*" or len(allele) != len(rref):
                            seqs = None
                            break
                        seqs.append(allele)
                    offset = rpos - (pos - flank)
                    span = range(offset, offset + len(rref))
                    if seqs is None or any(flank + i in span for i in range(len(ref))):
                        state[strain] = "U"
                        continue
                    het = 0 in alleles or len(nonref) > 1
                    options = ([rref] if 0 in alleles else []) + seqs
                    subs.setdefault(strain, []).append((offset, rref, options, het))

            unknown = sum(1 for s in state.values() if s == "U")
            if unknown > max_unknown:
                reject("too_unknown")
                continue

            # k-mers: the generic pair, then one per strain haplotype that
            # differs from it.
            kmers = {canonical(ref_probe): "R", canonical(alt_probe): "A"}
            ambiguous = False
            for strain, items in subs.items():
                if state.get(strain) == "U":
                    continue
                if sum(1 for item in items if item[3]) > args.max_het_neighbours:
                    state[strain] = "U"
                    continue
                centre_alleles = {"A": [alt], "H": [ref, alt]}.get(state.get(strain, "R"), [ref])
                for combo in itertools.product(*[item[2] for item in items]):
                    seq = list(window)
                    for (offset, rref, _, _), allele in zip(items, combo):
                        for i, base in enumerate(allele):
                            if 0 <= offset + i < len(seq):
                                seq[offset + i] = base
                    for centre in centre_alleles:
                        haplotype = "".join(seq[:flank]) + centre + "".join(seq[flank + len(ref):])
                        kmer = canonical(centred_kmer(haplotype, flank, len(centre), args.k))
                        label = "R" if centre == ref else "A"
                        if kmers.get(kmer, label) != label:
                            ambiguous = True
                        kmers[kmer] = label
            if ambiguous:
                # one strain's haplotype k-mer is another's of the other
                # allele: the site cannot be read
                reject("ambiguous")
                continue

            alt_strains = sorted(s for s, v in state.items() if v == "A")
            het_strains = sorted(s for s, v in state.items() if v == "H")
            unknown_strains = sorted(s for s, v in state.items() if v == "U")
            if not alt_strains and not het_strains:
                reject("uninformative")
                continue
            if len(alt_strains) + len(unknown_strains) >= n_strains:
                reject("uninformative")
                continue
            sites.append({"chrom": contig, "pos": pos, "ref": ref, "alt": alt,
                          "alt_strains": alt_strains, "het_strains": het_strains,
                          "unknown_strains": unknown_strains, "kmers": kmers,
                          "ref_kmer": canonical(ref_probe), "indel": indel,
                          "rare": rare, "carriers": len(alt_strains) + len(het_strains),
                          "traced": traced})
        parsed_cache.clear()
        sys.stderr.write("%s: %d candidate sites so far\n" % (contig, len(sites)))

    # ----------------------------------------------------------- uniqueness
    # The reference-window k-mer must occur once in the reference, every
    # other k-mer of the site never: a probe found elsewhere in the genome is
    # evidence about another locus.
    wanted = {}
    for site in sites:
        for kmer in site["kmers"]:
            wanted[kmer] = 0
    k = args.k
    for name in names:
        sequence = reference[name]
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            if kmer in wanted:
                wanted[kmer] += 1
            else:
                reverse = revcomp(kmer)
                if reverse in wanted:
                    wanted[reverse] += 1
    owner, collided = {}, set()
    survivors = []
    dropped_repeat = 0
    for index, site in enumerate(sites):
        ref_kmer = site["ref_kmer"]
        if wanted.get(ref_kmer, 0) != 1 or any(
                wanted.get(kmer, 0) for kmer in site["kmers"] if kmer != ref_kmer):
            dropped_repeat += 1
            if site["traced"]:
                trace_out.append((site["chrom"], site["pos"], "repeat"))
            continue
        for kmer in site["kmers"]:
            owner.setdefault(kmer, set()).add(index)
        survivors.append(index)
    # A short deletion's reference probe is centred on a deleted base, so it
    # is also the reference probe of a SNV at that base. The SNV is the
    # cleaner site and is kept; only the indels sharing its k-mer go. Any
    # other shared k-mer costs every site holding it.
    for holders in owner.values():
        if len(holders) < 2:
            continue
        snvs = [index for index in holders if not sites[index]["indel"]]
        collided.update(index for index in holders
                        if len(snvs) != 1 or sites[index]["indel"])
    survivors = [index for index in survivors if index not in collided]
    trace_out.extend((sites[i]["chrom"], sites[i]["pos"], "kmer_collision")
                     for i in collided if sites[i]["traced"])

    # ------------------------------------------------------------ selection
    def spread(site):
        return (site["pos"] * 2654435761) & 0xFFFFFFFF

    keep = set()
    by_strain = {}
    for index in survivors:
        site = sites[index]
        if site["rare"]:
            for strain in site["alt_strains"] + site["het_strains"]:
                by_strain.setdefault(strain, []).append(index)
    for strain, indices in by_strain.items():
        indices.sort(key=lambda i: (sites[i]["carriers"], spread(sites[i])))
        keep.update(indices[:args.per_strain_sites])
    common = [index for index in survivors if not sites[index]["rare"]]
    if len(common) > args.common_sites:
        step = len(common) / float(args.common_sites)
        common = [common[int(i * step)] for i in range(args.common_sites)]
    keep.update(common)
    kept = sorted(keep, key=lambda i: (names.index(sites[i]["chrom"]), sites[i]["pos"]))
    trace_out.extend((sites[i]["chrom"], sites[i]["pos"],
                      "kept" if i in keep else "not_selected")
                     for i in survivors if sites[i]["traced"])
    if trace is not None:
        with open(args.trace_out, "w", newline="") as handle:
            handle.write("chrom\tpos\tfate\n")
            for row in sorted(trace_out, key=lambda r: (names.index(r[0]) if r[0] in names
                                                        else -1, r[1])):
                handle.write("%s\t%d\t%s\n" % row)

    # ---------------------------------------------------------------- write
    for path in (args.out_sites, args.out_kmers, args.out_strains):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(args.out_strains, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["index", "strain", "display_name", "source_type", "source",
                         "clade", "ploidy", "zygosity", "origin", "n_sites_carried"])
        carried = [0] * n_strains
        for index in kept:
            for strain in sites[index]["alt_strains"] + sites[index]["het_strains"]:
                carried[strain] += 1
        for i, row in enumerate(strains):
            writer.writerow([i, row["strain"], row["display_name"], row["source_type"],
                             row["source"], row["clade"], row["ploidy"], row["zygosity"],
                             row["origin"], carried[i]])
    with open(args.out_sites, "w", newline="") as sites_handle, \
            open(args.out_kmers, "w", newline="") as kmers_handle:
        sites_writer = csv.writer(sites_handle, delimiter="\t", lineterminator="\n")
        kmers_writer = csv.writer(kmers_handle, delimiter="\t", lineterminator="\n")
        sites_writer.writerow(["site_id", "chrom", "pos", "ref", "alt", "n_alt", "n_het",
                               "n_unknown", "alt_strains", "het_strains",
                               "unknown_strains"])
        kmers_writer.writerow(["kmer", "site_id", "allele"])
        for number, index in enumerate(kept, start=1):
            site = sites[index]
            site_id = "s%07d" % number
            sites_writer.writerow([site_id, site["chrom"], site["pos"], site["ref"],
                                   site["alt"], len(site["alt_strains"]),
                                   len(site["het_strains"]), len(site["unknown_strains"]),
                                   ",".join(map(str, site["alt_strains"])),
                                   ",".join(map(str, site["het_strains"])),
                                   ",".join(map(str, site["unknown_strains"]))])
            for kmer, allele in sorted(site["kmers"].items()):
                kmers_writer.writerow([kmer, site_id, allele])

    # Strains no kept site separates: identical expected alleles everywhere
    # they are both known. Reported as groups, since the panel cannot tell
    # their members apart however deep the sample is sequenced.
    signature = [0] * n_strains
    for index in kept:
        site = sites[index]
        token = hash((site["chrom"], site["pos"]))
        for strain in site["alt_strains"]:
            signature[strain] ^= token
        for strain in site["het_strains"]:
            signature[strain] ^= token * 3
    groups = {}
    for i, value in enumerate(signature):
        groups.setdefault(value, []).append(strains[i]["strain"])
    tied = [members for members in groups.values() if len(members) > 1]

    n_kmers = sum(len(sites[i]["kmers"]) for i in kept)
    print("strains: %d" % n_strains)
    print("sites kept: %d (rare %d, common %d; SNV %d, indel %d), k-mers: %d, k=%d" % (
        len(kept), sum(1 for i in kept if sites[i]["rare"]),
        sum(1 for i in kept if not sites[i]["rare"]),
        sum(1 for i in kept if not sites[i]["indel"]),
        sum(1 for i in kept if sites[i]["indel"]), n_kmers, args.k))
    for key in sorted(stats):
        print("  %-16s %d" % (key, stats[key]))
    print("  %-16s %d" % ("repeat", dropped_repeat))
    print("  %-16s %d" % ("kmer_collision", len(collided)))
    if tied:
        print("strain groups with identical expected alleles at every kept site:")
        for members in sorted(tied, key=len, reverse=True):
            print("  " + ", ".join(sorted(members)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
