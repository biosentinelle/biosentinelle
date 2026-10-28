#!/usr/bin/env python3
"""
Genotype the native marker genes of a sample from its reads.

The reads are aligned to the whole reference genome (so mapping quality means
what it should, and a read from a paralog or a repeat does not land in a marker
by default), then only the marker loci are looked at. Three kinds of evidence
are combined, because the three ways a marker is altered look different:

  pileup      substitutions and small indels, at a fully covered gene.
              ura3-1, trp1-1, ade2-1, can1-100 and lys2-801 are all of this
              kind, and coverage alone would call every one of them "present".
  depth       a run of bases no read covers, inside a locus that is otherwise
              covered: a deletion. his3D1, ura3D0, leu2D0, lys2D0.
  junctions   reads soft-clipped at the same position from both sides: foreign
              sequence joined to the locus. A Ty1 in ura3-52, or a cassette.
              Clips at the edges of a depth gap refine the deletion boundaries.

The resulting events are classified by marker_effects.classify, exactly as the
assembly alleles in the catalogue were, and the protein the sample encodes is
looked up there:

  * a protein form seen in active strains, natural isolates mostly, is a
    neutral polymorphism of an active gene, whatever the number of SNPs;
  * a form recorded as a documented laboratory allele is named (ura3-1, ...);
  * loss of function (nonsense, frameshift, start loss, deletion, insertion)
    is called from the sequence itself, known or not;
  * a missense form seen nowhere is reported as uncertain. It is never called
    active or inactive on a guess.

Heterozygous calls (a diploid, a mixed culture) give two alleles, both
reported. A one-copy deletion in a diploid halves the depth and is not
detected; neither is a cassette that left the locus length unchanged.
"""

import argparse
import bisect
import csv
import hashlib
import json
import re
import statistics
import sys

from marker_effects import (classify, describe_event, genome_event_to_locus,
                            read_fasta, revcomp)

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def load_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def short_hash(text):
    return hashlib.sha1(text.encode()).hexdigest()[:8]


def parse_pileup(path):
    """(chrom, pos 1-based) -> dict(depth, ref, alt{}, ins{}, dels{})."""
    out = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, bases = fields[0], int(fields[1]), fields[4]
            ref, alt, ins, dels, placeholder = 0, {}, {}, {}, 0
            i = 0
            while i < len(bases):
                c = bases[i]
                if c == "^":
                    i += 2
                    continue
                if c == "$":
                    i += 1
                    continue
                if c in "+-":
                    j = i + 1
                    while j < len(bases) and bases[j].isdigit():
                        j += 1
                    n = int(bases[i + 1:j])
                    seq = bases[j:j + n].upper()
                    target = ins if c == "+" else dels
                    target[seq] = target.get(seq, 0) + 1
                    i = j + n
                    continue
                if c in ".,":
                    ref += 1
                elif c in "*#":
                    placeholder += 1
                elif c in "<>":
                    pass
                else:
                    alt[c.upper()] = alt.get(c.upper(), 0) + 1
                i += 1
            out[(chrom, pos)] = {"depth": ref + sum(alt.values()) + placeholder,
                                 "ref": ref, "alt": alt, "ins": ins, "dels": dels}
    return out


def parse_depth(path):
    out = {}
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 3:
                out[(fields[0], int(fields[1]))] = int(fields[2])
    return out


def parse_reads(path, min_clip, min_mapq):
    """Soft-clip junctions and read alignments from the SAM text.

    Clips: chrom -> list of (boundary, side, clipped sequence). Side "R" means
    the read aligns to the right of the boundary and is clipped on its left;
    "L" the reverse. Boundaries are 0-based, between bases.

    Reads: chrom -> list of (start, end, indels), 0-based half open, with
    indels as (position, "I" or "D", length) at the reference boundary where
    they sit. Used to count, for an indel, the reads that could have shown it.
    """
    clips, reads, lengths = {}, {}, []
    with open(path) as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.split("\t")
            if len(fields) < 10:
                continue
            flag, chrom, pos, mapq, cigar, seq = (int(fields[1]), fields[2], int(fields[3]),
                                                  int(fields[4]), fields[5], fields[9])
            if flag & 0x904 or mapq < min_mapq or cigar == "*":
                continue
            lengths.append(len(seq))
            ops = [(int(n), op) for n, op in CIGAR.findall(cigar)]
            start = pos - 1
            ref, indels = start, []
            for n, op in ops:
                if op in "M=X":
                    ref += n
                elif op == "D" or op == "N":
                    indels.append((ref, "D", n))
                    ref += n
                elif op == "I":
                    indels.append((ref, "I", n))
            reads.setdefault(chrom, []).append((start, ref, indels))
            if ops[0][1] == "S" and ops[0][0] >= min_clip:
                clips.setdefault(chrom, []).append((start, "R", seq[:ops[0][0]]))
            if ops[-1][1] == "S" and ops[-1][0] >= min_clip:
                clips.setdefault(chrom, []).append((ref, "L", seq[-ops[-1][0]:]))
    # Sorted by start, with the longest reference span, so that the reads over a
    # position can be found by bisection instead of a scan.
    indexed = {}
    for chrom, rows in reads.items():
        rows.sort(key=lambda r: r[0])
        indexed[chrom] = ([r[0] for r in rows], rows, max(r[1] - r[0] for r in rows))
    return clips, indexed, lengths


def indel_support(reads, event, margin=5):
    """Reads carrying an indel, and reads spanning its place without it.

    A read ending close to an indel is soft-clipped by local alignment rather
    than given a gap, so it shows the reference there whatever the truth. The
    pileup fraction is diluted by those reads, badly so with short reads: the
    two +1 G insertions of leu2-3,112 came out at 40 to 70 per cent with 35 bp
    reads. Only reads with margin aligned bases on both sides are counted.
    """
    start, end, alt = event
    kind, length = ("I", len(alt)) if alt else ("D", end - start)
    carry, without = 0, 0
    starts, rows, longest = reads
    first = bisect.bisect_left(starts, start - margin - longest)
    last = bisect.bisect_right(starts, start - margin)
    for r_start, r_end, indels in rows[first:last]:
        if r_end < end + margin:
            continue
        if (start, kind, length) in indels:
            carry += 1
        elif not any(abs(p - start) <= 3 for p, _, _ in indels):
            without += 1
    return carry, without


def continues_reference(clip, genome_seq, g0, max_mismatch=0.5):
    """Is a clipped read end just the reference with a few differences.

    Local alignment also clips a read at a dense cluster of SNPs (three in five
    bases in Y55 and SK1 LYS2), and many reads clipped at one place from both
    sides look exactly like an insertion junction. A real junction is followed
    by foreign sequence, which differs from the locus at three bases in four by
    chance; a divergence clip is followed by the locus itself (Y55 LYS2: three
    of eight bases differ). Half separates the two.
    """
    boundary, side, seq = clip
    offset = boundary - g0
    if side == "R":
        ref = genome_seq[max(0, offset - len(seq)):offset]
        seq = seq[len(seq) - len(ref):]
    else:
        ref = genome_seq[offset:offset + len(seq)]
        seq = seq[:len(ref)]
    if len(ref) < 5:
        return False
    mismatches = sum(1 for a, b in zip(seq.upper(), ref) if a != b)
    return mismatches <= max_mismatch * len(ref)


def cluster_junctions(clips, window=5):
    """Group nearby boundaries; returns clusters with counts per side."""
    clusters = []
    for boundary, side, seq in sorted(clips):
        if clusters and boundary - clusters[-1]["last"] <= window:
            c = clusters[-1]
        else:
            c = {"positions": [], "L": 0, "R": 0, "seqs": [], "last": boundary}
            clusters.append(c)
        c["positions"].append(boundary)
        c[side] += 1
        c["seqs"].append(seq)
        c["last"] = boundary
    for c in clusters:
        c["position"] = int(statistics.median(c["positions"]))
        c["reads"] = c["L"] + c["R"]
    return clusters


def kmers(seq, k):
    return set(min(seq[i:i + k], revcomp(seq[i:i + k])) for i in range(len(seq) - k + 1))


def name_clips(seqs, tools, k=15):
    """Foreign element that the clipped sequences belong to, if any."""
    found = {}
    for seq in seqs:
        hits = set()
        for kmer in kmers(seq, k):
            for name, element in tools.items():
                if kmer in element:
                    hits.add(name)
        for name in hits:
            found[name] = found.get(name, 0) + 1
    if not found:
        return ""
    name, reads = max(found.items(), key=lambda kv: kv[1])
    return "%s (%d clipped reads)" % (name, reads) if reads >= 2 else ""


def call_small(pileup, reads, chrom, g0, g1, min_depth, hom, het, min_het_reads):
    """Substitutions and small indels in [g0, g1), split by zygosity.

    Also returns the positions where too few good reads exist to call anything:
    those are not the reference, they are unknown.
    """
    homs, hets, unreadable = [], [], []
    for pos in range(g0 + 1, g1 + 1):
        site = pileup.get((chrom, pos))
        if not site or site["depth"] < min_depth:
            unreadable.append(pos - 1)
            continue
        depth = float(site["depth"])
        candidates = []
        for base, n in site["alt"].items():
            if base in "ACGT":
                candidates.append((n, depth, (pos - 1, pos, base)))
        indels = [(pos, pos, seq) for seq in site["ins"]] + \
                 [(pos, pos + len(seq), "") for seq in site["dels"]]
        for event in indels:
            carry, without = indel_support(reads, event)
            if carry + without >= min_depth:
                candidates.append((carry, float(carry + without), event))
        for n, total, event in candidates:
            fraction = n / total
            if fraction >= hom:
                homs.append(event)
            elif fraction >= het and n >= min_het_reads:
                hets.append(event)
    return homs, hets, unreadable


def intervals_of(positions):
    out = []
    for p in sorted(positions):
        if out and p == out[-1][1]:
            out[-1] = (out[-1][0], p + 1)
        else:
            out.append((p, p + 1))
    return out


def in_mask(start, end, mask, fraction=0.5):
    if end <= start:
        return any(a <= start < b for a, b in mask)
    inside = sum(max(0, min(end, b) - max(start, a)) for a, b in mask)
    return inside >= fraction * (end - start)


def catalogue_events(row):
    """A catalogue row's events, with insertions of untranscribed sequence as None."""
    return [(a, b, None if (alt == "?" or alt.startswith("INS")) else alt)
            for a, b, alt in decode_events(row["events"])]


def hidden(event, missing, slack=10):
    """Could the reads have missed this event, given where they cannot see."""
    start = event[0]
    if event[1] == start:  # an insertion shows as a dip around its site
        return any(a - slack <= start <= b + slack for a, b in missing)
    return any(a <= start < b for a, b in missing)


def consistent_forms(rows, events, missing, orf_start, orf_end):
    """Catalogued alleles the reads agree with wherever they can be read.

    A divergent stretch (three SNPs in six bases in SK1 LYS2) stops short reads
    from mapping there at all, so the reads see the strain's other variants and
    nothing at that stretch. A Ty1 insertion leaves the locus sequence intact
    and shows, if at all, as a short dip around its target site. A catalogued
    allele is consistent when it carries every variant the reads show in the
    ORF, and every variant it has that the reads do not show sits where the
    reads cannot see.
    """
    def coding(evts):
        return set(e for e in evts if e[1] > orf_start - 1 and e[0] < orf_end)

    mine = coding(events)
    out = {}
    for r in rows:
        theirs = coding(catalogue_events(r))
        if not mine <= theirs:
            continue
        extra = theirs - mine
        if extra and all(hidden(e, missing) for e in extra):
            out.setdefault(r["protein_id"] or r["allele_id"], []).append(r)
    return out


def depth_gaps(depth, chrom, g0, g1, floor, min_len):
    """Runs of [start, end) genome positions with depth at or below floor."""
    gaps, run_start = [], None
    for g in range(g0, g1):
        low = depth.get((chrom, g + 1), 0) <= floor
        if low and run_start is None:
            run_start = g
        if not low and run_start is not None:
            if g - run_start >= min_len:
                gaps.append((run_start, g))
            run_start = None
    if run_start is not None and g1 - run_start >= min_len:
        gaps.append((run_start, g1))
    return gaps


def percentile(values, q):
    if not values:
        return 0.0
    values = sorted(values)
    return float(values[min(len(values) - 1, int(q * len(values)))])


def decode_events(text):
    out = []
    for item in filter(None, (text or "").split(";")):
        start, end, alt = item.split(":", 2)
        out.append((int(start), int(end), alt))
    return out


def overlap_fraction(a, b):
    inter = max(0, min(a[1], b[1]) - max(a[0], b[0]))
    return min(inter / float(a[1] - a[0]), inter / float(b[1] - b[0])) if inter else 0.0


def matches_structure(event, deletions, insertions, tolerance=20):
    """Does one sample event correspond to a catalogue deletion or insertion.

    A catalogue interval may be wider than the deletion itself, when the
    assembly alignment also lost some of the flank, so it only has to cover
    half of its own length with the sample deletion; the sample deletion must
    fall at least 80 per cent inside it.
    """
    start, end, _ = event
    if end > start:
        for a, b in deletions:
            inter = max(0, min(end, b) - max(start, a))
            if inter >= 0.8 * (end - start) and inter >= 0.5 * (b - a):
                return True
        return False
    return any(abs(start - a) <= tolerance for a in insertions)


def snap_deletions(events, rows, tolerance=15):
    """Replace depth-estimated deletions by an exact catalogue deletion.

    Depth only places a breakpoint to within a few bases, since reads overlap
    the edge by a variable amount and local alignment trims mismatching ends.
    When a catalogued deletion of the same gene has both ends within that
    error, it is taken as the exact event: the reads are consistent with it,
    and its sequence, hence its protein, is known exactly.
    """
    out, notes = [], []
    catalogued = {}
    for r in rows:
        for a, b, alt in decode_events(r["events"]):
            if b - a >= 30 and not alt.startswith("INS") and len(alt) < 30:
                catalogued.setdefault((a, b, alt), set()).add(r["strain"])
    for event in events:
        start, end, alt = event
        best = None
        if end > start and alt == "":
            for (a, b, calt), strains in catalogued.items():
                shift = max(abs(a - start), abs(b - end))
                if shift <= tolerance and (best is None or shift < best[0]):
                    best = (shift, (a, b, calt), strains)
        if best:
            out.append(best[1])
            notes.append("deletion of %d bp from depth matches the catalogued deletion of "
                         "%d bp in %s within %d bp at both ends; its exact breakpoints are used"
                         % (end - start, best[1][1] - best[1][0],
                            ", ".join(sorted(best[2])[:4]), best[0]))
        else:
            out.append(event)
    return out, notes


def protein_id(gene, protein):
    return "%s-p%s" % (gene, short_hash(protein)) if protein else ""


def known_forms(catalogue, gene, effect, protein, locus_events, orf):
    """What the catalogue knows about this allele.

    Returns carriers (assembly names), named documented alleles, and the
    catalogue's verdict on the protein form.
    """
    rows = [r for r in catalogue if r["gene"] == gene]
    pid = protein_id(gene, protein)
    carriers = [r for r in rows if pid and r["protein_id"] == pid]
    # Only the events that touch the ORF decide the allele; a structural
    # variant further upstream is reported but not required to match.
    big = [e for e in locus_events
           if (e[2] is None or (e[1] - e[0]) >= 30 or len(e[2] or "") >= 30)
           and e[1] >= orf[0] and e[0] <= orf[1]]
    if not carriers and big:
        # A deletion or an insertion has no protein worth comparing, or the
        # protein depends on breakpoints the reads only give to a few bases.
        # The structure is compared instead.
        for r in rows:
            if r["activity"] != "inactive":
                continue
            theirs = [(a, b) for a, b, alt in decode_events(r["events"]) if b - a >= 30]
            theirs += [tuple(map(int, x.split("-"))) for x in
                       filter(None, (r.get("uncovered") or "").split(";"))]
            inserts = [a for a, b, alt in decode_events(r["events"])
                       if alt == "?" or alt.startswith("INS") or len(alt) >= 100]
            if all(matches_structure(e, theirs, inserts) for e in big):
                carriers.append(r)
    named = sorted(set("%s (%s)" % (r["documented_allele"], r["strain"]) for r in carriers
                       if r["documented_allele"] and r["documented_allele"] != gene.upper()))
    wild = sorted(set(r["strain"] for r in carriers
                      if r["documented_allele"] == gene.upper()))
    return carriers, named, wild


def judge(effect_result, carriers, named, wild):
    """Final activity, combining the sequence effect with the catalogue."""
    activity = effect_result["activity"]
    effect = effect_result["effect"]
    if activity == "inactive":
        return "inactive", "loss of function from the sequence itself (%s)" % effect
    if activity == "active":
        return "active", ("reference protein" if effect == "reference"
                          else "same protein as the reference")
    # Uncertain from sequence: the catalogue decides, or nothing does.
    if named and not wild:
        return "inactive", "protein form of the documented allele %s" % ", ".join(named)
    if wild and not named:
        return "active", "protein form carried by documented wild-type %s" % ", ".join(wild)
    if named and wild:
        return "uncertain", ("protein form shared by documented mutant %s and wild-type %s"
                             % (", ".join(named), ", ".join(wild)))
    if carriers:
        active_like = [r for r in carriers if r["activity"] != "inactive"]
        return ("active",
                "natural form: same protein in %d assemblies (%s%s), none documented as a "
                "mutant; isolates are overwhelmingly prototrophic, so presumed active"
                % (len(active_like), ", ".join(sorted(set(r["strain"] for r in active_like))[:6]),
                   ", ..." if len(active_like) > 6 else ""))
    return "uncertain", "%s not seen in any catalogued strain" % effect


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--loci", required=True)
    parser.add_argument("--loci-fasta", required=True)
    parser.add_argument("--alleles", required=True, help="marker_alleles.tsv catalogue")
    parser.add_argument("--pileup", required=True, help="samtools mpileup over the loci")
    parser.add_argument("--depth", required=True, help="samtools depth -a, all mapping qualities")
    parser.add_argument("--sam", required=True, help="reads over the loci, SAM text")
    parser.add_argument("--documented", help="documented_genotypes.tsv")
    parser.add_argument("--tools", help="tools.fasta, to name what clipped reads belong to")
    parser.add_argument("--tools-manifest", help="tools.tsv, only foreign elements are used")
    parser.add_argument("--min-depth", type=int, default=4)
    parser.add_argument("--hom-fraction", type=float, default=0.8)
    parser.add_argument("--het-fraction", type=float, default=0.25)
    parser.add_argument("--min-gap", type=int, default=30,
                        help="shortest uncovered run called a deletion (default 30)")
    parser.add_argument("--min-clip", type=int, default=8)
    parser.add_argument("--min-het-reads", type=int, default=3,
                        help="reads needed for a heterozygous call (default 3). Two is "
                             "what sequencing error alone gives at a few positions per locus")
    parser.add_argument("--repeats", help="marker_repeats.tsv, flank intervals repeated in the genome")
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--out", required=True, help="JSON output")
    parser.add_argument("--tsv", help="one line per gene")
    args = parser.parse_args()

    loci = load_tsv(args.loci)
    seqs = dict((n.split()[0], s) for n, s in read_fasta(args.loci_fasta))
    catalogue = load_tsv(args.alleles)
    pileup = parse_pileup(args.pileup)
    depth = parse_depth(args.depth)
    clips, reads, read_lengths = parse_reads(args.sam, args.min_clip, args.min_mapq)
    read_length = statistics.median(read_lengths) if read_lengths else 0
    repeats = {}
    if args.repeats:
        for row in load_tsv(args.repeats):
            repeats.setdefault(row["locus_id"], []).append((int(row["start"]), int(row["end"])))
    documented = {}
    if args.documented:
        genes = [l["gene"] for l in loci]
        for row in load_tsv(args.documented):
            tokens = row["documented_genotype"].split()
            alleles = {}
            for gene in genes:
                mutant = [t for t in tokens if t.lower().startswith(gene.lower())
                          and not t.startswith(gene.upper())]
                # A diploid heterozygous for a marker (trp1D63/TRP1 in FY1679)
                # carries one working copy, and no active/inactive call
                # compares with it.
                if mutant and gene.upper() in mutant[0].split("/"):
                    alleles[gene] = None
                else:
                    alleles[gene] = mutant[0].split("/")[0] if mutant else gene.upper()
            documented[row["strain"]] = alleles

    tools = {}
    if args.tools:
        foreign = None
        if args.tools_manifest:
            foreign = set(r["element_id"] for r in load_tsv(args.tools_manifest)
                          if r.get("origin") == "foreign")
        for name, seq in read_fasta(args.tools):
            if foreign is None or name.split()[0] in foreign:
                tools[name.split()[0].split("__")[0]] = kmers(seq, 15)

    for locus in loci:
        for key in ("genome_start", "genome_end", "orf_start", "orf_end"):
            locus[key] = int(locus[key])
        g0, g1 = locus["genome_start"], locus["genome_end"]
        values = [depth.get((locus["chrom"], g + 1), 0) for g in range(g0, g1)]
        locus["_p90"] = percentile(values, 0.9)
    baseline = statistics.median([l["_p90"] for l in loci]) if loci else 0.0

    results = []
    for locus in loci:
        gene, chrom = locus["gene"], locus["chrom"]
        g0, g1 = locus["genome_start"], locus["genome_end"]
        locus_seq = seqs[locus["locus_id"]]
        orf_start, orf_end = locus["orf_start"], locus["orf_end"]
        local = locus["_p90"] if locus["_p90"] >= 0.3 * baseline else baseline
        record = {
            "gene": gene, "systematic_id": locus["systematic_id"],
            "loss_of_function_phenotype": locus["loss_of_function_phenotype"],
            "locus": "%s:%d-%d (%s)" % (chrom, g0 + 1, g1, locus["strand"]),
            "depth_p90": locus["_p90"], "baseline_depth": baseline,
        }
        if baseline < args.min_depth:
            record.update({"status": "no_data", "explanation":
                           "the marker loci are covered at %.0fx, below the %dx needed"
                           % (baseline, args.min_depth), "alleles": []})
            results.append(record)
            continue
        if locus["_p90"] < 0.1 * baseline:
            orf_values = [depth.get((chrom, g + 1), 0) for g in range(g0, g1)]
            record["orf_median_depth"] = statistics.median(orf_values)

        notes = []
        # Deletions, from depth, with boundaries snapped to clip clusters.
        floor = max(1.0, 0.1 * local)
        gaps = depth_gaps(depth, chrom, g0, g1, floor, args.min_gap)
        # A flank repeat (a Ty LTR or delta) gets its reads shared out with every
        # other copy in the genome, so a dip there is not a deletion.
        mask = repeats.get(locus["locus_id"], [])
        genome_mask = [(g0 + a, g0 + b) if locus["strand"] == "+" else (g1 - b, g1 - a)
                       for a, b in mask]
        kept_gaps = []
        for s_, e_ in gaps:
            if in_mask(s_, e_, genome_mask):
                notes.append("low depth at %s:%d-%d lies in a repeated flank, not called"
                             % (chrom, s_ + 1, e_))
            else:
                kept_gaps.append((s_, e_))
        gaps = kept_gaps
        genome_seq = locus_seq if locus["strand"] == "+" else revcomp(locus_seq)
        junction_clips = [c for c in clips.get(chrom, [])
                          if g0 <= c[0] <= g1 and not continues_reference(c, genome_seq, g0)]
        clusters = [c for c in cluster_junctions(junction_clips)
                    if c["reads"] >= max(3, 0.15 * local)]
        used = set()
        structural = []
        large = []
        for s, e in gaps:
            left = [i for i, c in enumerate(clusters) if c["L"] and abs(c["position"] - s) <= 20]
            right = [i for i, c in enumerate(clusters) if c["R"] and abs(c["position"] - e) <= 20]
            if left:
                s = clusters[left[0]]["position"]
                used.add(left[0])
            if right:
                e = clusters[right[0]]["position"]
                used.add(right[0])
            names = [name_clips(clusters[i]["seqs"], tools) for i in left + right]
            names = sorted(set(n for n in names if n))
            large.append((s, e, ""))
            reaches = []
            if s <= g0 + 20:
                reaches.append("upstream")
            if e >= g1 - 20:
                reaches.append("downstream")
            structural.append("no read covers %s:%d-%d (%d bp)%s%s" % (
                chrom, s + 1, e, e - s,
                ", the gap runs past the %s end of the locus" % " and ".join(reaches)
                if reaches else "",
                "; clipped reads at the edges belong to %s" % ", ".join(names) if names else ""))
        for i, c in enumerate(clusters):
            if i in used:
                continue
            name = name_clips(c["seqs"], tools)
            if c["L"] and c["R"]:
                large.append((c["position"], c["position"], None))
                structural.append("foreign sequence joined at %s:%d, %d reads clipped from "
                                  "both sides%s" % (chrom, c["position"], c["reads"],
                                                    ", it belongs to " + name if name else ""))
            else:
                notes.append("one-sided junction at %s:%d (%d reads)%s, not called" % (
                    chrom, c["position"], c["reads"], ", " + name if name else ""))

        homs, hets, unreadable = call_small(pileup, reads.get(chrom, ([], [], 0)), chrom, g0, g1,
                                            args.min_depth, args.hom_fraction,
                                            args.het_fraction, args.min_het_reads)
        # Small calls inside a deletion are noise from the few reads left there.
        def outside(event):
            return not any(s <= event[0] < e for s, e, _ in large if e > s)
        homs = [e for e in homs if outside(e) and not in_mask(e[0], e[1], genome_mask)]
        hets = [e for e in hets if outside(e) and not in_mask(e[0], e[1], genome_mask)]
        unreadable = [p for p in unreadable
                      if not any(s <= p < e for s, e, _ in large if e > s)]

        def to_locus(events):
            return sorted(genome_event_to_locus(e, locus) for e in events)

        gene_rows = [r for r in catalogue if r["gene"] == gene]
        large_locus, snap_notes = snap_deletions(to_locus(large), gene_rows)
        notes.extend(snap_notes)
        haplotypes = [("hom", sorted(to_locus(homs) + large_locus))]
        if hets:
            haplotypes.append(("het", sorted(to_locus(homs + hets) + large_locus)))

        missing = [iv for iv in intervals_of(
            [genome_event_to_locus((p, p + 1, ""), locus)[0] for p in unreadable])
            if iv[1] > orf_start and iv[0] < orf_end]
        missing_bp = sum(min(e, orf_end) - max(s, orf_start) for s, e in missing)
        alleles = []
        for label, events in haplotypes:
            effect = classify(locus_seq, orf_start, orf_end, events)
            carriers, named, wild = known_forms(catalogue, gene, effect["effect"],
                                                effect.get("protein"), events,
                                                (orf_start, orf_end))
            activity, reason = judge(effect, carriers, named, wild)
            if missing_bp and activity != "inactive":
                # The reads cannot see part of the ORF. What matters is whether
                # anything known could be hiding there. At 12x a few bases fall
                # below the depth floor by chance alone in every gene; if no
                # catalogued allele differs from the reference at those bases,
                # they are read as the reference and the call stands.
                forms = consistent_forms(gene_rows, events, missing, orf_start, orf_end)
                outside_orf = [e for e in events if e[1] <= orf_start or e[0] >= orf_end]
                verdicts = []
                for key, rows in sorted(forms.items()):
                    theirs = [e for e in catalogue_events(rows[0])
                              if e[1] > orf_start - 1 and e[0] < orf_end]
                    form_effect = classify(locus_seq, orf_start, orf_end,
                                           sorted(set(theirs + outside_orf),
                                                  key=lambda e: (e[0], e[1])))
                    c2, n2, w2 = known_forms(catalogue, gene, form_effect["effect"],
                                             form_effect.get("protein"), theirs,
                                             (orf_start, orf_end))
                    verdicts.append((judge(form_effect, c2, n2, w2), form_effect, c2, n2, w2, rows))
                activities = set(v[0][0] for v in verdicts)
                if missing_bp > 0.1 * (orf_end - orf_start):
                    activity = "uncertain"
                    reason = ("%d of %d bp of the ORF cannot be read from these reads"
                              % (missing_bp, orf_end - orf_start))
                elif not verdicts:
                    notes.append("%d bp of the ORF are below the depth floor; no catalogued "
                                 "allele differs from the reference there, so they are read "
                                 "as reference" % missing_bp)
                elif len(activities) == 1 and (activity == "uncertain"
                                               or activities == set([activity])):
                    # Every allele that fits agrees, and either agrees with what
                    # is read or resolves a form the reads alone left unknown.
                    (activity, why), effect, carriers, named, wild, rows = verdicts[0]
                    reason = ("%d bp of the ORF cannot be read from these reads; everything "
                              "that can be read matches the catalogued allele of %s, which "
                              "differs from what is read only there. %s"
                              % (missing_bp,
                                 ", ".join(sorted(set(r["strain"] for r in rows))[:4]), why))
                else:
                    def carriers_of(rows):
                        named_first = sorted(set(r["documented_allele"] for r in rows
                                                 if r["documented_allele"]))
                        return named_first + sorted(set(r["strain"] for r in rows))
                    others = sorted(set("%s (%s)" % (v[0][0], ", ".join(carriers_of(v[5])[:3]))
                                        for v in verdicts))
                    activity = "uncertain"
                    reason = ("%d bp of the ORF cannot be read from these reads, and catalogued "
                              "alleles that differ there would change the call: %s"
                              % (missing_bp, "; ".join(others)))
                if effect["effect"] == "partial":
                    effect = classify(locus_seq, orf_start, orf_end, events)
            exact = sorted(set(r["strain"] for r in carriers
                               if r["allele_id"] == "%s-%s" % (gene, short_hash(
                                   ";".join("%d:%d:%s" % (s, e, a if a is not None else "?")
                                            for s, e, a in events)) if events else "ref")))
            alleles.append({
                "copy": "all reads" if label == "hom" and not hets else
                        ("first allele (homozygous calls only)" if label == "hom"
                         else "second allele (with heterozygous calls)"),
                "activity": activity, "reason": reason,
                "effect": effect["effect"], "protein_change": effect.get("protein_change", ""),
                "coding_events": effect.get("coding", []),
                "upstream_events": effect.get("upstream", []),
                "downstream_events": effect.get("downstream", []),
                "protein_form": protein_id(gene, effect.get("protein")),
                "known_as": named,
                "documented_wild_type_carriers": wild,
                "carriers": len(set(r["assembly"] for r in carriers)),
                "carrier_examples": sorted(set(r["strain"] for r in carriers))[:10],
                "identical_nucleotide_allele_in": exact[:10],
            })
        statuses = set(a["activity"] for a in alleles)
        status = alleles[0]["activity"] if len(statuses) == 1 else "mixed"
        record.update({
            "status": status,
            "alleles": alleles,
            "structural": structural,
            "notes": notes,
            "heterozygous_calls": [describe_event(e, locus_seq, orf_start, orf_end)
                                   for e in to_locus(hets)],
        })
        results.append(record)

    # Which documented strains this marker genotype is consistent with. Only
    # active against inactive is compared: documented genotypes name alleles,
    # not protein forms, and "uncertain" or "mixed" calls are left out.
    matches = []
    for strain, alleles in documented.items():
        agree, compared, differ = 0, 0, []
        for r in results:
            if r["status"] not in ("active", "inactive"):
                continue
            allele = alleles.get(r["gene"], r["gene"].upper())
            if allele is None:
                continue
            expected = "active" if allele == r["gene"].upper() else "inactive"
            compared += 1
            if expected == r["status"]:
                agree += 1
            else:
                differ.append("%s %s, documented %s" % (r["gene"], r["status"], allele))
        if compared:
            matches.append({"strain": strain, "agree": agree, "compared": compared,
                            "differences": differ})
    matches.sort(key=lambda m: (-m["agree"], m["strain"]))

    def genotype_token(record):
        if record["status"] == "no_data":
            return "%s(?)" % record["gene"]
        allele = record["alleles"][0]
        if record["status"] == "active":
            return record["gene"].upper()
        if record["status"] == "inactive":
            # Several documented alleles can encode the same protein (ade2-1
            # and ade2-101 are the same ochre codon); the reads cannot choose.
            name = "/".join(sorted(set(k.split(" ")[0] for k in allele["known_as"]))) \
                if allele["known_as"] else \
                "%s[%s]" % (record["gene"].lower(), allele["effect"].replace("_", " "))
            return name
        return "%s(%s)" % (record["gene"].lower(), record["status"])

    limits = []
    if read_length and read_length < 50:
        limits.append(
            "reads are %d bp: an insertion that leaves the locus sequence intact, a Ty1 as in "
            "ura3-52 or a cassette inserted without deletion, is seen only through reads "
            "crossing its junction, and too few reads this short cross it to call it. Such a "
            "gene can be reported active when it is not." % read_length)
    limits.append("a deletion of one copy in a diploid halves the depth and is not called")
    report = {
        "baseline_depth": baseline,
        "read_length": read_length,
        "limits": limits,
        "genotype": " ".join(genotype_token(r) for r in results),
        "markers": results,
        "documented_genotype_matches": matches,
        "catalogue_assemblies": len(set(r["assembly"] for r in catalogue)),
    }
    with open(args.out, "w") as handle:
        json.dump(report, handle, indent=1)
    if args.tsv:
        with open(args.tsv, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene", "status", "effect", "protein_change", "known_as",
                             "reason", "structural"])
            for r in results:
                a = r["alleles"][0] if r["alleles"] else {}
                writer.writerow([r["gene"], r["status"], a.get("effect", ""),
                                 a.get("protein_change", ""), ";".join(a.get("known_as", [])),
                                 a.get("reason", r.get("explanation", "")),
                                 " | ".join(r.get("structural", []))])
    sys.stderr.write("genotype: %s\n" % report["genotype"])
    for r in results:
        a = r["alleles"][0] if r["alleles"] else {}
        sys.stderr.write("  %-6s %-9s %-14s %s %s\n" % (
            r["gene"], r["status"], a.get("effect", ""), a.get("protein_change", ""),
            a.get("reason", r.get("explanation", ""))))
    for m in matches[:4]:
        sys.stderr.write("documented %-10s %d/%d %s\n" % (
            m["strain"], m["agree"], m["compared"], "; ".join(m["differences"])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
