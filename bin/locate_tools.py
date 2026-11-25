#!/usr/bin/env python3
"""
Where the tool elements present in a sample sit in its genome.

The reads are aligned to the whole reference genome in local mode by the
MarkerAlign process, so a read crossing the junction between the genome and an
inserted construct aligns on its genomic part and is soft-clipped where the
construct begins. The evidence read here:

  clips    reads whose genomic part is well mapped (mapping quality 20 or more)
           with 15 clipped bases or more, the clipped bases being found in an
           element (15-mers). The edge of the clip is the junction, to the base
           (give or take the few bases the genome and the element may share
           there). It needs reads long enough to hold 15 clipped bases beside a
           mappable part: about 50 bases or more, or the pieces of long
           sequences.
  mates    reads found in an element whose mate, from the same fragment, is well
           mapped on the genome: the junction lies within a fragment length of
           the mate, on the side the mate points to. Paired reads only; it is
           what places a construct when the reads are too short to be clipped.
  pieces   for long sequences cut into pieces (an assembly, long reads), the
           pieces of the same sequence near a junction piece: the rest of the
           construct, which tells apart cassettes that share their ends (kanMX,
           natMX6, hphMX6 and HIS3MX6 all start with the same TEF promoter and
           end with the same TEF terminator).
  chunks   50 base pieces of each element aligned to the reference: where the
           reference itself carries the element, i.e., the native locus of a
           yeast promoter, terminator or gene used in cassettes. A read anchored
           there cannot tell the native copy from a copy in a construct, so these
           anchors do not place an element; they link the native element to the
           foreign one joined to it (the ADH1 terminator between a GFP tag and
           its marker).

Junctions on one chromosome within --max-site-span bases make one insertion
site: the left junction (genome, then construct) and the right junction
(construct, then genome) of one insertion. Each piece of evidence of a site is
matched to the elements; an element whose matches are a strict subset of
another's is not at the site (they are sequence it shares with that one),
unless it is contained in it (a promoter inside its cassette). Elements left
with the same matches are alternatives: the reads do not tell which one is
there, short reads crossing the TEF promoter of a cassette for instance, and
those not absent from the sample are reported as such. A site lists the
elements left, the elements contained in them, and the native elements joined
to one of them. An element present with no site is either native (its native
locus is given) or not joined to the genome in the reads (a plasmid, an
element inside a longer construct whose ends are not in the database, or reads
too short and unpaired).

Only the elements present (score >= 99, green in the report, the score of
score_strains.py) are written.
"""

import argparse
import csv
import gzip
import json
import re
import statistics
import sys

from score_strains import score_tools

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
K = 15


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def load_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path):
    """{first word of the header: sequence, upper case}."""
    out, name, chunks = {}, None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(chunks).upper()
                name, chunks = line[1:].split()[0] if line[1:].split() else "seq", []
            elif name is not None:
                chunks.append(line)
    if name is not None:
        out[name] = "".join(chunks).upper()
    return out


def chromosome_names(path):
    """Display name of each reference sequence, from the FASTA headers:
    [chromosome=IV] -> chrIV, [location=mitochondrion] -> chrM, else the first
    word. Only the header lines are read."""
    names = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            fields = line[1:].split()
            if not fields:
                continue
            chrom = re.search(r"\[chromosome=([^\]\s]+)\]", line)
            if chrom:
                names[fields[0]] = "chr" + chrom.group(1)
            elif "[location=mitochondrion]" in line:
                names[fields[0]] = "chrM"
            else:
                names[fields[0]] = fields[0]
    return names


def reference_span(cigar):
    return sum(int(n) for n, op in CIGAR.findall(cigar) if op in "MD=XN")


def strip_mate(name):
    return name[:-2] if name.endswith(("/1", "/2")) else name


def piece_origin(name):
    """(sequence, start) of a piece of a long sequence, named <sequence>_<start>
    by prepare_reads.py, or None."""
    base, _, start = name.rpartition("_")
    return (base, int(start)) if base and start.isdigit() else None


def short(eid):
    return eid.split("__")[0]


def majority(values):
    counts = {}
    for v in values:
        counts[v] = counts.get(v, 0) + 1
    return max(sorted(counts), key=lambda v: counts[v]) if counts else ""


# --------------------------------------------------------------------------- #
# element sequences
# --------------------------------------------------------------------------- #

def kmer_index(elements):
    """15-mer -> [(element_id, strand, position on that strand)], from both
    strands of every element. Low-complexity k-mers (fewer than three different
    bases: poly-A, AT repeats) are left out: reads of the genome carry them
    everywhere."""
    index = {}
    for eid, seq in elements.items():
        for strand, s in (("+", seq), ("-", revcomp(seq))):
            for i in range(len(s) - K + 1):
                kmer = s[i:i + K]
                if "N" in kmer or len(set(kmer)) < 3:
                    continue
                index.setdefault(kmer, []).append((eid, strand, i))
    return index


def match(seq, index, min_hits=3):
    """Elements a sequence comes from.

    Returns {element_id: (strand, diagonal, hits)}: strand of the element the
    sequence reads, diagonal = position on that strand of the first base of the
    sequence (negative when the sequence starts before the element, a linker
    for instance), hits = 15-mers found. An element is kept with at least
    min_hits k-mers (all of them when the sequence holds fewer) and 80 % of the
    hits of the best element: a sequence of the TEF promoter alone matches every
    cassette built on it, a longer one only the cassette it comes from.
    """
    votes = {}
    n = len(seq) - K + 1
    for i in range(max(n, 0)):
        for eid, strand, pos in index.get(seq[i:i + K], ()):
            votes.setdefault((eid, strand), []).append(pos - i)
    if not votes:
        return {}
    best = {}
    for (eid, strand), diagonals in votes.items():
        if eid not in best or len(diagonals) > best[eid][2]:
            best[eid] = (strand, int(statistics.median(diagonals)), len(diagonals))
    need = min(min_hits, n)
    top = max(v[2] for v in best.values())
    return {eid: v for eid, v in best.items() if v[2] >= need and v[2] >= 0.8 * top}


def containment(elements, fraction=0.9):
    """{element: set of the longer elements holding at least 90 % of its
    15-mers}: a promoter inside a cassette, a resistance gene inside its
    cassette."""
    sets = {}
    for eid, seq in elements.items():
        s = set()
        for i in range(len(seq) - K + 1):
            kmer = seq[i:i + K]
            s.add(min(kmer, revcomp(kmer)))
        sets[eid] = s
    out = {}
    for a, ka in sets.items():
        out[a] = set(b for b, kb in sets.items()
                     if ka and b != a and len(kb) > len(ka) and len(ka & kb) >= fraction * len(ka))
    return out


# --------------------------------------------------------------------------- #
# native loci, from the element chunks aligned to the reference
# --------------------------------------------------------------------------- #

def native_loci(chunks_path, elements, min_score=80):
    """Native locus of each element carried by the reference, and the stretches
    of the reference that any element matches (the mask).

    Chunks are 50 bases long, one every 10 bases (tool_chunks.tsv: name
    element:start, flag, chromosome, position, MAPQ, CIGAR, alignment score). A
    chunk counts with a local alignment score of 80 or more (40 bases or more
    matching, a few mismatches). The chunks of one copy line up along the
    reference; the copy holding the most chunks is the native locus when it
    holds half of the element or more.
    """
    total = {eid: (len(seq) - 50) // 10 + 1 + (1 if (len(seq) - 50) % 10 else 0)
             for eid, seq in elements.items() if len(seq) >= 50}
    hits, mask = {}, {}
    with open(chunks_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7 or fields[2] == "*":
                continue
            eid, _, offset = fields[0].rpartition(":")
            try:
                score = int(fields[6])
                offset = int(offset) - 1
            except ValueError:
                continue
            if score < min_score:
                continue
            flag, chrom, pos, mapq = int(fields[1]), fields[2], int(fields[3]) - 1, int(fields[4])
            end = pos + reference_span(fields[5])
            mask.setdefault(chrom, []).append((pos, end))
            strand = "-" if flag & 16 else "+"
            # Reference position of the element's first base, the same for every
            # chunk of one copy (give or take an indel).
            anchor = pos - offset if strand == "+" else end + offset
            hits.setdefault(eid, []).append((chrom, strand, anchor, pos, end, mapq))
    loci = {}
    for eid, rows in hits.items():
        rows.sort(key=lambda r: (r[0], r[1], r[2]))
        best, group = [], []
        for row in rows:
            if group and (row[0], row[1]) == (group[-1][0], group[-1][1]) \
                    and row[2] - group[-1][2] <= 200:
                group.append(row)
            else:
                group = [row]
            if len(group) > len(best):
                best = list(group)
        if eid in total and len(best) >= 0.5 * total[eid]:
            loci[eid] = {
                "chrom": best[0][0], "strand": best[0][1],
                "start": min(r[3] for r in best), "end": max(r[4] for r in best),
                "fraction": round(min(1.0, len(best) / float(total[eid])), 2),
                "repeated": statistics.median(r[5] for r in best) < 10,
            }
    merged = {}
    for chrom, spans in mask.items():
        spans.sort()
        out = [list(spans[0])]
        for s, e in spans[1:]:
            if s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        merged[chrom] = out
    return loci, merged


def masked_bases(mask, chrom, start, end):
    return sum(max(0, min(end, e) - max(start, s)) for s, e in mask.get(chrom, ()))


def near_native(eid, chrom, pos, loci, distance):
    locus = loci.get(eid)
    return bool(locus) and locus["chrom"] == chrom and \
        locus["start"] - distance <= pos <= locus["end"] + distance


def natives_at(loci, chrom, start, end):
    return [eid for eid, locus in loci.items()
            if locus["chrom"] == chrom and locus["start"] < end and start < locus["end"]]


# --------------------------------------------------------------------------- #
# evidence
# --------------------------------------------------------------------------- #

def joined_end(strand, coordinate, length, slack=50):
    """Which end of the element a junction joins, from the coordinate, on the
    strand read, of the element base next to the genome."""
    if coordinate <= slack:
        return "start" if strand == "+" else "end"
    if coordinate >= length - slack:
        return "end" if strand == "+" else "start"
    return "internal"


def read_clips(path, index, elements, mask, loci, near):
    """Junction reads, and the links between a native element and the element
    joined to it.

    tool_clips.tsv.gz: chromosome, boundary (0-based, between two bases), side,
    start and end of the aligned part, clipped bases, read name. Side L: the
    read aligns on the left of the boundary and the construct follows; R: the
    construct precedes. A read whose aligned part lies within an element-like
    stretch of the reference (a native locus) is set apart as a link between
    that native element and the clipped one.
    """
    clips, links = [], []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            chrom, boundary, side = fields[0], int(fields[1]), fields[2]
            a0, a1, seq = int(fields[3]), int(fields[4]), fields[5].upper()
            name = fields[6] if len(fields) > 6 else ""
            hit = {}
            for eid, (strand, diagonal, _) in match(seq, index).items():
                if near_native(eid, chrom, boundary, loci, near):
                    # The genome beside the element's own native locus: a local
                    # change of that locus, not an insertion.
                    continue
                coordinate = diagonal if side == "L" else diagonal + len(seq)
                hit[eid] = (strand, joined_end(strand, coordinate, len(elements[eid])))
            if not hit:
                continue
            if (a1 - a0) - masked_bases(mask, chrom, a0, a1) < 20:
                for native in natives_at(loci, chrom, a0, a1):
                    links.append((native, frozenset(eid for eid in hit if eid != native)))
                continue
            clips.append({"chrom": chrom, "boundary": boundary, "side": side, "hit": hit,
                          "name": name})
    return clips, links


def read_mates(reads_path, mates_path, index, mask, loci, near):
    """Read pairs with one read in an element and its mate on the genome, and
    the links between a native element and the element joined to it.

    tool_reads.tsv: name, flag, element, position, CIGAR, sequence of the reads
    aligned end to end to the elements. tool_mates.tsv: the same columns for the
    genome alignments (mapping quality 20 or more) of the reads of the same
    names, mates included. A genome read whose own sequence is found in an
    element is the element read itself, or another read of the construct: it
    places nothing, but in a native locus it links that native element to the
    element its mate is in.
    """
    element_reads = {}
    with open(reads_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            seq = fields[5].upper()
            original = revcomp(seq) if int(fields[1]) & 16 else seq
            element_reads.setdefault(strip_mate(fields[0]), []).append(original)
    pairs, links = [], []
    with open(mates_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            name = strip_mate(fields[0])
            if name not in element_reads:
                continue
            flag, chrom, pos = int(fields[1]), fields[2], int(fields[3]) - 1
            seq = fields[5].upper()
            end = pos + reference_span(fields[4])
            original = revcomp(seq) if flag & 16 else seq
            forward = not flag & 16
            others = [r for r in element_reads[name] if r != original]
            if not others:
                continue
            own = match(seq, index, min_hits=1)
            in_mask = masked_bases(mask, chrom, pos, end) >= 0.5 * (end - pos)
            if (own and max(v[2] for v in own.values()) >= 0.3 * max(1, len(seq) - K + 1)) \
                    or in_mask:
                if in_mask:
                    for n in natives_at(loci, chrom, pos, end):
                        for other in others:
                            links.append((n, frozenset(
                                eid for eid in match(other, index)
                                if eid != n and not near_native(eid, chrom, pos, loci, near))))
                continue
            for other in others:
                # The mates point at each other: the junction is on the right of
                # a forward mate (side L), on the left of a reverse one (side R),
                # and the element read, seen on the strand of the genome, is the
                # reverse complement of its sequence for a forward mate.
                framed = revcomp(other) if forward else other
                hit = {eid: (v[0], "") for eid, v in match(framed, index).items()
                       if not near_native(eid, chrom, pos, loci, near)}
                if hit:
                    pairs.append({"chrom": chrom, "start": pos, "end": end,
                                  "side": "L" if forward else "R", "hit": hit})
    return pairs, links


# --------------------------------------------------------------------------- #
# junctions and sites
# --------------------------------------------------------------------------- #

def clip_junctions(clips, min_reads, window=3):
    """Clips clustered by chromosome, side and boundary (within 3 bases), kept
    with min_reads clips or more. The junction is at the boundary most reads
    share."""
    out, group = [], []
    clips.sort(key=lambda c: (c["chrom"], c["side"], c["boundary"]))

    def close(group):
        if len(group) >= min_reads:
            out.append({"chrom": group[0]["chrom"], "side": group[0]["side"],
                        "position": majority([c["boundary"] for c in group]),
                        "exact": True, "clips": group, "pairs": []})

    for c in clips:
        if group and (c["chrom"], c["side"]) == (group[-1]["chrom"], group[-1]["side"]) \
                and c["boundary"] - group[-1]["boundary"] <= window:
            group.append(c)
        else:
            if group:
                close(group)
            group = [c]
    if group:
        close(group)
    return out


def pair_junctions(pairs, junctions, min_pairs, window):
    """Read pairs added to the clipped junction they point at (same side,
    within a fragment length), or else grouped into a junction of their own,
    kept with min_pairs pairs or more and placed at the mate nearest to it: the
    end of the last forward mate (side L), the start of the first reverse one
    (side R)."""
    loose = []
    for p in pairs:
        target = None
        for j in junctions:
            if j["chrom"] == p["chrom"] and j["side"] == p["side"] and j["exact"] and (
                    (p["side"] == "L" and 0 <= j["position"] - p["start"] <= window) or
                    (p["side"] == "R" and 0 <= p["end"] - j["position"] <= window)):
                target = j
                break
        if target is None:
            loose.append(p)
        else:
            target["pairs"].append(p)
    loose.sort(key=lambda p: (p["chrom"], p["side"], p["start"]))
    groups, group = [], []
    for p in loose:
        if group and (p["chrom"], p["side"]) == (group[-1]["chrom"], group[-1]["side"]) \
                and p["start"] - group[-1]["start"] <= window:
            group.append(p)
        else:
            if group:
                groups.append(group)
            group = [p]
    if group:
        groups.append(group)
    for group in groups:
        if len(group) < min_pairs:
            continue
        side = group[0]["side"]
        junctions.append({
            "chrom": group[0]["chrom"], "side": side, "exact": False, "clips": [],
            "pairs": group,
            "position": max(p["end"] for p in group) if side == "L" else min(p["start"] for p in group)})
    return junctions


def make_sites(junctions, span):
    """Junctions of one chromosome within span bases of each other, in order."""
    junctions.sort(key=lambda j: (j["chrom"], j["position"]))
    sites, site = [], []
    for j in junctions:
        if site and j["chrom"] == site[-1]["chrom"] and j["position"] - site[-1]["position"] <= span:
            site.append(j)
        else:
            if site:
                sites.append(site)
            site = [j]
    if site:
        sites.append(site)
    return sites


def site_pieces(site, pieces, index, loci, near, window=5000):
    """For long sequences cut into pieces: the matches of the pieces aligned to
    an element that come from the same sequence as a junction piece of the
    site, within window bases of it."""
    places = {}
    for j in site:
        for c in j["clips"]:
            origin = piece_origin(c["name"])
            if origin:
                places.setdefault(origin[0], []).append(origin[1])
    out = []
    for source, starts in places.items():
        lo, hi = min(starts) - window, max(starts) + window
        for start, seq in pieces.get(source, ()):
            if lo <= start <= hi:
                hit = {eid: (v[0], "") for eid, v in match(seq, index).items()
                       if not near_native(eid, site[0]["chrom"], site[0]["position"], loci, near)}
                if hit:
                    out.append(hit)
    return out


def resolve(evidence, inside):
    """Elements at a site from the matches of its evidence ([{element: ...}]).

    An element whose evidence is a strict subset of another's is dropped: its
    matches are sequence it shares with that one (natMX6 against the pieces of
    a kanMX cassette, which also cover the kan gene), unless it is contained in
    it (the TEF promoter inside its cassette). Elements left with the same
    evidence, neither containing the other, are alternatives.
    """
    ev = {}
    for i, hit in enumerate(evidence):
        for eid in hit:
            ev.setdefault(eid, set()).add(i)
    kept = {a: ea for a, ea in ev.items()
            if not any(ea < eb and b not in inside.get(a, ()) for b, eb in ev.items() if b != a)}
    alternatives = {}
    for a, ea in kept.items():
        alternatives[a] = sorted(b for b, eb in kept.items()
                                 if b != a and eb == ea and b not in inside.get(a, ())
                                 and a not in inside.get(b, ()))
    return kept, alternatives


def describe_site(site, names):
    """Place of a site, as reference coordinates (1-based, both ends included),
    a label and a sentence.

    The left end is the outermost junction where the construct follows the
    genome (side L), the right end the outermost one where the genome follows
    the construct (side R). A side seen by read pairs only is placed at the
    nearest mate, the junction lying beyond it within about a fragment length.
    """
    chrom = site[0]["chrom"]
    display = names.get(chrom, chrom)

    def side_position(side):
        for exact in (True, False):
            rows = [j["position"] for j in site if j["side"] == side and j["exact"] == exact]
            if rows:
                return (min(rows) if side == "L" else max(rows)), exact
        return None, None

    a, a_exact = side_position("L")
    b, b_exact = side_position("R")
    reads = sum(len(j["clips"]) for j in site)
    pairs = sum(len(j["pairs"]) for j in site)
    evidence = []
    if reads:
        evidence.append("{} read{} crossing {}".format(
            reads, "s" if reads > 1 else "",
            "the two junctions" if a_exact and b_exact else "the junction"))
    if pairs:
        evidence.append("{} read pair{} with one read in the construct".format(
            pairs, "s" if pairs > 1 else ""))
    evidence = " and ".join(evidence)
    approx = not ((a is None or a_exact) and (b is None or b_exact))
    if a is not None and b is not None and not approx:
        if b > a:
            start, end = a + 1, b
            text = "replaces {}:{}-{} ({} bp of the reference)".format(display, start, end, b - a)
        elif b == a:
            start, end = a, a + 1
            text = "inserted between {}:{} and {}".format(display, a, a + 1)
        else:
            start, end = b + 1, a
            text = ("inserted at {}:{}-{}, these {} bp of the reference being found on both "
                    "sides of it".format(display, start, end, a - b))
        kind = "insertion"
        label = "{}:{}-{}".format(display, start, end)
    elif a is not None and b is not None:
        if b > a:
            start, end = a + 1, b
            text = ("inserted between about {}:{} and {}, replacing part or none of the "
                    "reference between them".format(display, a, b + 1))
        elif b == a:
            start, end = a, a + 1
            text = "inserted between about {}:{} and {}".format(display, a, a + 1)
        else:
            start, end = b + 1, a
            text = "inserted at about {}:{}-{}".format(display, start, end)
        kind = "approximate"
        label = "~{}:{}-{}".format(display, start, end)
    else:
        exact = a_exact if a is not None else b_exact
        if a is not None:
            start = end = a
            text = "after {}{}:{}, the genome then the construct".format(
                "" if exact else "about ", display, a)
        else:
            start = end = b + 1
            text = "before {}{}:{}, the construct then the genome".format(
                "" if exact else "about ", display, b + 1)
        if exact:
            text += "; the junction on the other side was not seen"
        kind = "one junction" if exact else "approximate"
        label = "{}{}:{}".format("" if exact else "~", display, start)
    if kind == "approximate":
        text += ("; a side seen by read pairs only is placed at the nearest mate, the junction "
                 "lying within about a fragment length beyond it")
    return {"chrom": chrom, "chromosome": display, "start": start, "end": end, "kind": kind,
            "label": label, "text": text + ("; " + evidence if evidence else ""),
            "reads": reads, "pairs": pairs}


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--clips", required=True, help="tool_clips.tsv.gz of MarkerAlign")
    parser.add_argument("--reads", required=True, help="tool_reads.tsv of MarkerAlign")
    parser.add_argument("--mates", required=True, help="tool_mates.tsv of MarkerAlign")
    parser.add_argument("--chunks", required=True, help="tool_chunks.tsv of MarkerAlign")
    parser.add_argument("--read-length", help="read_length.txt of MarkerAlign (mean read length)")
    parser.add_argument("--layout", help="read_layout.txt of MarkerAlign: paired or single, then "
                                         "pieces when long sequences were cut into pieces")
    parser.add_argument("--tools", required=True, help="tools.fasta")
    parser.add_argument("--manifest", required=True, help="tools.tsv")
    parser.add_argument("--coverage", required=True, help="tool_coverage.tsv of ToolScore")
    parser.add_argument("--reference", required=True, help="reference genome FASTA, for the "
                                                           "chromosome names")
    parser.add_argument("--seq-kind", choices=("reads", "assembly"), default="reads")
    parser.add_argument("--min-reads", type=int, default=2,
                        help="clipped reads needed for a junction (1 for an assembly, whose "
                             "pieces cross each junction once or twice)")
    parser.add_argument("--min-pairs", type=int, default=3,
                        help="read pairs needed for a junction seen by pairs alone")
    parser.add_argument("--fragment", type=int, default=1000,
                        help="longest fragment expected: distance within which a mate points "
                             "at its junction")
    parser.add_argument("--max-site-span", type=int, default=20000,
                        help="junctions closer than this on one chromosome are one insertion "
                             "site (the longest yeast gene that a cassette replaces is 15 kb)")
    parser.add_argument("--near-native", type=int, default=10000,
                        help="a junction of a native element closer than this to its native "
                             "locus is a change of that locus, not an insertion")
    parser.add_argument("--out", required=True, help="tool_positions.json")
    parser.add_argument("--tsv", required=True, help="tool_positions.tsv")
    args = parser.parse_args()

    tools = score_tools(load_tsv(args.manifest), load_tsv(args.coverage), 1)
    colour = {t["element_id"]: t["colour"] for t in tools}
    plasmid = {t["element_id"]: t.get("source_plasmid") or "" for t in tools}

    def other_name(eid, own):
        """Name of an alternative element, with its plasmid when it shares the
        name of the element it is an alternative to (two sequences of natMX6)."""
        if short(eid) == short(own) and plasmid.get(eid):
            return "{} of {}".format(short(eid), plasmid[eid])
        return short(eid)

    present = [t for t in tools if t["colour"] == "green"]
    present_ids = set(t["element_id"] for t in present)
    elements = read_fasta(args.tools)
    names = chromosome_names(args.reference)
    read_length = None
    if args.read_length:
        try:
            with open(args.read_length) as handle:
                read_length = float(handle.read().strip())
        except (OSError, ValueError):
            read_length = None
    words = []
    if args.layout:
        with open(args.layout) as handle:
            words = handle.read().split()
    layout = words[0] if words else "single"
    cut = "pieces" in words or args.seq_kind == "assembly"
    min_reads = 1 if args.seq_kind == "assembly" else args.min_reads

    index = kmer_index(elements)
    inside = containment(elements)
    loci, mask = native_loci(args.chunks, elements)
    clips, links = read_clips(args.clips, index, elements, mask, loci, args.near_native)
    junctions = clip_junctions(clips, min_reads)
    pairs = []
    if layout == "paired":
        pairs, pair_links = read_mates(args.reads, args.mates, index, mask, loci, args.near_native)
        links += pair_links
        junctions = pair_junctions(pairs, junctions, args.min_pairs, args.fragment)
    pieces = {}
    if cut:
        with open(args.reads) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                origin = piece_origin(fields[0]) if len(fields) >= 6 else None
                if origin:
                    pieces.setdefault(origin[0], []).append((origin[1], fields[5].upper()))

    sites = []
    for site in make_sites(junctions, args.max_site_span):
        evidence = [c["hit"] for j in site for c in j["clips"]] + \
                   [p["hit"] for j in site for p in j["pairs"]]
        if cut:
            evidence += site_pieces(site, pieces, index, loci, args.near_native)
        kept, alternatives = resolve(evidence, inside)
        strands = {}
        for hit in evidence:
            for eid, (strand, _) in hit.items():
                strands.setdefault(eid, []).append(strand)
        listed = {}
        for eid in kept:
            listed[eid] = {"how": "junction", "strand": majority(strands.get(eid, [])),
                           "or": [b for b in alternatives[eid] if colour.get(b) != "red"]}
        for small, holders in inside.items():
            hosts = sorted(h for h in holders if h in kept)
            if small in listed or not hosts:
                continue
            # Part of the construct when every element that may be there holds it.
            maybe = set(hosts) | set(b for h in hosts for b in listed[h]["or"])
            listed[small] = {"how": "inside", "hosts": hosts, "strand": listed[hosts[0]]["strand"],
                             "or": [], "unsure": not maybe <= holders}
        sites.append({"site": site, "listed": listed})

    # A native element joined to an element placed at a single site, without
    # alternatives there, is a copy in that construct. Only the reads that
    # match one construct count: a read of the TEF promoter, which every
    # cassette built on it shares, does not tell which one the native element
    # is joined to.
    placed = {}
    for s in sites:
        for eid, row in s["listed"].items():
            if row["how"] == "junction":
                placed.setdefault(eid, []).append((s, row))
    link_count = {}
    for native, hit in links:
        tops = [e for e in hit if not any(h in hit for h in inside.get(e, ()))]
        if len(tops) == 1:
            link_count[(native, tops[0])] = link_count.get((native, tops[0]), 0) + 1
    for (native, eid), n in sorted(link_count.items()):
        where = placed.get(eid, [])
        if n >= min_reads and len(where) == 1 and not where[0][1]["or"] \
                and native not in where[0][0]["listed"]:
            where[0][0]["listed"][native] = {"how": "joined", "hosts": [eid], "strand": "",
                                             "or": []}

    out_sites = []
    for s in sites:
        shown = {eid: row for eid, row in s["listed"].items() if eid in present_ids}
        if not shown:
            continue
        place = describe_site(s["site"], names)
        place["id"] = len(out_sites) + 1
        place["junctions"] = [{"position": j["position"], "side": j["side"], "exact": j["exact"],
                               "reads": len(j["clips"]), "pairs": len(j["pairs"])}
                              for j in s["site"]]
        place["elements"] = shown
        out_sites.append(place)

    notes = []
    if args.seq_kind == "reads" and not cut and read_length is not None and read_length < 50:
        notes.append("reads of {:.0f} bases on average: too short to hold the 15 clipped bases of "
                     "a junction beside a mappable part, so {}".format(
                         read_length, "only the read pairs place the elements, within about a "
                         "fragment length" if layout == "paired" else
                         "the elements cannot be placed (unpaired reads)"))
    out_elements, rows = {}, []
    for t in present:
        eid = t["element_id"]
        positions = []
        for site in out_sites:
            row = site["elements"].get(eid)
            if row is None:
                continue
            detail = site["text"]
            alternatives = row.get("or") or []
            if row["how"] == "inside":
                detail = "part of {}{}; ".format(
                    " or ".join(short(h) for h in row["hosts"]),
                    ", if that cassette is the one here" if row.get("unsure") else "") + detail
            elif row["how"] == "joined":
                detail = "a copy joined to {} here; ".format(short(row["hosts"][0])) + detail
            if alternatives:
                detail += ("; or {}, in its place: the reads seen here are in sequence they "
                           "share".format(", ".join(other_name(b, eid) for b in alternatives)))
            if row.get("strand"):
                detail += "; element on the {} strand of the reference".format(
                    "forward" if row["strand"] == "+" else "reverse")
            positions.append({"kind": site["kind"], "site": site["id"], "chrom": site["chrom"],
                              "chromosome": site["chromosome"], "start": site["start"],
                              "end": site["end"], "strand": row.get("strand") or "",
                              "label": site["label"], "reads": site["reads"],
                              "pairs": site["pairs"],
                              "alternatives": [other_name(b, eid) for b in alternatives],
                              "detail": detail})
        locus = loci.get(eid)
        if locus:
            display = names.get(locus["chrom"], locus["chrom"])
            positions.append({
                "kind": "native", "site": None, "chrom": locus["chrom"], "chromosome": display,
                "start": locus["start"] + 1, "end": locus["end"], "strand": locus["strand"],
                "label": "{}:{}-{}".format(display, locus["start"] + 1, locus["end"]),
                "reads": 0, "pairs": 0, "alternatives": [],
                "detail": "native locus: the reference genome carries this sequence there{}{}".format(
                    " (over {:.0f} % of the element)".format(100 * locus["fraction"])
                    if locus["fraction"] < 0.95 else "",
                    ", and elsewhere too (repeated sequence)" if locus["repeated"] else "")})
        note = ""
        if not positions:
            if args.seq_kind == "assembly":
                note = ("no piece of the assembly joins it to the genome: it may be on a contig of "
                        "its own (a plasmid, or the assembly broken at the construct)")
            elif read_length is not None and read_length < 50 and layout != "paired" and not cut:
                note = "not placed: the reads are too short to cross a junction, and unpaired"
            else:
                note = ("no read{} joins it to the genome: it may be on a plasmid, or inside a "
                        "longer construct whose ends are not in the database".format(
                            " nor read pair" if layout == "paired" else ""))
        out_elements[eid] = {"name": t["name"], "origin": t["origin"], "positions": positions,
                             "note": note}
        for p in positions or [{}]:
            rows.append([eid, t["name"], t["origin"], p.get("kind", "not found"),
                         p.get("chromosome", ""), p.get("chrom", ""), p.get("start", ""),
                         p.get("end", ""), p.get("strand", ""), p.get("reads", ""),
                         p.get("pairs", ""), p.get("label", ""),
                         ",".join(p.get("alternatives") or []), p.get("detail", note)])

    for site in out_sites:
        site["elements"] = {eid: {"how": row["how"], "strand": row.get("strand") or "",
                                  "hosts": row.get("hosts") or [],
                                  "alternatives": [short(b) for b in row.get("or") or []]}
                            for eid, row in site["elements"].items()}
    result = {"searched": True, "seq_kind": args.seq_kind, "read_layout": layout,
              "read_length": read_length, "notes": notes, "elements": out_elements,
              "sites": out_sites}
    with open(args.out, "w") as handle:
        json.dump(result, handle, indent=1)
    with open(args.tsv, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["element_id", "name", "origin", "position_kind", "chromosome",
                         "reference_sequence", "start", "end", "strand", "junction_reads",
                         "read_pairs", "position", "or_instead", "detail"])
        writer.writerows(rows)

    print("{} element(s) present, {} insertion site(s) holding one, {} clipped junction read(s), "
          "{} read pair(s) with one read in an element".format(
              len(present), len(out_sites), len(clips), len(pairs)))
    for eid, e in out_elements.items():
        print("{}: {}".format(eid, "; ".join(
            p["label"] + " (" + p["kind"] + (", or " + "/".join(p["alternatives"])
                                             if p["alternatives"] else "") + ")"
            for p in e["positions"]) or e["note"]))
    for note in notes:
        print("note: " + note)
    return 0


if __name__ == "__main__":
    sys.exit(main())
