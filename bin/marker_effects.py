#!/usr/bin/env python3
"""
Shared logic for native marker gene genotyping.

A marker gene is held as a locus: upstream flank, ORF, downstream flank, always
written in the gene's coding orientation, so the ORF sits at the same offset
whatever strand the gene is on. Every difference from the reference, whether it
was read off a strain assembly or called from sample reads, is expressed as an
event on that locus,

    (start, end, alt)    locus[start:end] is replaced by alt, 0-based, half open

so a substitution is (p, p + 1, "A"), a deletion is (s, e, ""), an insertion
of known sequence is (p, p, "ACG"), and an insertion whose sequence is not
known, which is what short reads give at a transposon or cassette junction, is
(p, p, None).

The same events then go through the same classification, so the database and
the sample are judged by identical rules. The classification asks one
question: what protein does the mutated locus encode, compared with the
reference one. That is what separates the three situations a coverage test
cannot tell apart:

  * a functional gene carrying the strain's own neutral polymorphisms
    (synonymous, or missense forms also seen in active strains),
  * a gene inactivated in place by a point mutation or a small indel
    (nonsense, frameshift, start loss), with full coverage,
  * a gene deleted outright, or interrupted by an insertion.

A missense change is never declared inactivating or harmless on its own: the
protein form has to be recognised in the allele catalogue for that.

Python 3.7 compatible on purpose: this runs in the pipeline's python container.
"""

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

CODON_TABLE = {}
_BASES = "TCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _a in enumerate(_BASES):
    for _j, _b in enumerate(_BASES):
        for _k, _c in enumerate(_BASES):
            CODON_TABLE[_a + _b + _c] = _AMINO[16 * _i + 4 * _j + _k]

# Effects, from the most to the least certain loss of function.
INACTIVE_EFFECTS = ("absent", "orf_deleted", "insertion", "large_deletion",
                    "start_lost", "frameshift", "nonsense")
UNCERTAIN_EFFECTS = ("missense", "inframe_indel", "late_truncation", "stop_lost",
                     "partial")
ACTIVE_EFFECTS = ("reference", "synonymous")

# A premature stop this close to the natural one may leave a working protein.
LATE_TRUNCATION_FRACTION = 0.95

# An in-frame deletion removing at least this share of the protein is called
# inactivating. Smaller ones stay uncertain: a few residues lost from a loop can
# be harmless, a fifth of the protein essentially never is.
LARGE_INFRAME_DELETION_FRACTION = 0.2


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def translate(seq, to_stop=True):
    """Translate from the first base. Stops at the first stop codon (kept as *)."""
    protein = []
    for i in range(0, len(seq) - 2, 3):
        amino = CODON_TABLE.get(seq[i:i + 3].upper(), "X")
        protein.append(amino)
        if to_stop and amino == "*":
            break
    return "".join(protein)


def read_fasta(path):
    name, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name, chunks = line[1:], []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def genome_event_to_locus(event, locus):
    """Convert a (start, end, alt) event in genome coordinates to locus ones.

    locus carries genome_start (0-based), genome_end and strand ("+" or "-").
    """
    start, end, alt = event
    g0, g1 = locus["genome_start"], locus["genome_end"]
    if locus["strand"] == "+":
        return (start - g0, end - g0, alt)
    return (g1 - end, g1 - start, revcomp(alt) if alt else alt)


def cdna_position(locus_pos, orf_start, orf_end):
    """HGVS-like coding position of a 0-based locus position (no c.0)."""
    if locus_pos < orf_start:
        return "-%d" % (orf_start - locus_pos)
    if locus_pos >= orf_end:
        return "*%d" % (locus_pos - orf_end + 1)
    return str(locus_pos - orf_start + 1)


def describe_event(event, locus_seq, orf_start, orf_end):
    """Coding DNA notation of one event, close to HGVS."""
    start, end, alt = event

    def pos(p):
        return cdna_position(p, orf_start, orf_end)

    ref = locus_seq[start:end]
    if alt is None:
        return "c.%s_%sins(unknown)" % (pos(start - 1), pos(start))
    if end - start == 1 and len(alt) == 1:
        return "c.%s%s>%s" % (pos(start), ref, alt)
    if not alt:
        if end - start == 1:
            return "c.%sdel" % pos(start)
        return "c.%s_%sdel(%dbp)" % (pos(start), pos(end - 1), end - start)
    if start == end:
        label = alt if len(alt) <= 12 else "%dbp" % len(alt)
        return "c.%s_%sins%s" % (pos(start - 1), pos(start), label)
    label = alt if len(alt) <= 12 else "%dbp" % len(alt)
    return "c.%s_%sdelins%s" % (pos(start), pos(end - 1), label)


def apply_events(seq, events):
    """Apply non-overlapping events. Returns the new sequence and an offset map."""
    out, cursor = [], 0
    for start, end, alt in sorted(events, key=lambda e: (e[0], e[1])):
        if start < cursor:
            # Overlapping calls cannot both be true; keep the first.
            continue
        out.append(seq[cursor:start])
        out.append(alt or "")
        cursor = end
    out.append(seq[cursor:])
    return "".join(out)


def shifted(position, events):
    """Where a reference locus position lands after the events are applied."""
    delta = 0
    for start, end, alt in events:
        if end <= position:
            delta += len(alt or "") - (end - start)
    return position + delta


def protein_changes(ref_protein, alt_protein):
    """Substitutions between two proteins of equal length, as p.G234D."""
    return ["p.%s%d%s" % (r, i + 1, a)
            for i, (r, a) in enumerate(zip(ref_protein, alt_protein)) if r != a]


def classify(locus_seq, orf_start, orf_end, events, uncovered=None):
    """Effect of a set of events on the ORF.

    uncovered is an optional list of (start, end) locus intervals for which no
    information exists; an ORF that is partly uncovered cannot be called
    reference with confidence and is reported as "partial".

    Returns a dict with effect, activity, protein (the encoded protein, up to and
    excluding the stop), protein_change, and the events split by region.
    """
    ref_orf = locus_seq[orf_start:orf_end]
    ref_protein = translate(ref_orf).rstrip("*")
    orf_len = orf_end - orf_start
    events = sorted(events, key=lambda e: (e[0], e[1]))

    upstream = [e for e in events if e[1] <= orf_start and e[0] < orf_start]
    downstream = [e for e in events if e[0] >= orf_end]
    in_orf = [e for e in events if e not in upstream and e not in downstream]

    result = {
        "effect": None, "activity": None, "protein": None, "protein_change": "",
        "upstream": [describe_event(e, locus_seq, orf_start, orf_end) for e in upstream],
        "downstream": [describe_event(e, locus_seq, orf_start, orf_end) for e in downstream],
        "coding": [describe_event(e, locus_seq, orf_start, orf_end) for e in in_orf],
        "orf_length": orf_len,
    }

    def done(effect, change=""):
        result["effect"] = effect
        result["protein_change"] = change
        if effect in ACTIVE_EFFECTS:
            result["activity"] = "active"
        elif effect in INACTIVE_EFFECTS:
            result["activity"] = "inactive"
        else:
            result["activity"] = "uncertain"
        return result

    # Bases of the ORF removed. A large deletion counts in full whatever was put
    # in its place, a 25 bp scar or a 1.5 kb cassette: new sequence does not give
    # the ORF back. A small replacement is left to the protein comparison.
    deleted = 0
    for start, end, alt in in_orf:
        if alt is not None and (end - start >= 30 or not alt):
            deleted += max(0, min(end, orf_end) - max(start, orf_start))
    if deleted >= 0.5 * orf_len:
        return done("orf_deleted", "%d of %d bp of the ORF deleted" % (deleted, orf_len))

    for start, end, alt in in_orf:
        if alt is None or len(alt) - (end - start) >= 100:
            size = "unknown size" if alt is None else "%d bp" % len(alt)
            return done("insertion", "ORF interrupted at c.%s by an insertion of %s"
                        % (cdna_position(start, orf_start, orf_end), size))

    for start, end, alt in in_orf:
        if start < orf_start + 3 and end > orf_start:
            new_start = shifted(orf_start, upstream)
            mutated_start = apply_events(locus_seq, events)[new_start:new_start + 3]
            if mutated_start != "ATG":
                return done("start_lost", "start codon changed to %s" % (mutated_start or "nothing"))

    if uncovered:
        missing = 0
        for start, end in uncovered:
            missing += max(0, min(end, orf_end) - max(start, orf_start))
        if missing:
            result["missing_bp"] = missing

    if not in_orf:
        if result.get("missing_bp"):
            return done("partial", "%d bp of the ORF have no data" % result["missing_bp"])
        result["protein"] = ref_protein
        return done("reference")

    mutated = apply_events(locus_seq, events)
    new_start = shifted(orf_start, upstream)
    alt_full = translate(mutated[new_start:])
    stopped = alt_full.endswith("*")
    alt_protein = alt_full.rstrip("*")
    result["protein"] = alt_protein

    net = sum(len(alt or "") - (end - start) for start, end, alt in in_orf)
    first_diff = next((i for i, (r, a) in enumerate(zip(ref_protein, alt_protein)) if r != a),
                      min(len(ref_protein), len(alt_protein)))

    if not stopped:
        return done("stop_lost", "no stop codon before the end of the locus")
    if net % 3:
        if first_diff >= LATE_TRUNCATION_FRACTION * len(ref_protein):
            return done("late_truncation", "frameshift at p.%d of %d"
                        % (first_diff + 1, len(ref_protein)))
        return done("frameshift", "p.%s%dfs, protein of %d instead of %d aa"
                    % (ref_protein[first_diff] if first_diff < len(ref_protein) else "*",
                       first_diff + 1, len(alt_protein), len(ref_protein)))
    if net == 0 and len(alt_protein) < len(ref_protein):
        if len(alt_protein) >= LATE_TRUNCATION_FRACTION * len(ref_protein):
            return done("late_truncation", "p.%s%d*, %d of %d aa kept"
                        % (ref_protein[len(alt_protein)], len(alt_protein) + 1,
                           len(alt_protein), len(ref_protein)))
        return done("nonsense", "p.%s%d*, %d of %d aa kept"
                    % (ref_protein[len(alt_protein)], len(alt_protein) + 1,
                       len(alt_protein), len(ref_protein)))
    if net == 0 and len(alt_protein) > len(ref_protein):
        return done("stop_lost", "protein extended to %d aa" % len(alt_protein))
    if net:
        expected = len(ref_protein) + net // 3
        if len(alt_protein) < expected - 1:
            return done("nonsense", "in-frame indel creating a stop, %d of %d aa kept"
                        % (len(alt_protein), len(ref_protein)))
        if -net // 3 >= LARGE_INFRAME_DELETION_FRACTION * len(ref_protein):
            return done("large_deletion", "in-frame deletion of %d of %d aa"
                        % (-net // 3, len(ref_protein)))
        return done("inframe_indel", "%+d aa" % (net // 3))
    changes = protein_changes(ref_protein, alt_protein)
    if not changes:
        return done("synonymous")
    return done("missense", ",".join(changes))
