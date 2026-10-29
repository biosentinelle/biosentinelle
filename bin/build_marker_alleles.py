#!/usr/bin/env python3
"""
Build the native marker allele catalogue from strain assemblies.

For every assembly, each marker locus (build_marker_loci.py) is aligned with
minimap2 and the strain's allele is read off the alignment: substitutions and
small indels from the cs tag, and the larger events from the way the alignment
breaks up,

    locus gap, no strain sequence     deletion (BY4742 ura3D0, his3D1)
    strain sequence, no locus gap     insertion (YPH499 ura3-52, a Ty1)
    both                              replacement (a cassette swapped in)
    one flank only                    deletion running past the locus

Every allele then goes through marker_effects.classify, which judges it by the
protein it encodes. This is what lets the pipeline tell apart, at a fully
covered gene, a strain's own neutral polymorphisms from a mutation put there
to inactivate it.

The catalogue records, per gene, every protein form seen and which strains
carry it. Most of the assemblies are natural and industrial isolates, which
are overwhelmingly prototrophic, so a missense form seen in them is evidence
that the form is a natural active variant. The documented genotypes
(documented_genotypes.tsv, each with its citation) name the laboratory
alleles, and are checked against what the assembly actually carries: a
disagreement is reported, never silently resolved in favour of either.

Assemblies are chosen one per strain: asm01 HP0 when it exists, otherwise both
haplotypes of a phased assembly, otherwise the collapsed one.
"""

import argparse
import concurrent.futures
import csv
import hashlib
import os
import re
import shlex
import subprocess
import sys

from marker_effects import classify, read_fasta, revcomp

CS_TOKEN = re.compile(r"(:[0-9]+|\*[a-z][a-z]|[+-][a-z]+|~[a-z]{2}[0-9]+[a-z]{2})")


def select_assemblies(directory):
    """One assembly per strain, as (label, path)."""
    by_strain = {}
    for name in sorted(os.listdir(directory)):
        if not name.endswith(".fa") and not name.endswith(".fasta"):
            continue
        parts = name.split(".")
        if len(parts) < 4:
            continue
        strain, asm, hap = parts[0], parts[1], parts[2]
        by_strain.setdefault(strain, []).append((asm, hap, os.path.join(directory, name)))
    chosen = []
    for strain, files in sorted(by_strain.items()):
        files.sort()
        first_asm = files[0][0]
        files = [f for f in files if f[0] == first_asm]
        haps = {hap: path for _, hap, path in files}
        if "HP0" in haps:
            chosen.append((strain, haps["HP0"]))
        elif "HP1" in haps and "HP2" in haps:
            chosen.append((strain + ".HP1", haps["HP1"]))
            chosen.append((strain + ".HP2", haps["HP2"]))
        elif "HP" in haps:
            chosen.append((strain, haps["HP"]))
        elif "collapsed" in haps:
            chosen.append((strain, haps["collapsed"]))
    return chosen


def run_minimap2(command, loci_fasta, assembly, paf):
    if os.path.exists(paf):  # written atomically, so an empty file is a real result
        return paf
    tmp = paf + ".tmp"
    with open(tmp, "w") as handle:
        subprocess.run(shlex.split(command) + ["-c", "-x", "asm5", "--cs", "-t", "2",
                                               loci_fasta, assembly],
                       stdout=handle, stderr=subprocess.DEVNULL, check=True)
    os.replace(tmp, paf)
    return paf


def parse_paf(path):
    alignments = []
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            tags = dict((f[:2], f[5:]) for f in fields[12:])
            alignments.append({
                "qname": fields[0], "qlen": int(fields[1]),
                "qs": int(fields[2]), "qe": int(fields[3]), "strand": fields[4],
                "tname": fields[5], "tlen": int(fields[6]),
                "ts": int(fields[7]), "te": int(fields[8]),
                "mapq": int(fields[11]), "cs": tags.get("cs", ""),
                "primary": tags.get("tp", "P") == "P",
            })
    return alignments


def cs_events(alignment):
    """Small events of one alignment, in locus coordinates."""
    events, pos = [], alignment["ts"]
    for token in CS_TOKEN.findall(alignment["cs"]):
        kind = token[0]
        if kind == ":":
            pos += int(token[1:])
        elif kind == "*":
            events.append((pos, pos + 1, token[2].upper()))
            pos += 1
        elif kind == "+":
            events.append((pos, pos, token[1:].upper()))
        elif kind == "-":
            events.append((pos, pos + len(token) - 1, ""))
            pos += len(token) - 1
    return events


def merge_intervals(intervals):
    out = []
    for start, end in sorted(intervals):
        if out and start <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], end))
        else:
            out.append((start, end))
    return out


def unique_bases(block, mask):
    """Locus bases of one alignment block outside the repeat mask."""
    covered = block["te"] - block["ts"]
    for start, end in mask:
        covered -= max(0, min(end, block["te"]) - max(start, block["ts"]))
    return covered


def repeat_mask(reference_alignments, loci):
    """Locus intervals that also occur elsewhere in the reference genome.

    The reference is aligned to the loci like any assembly. The block on the
    locus' own chromosome and coordinates is the locus itself; every other
    block is a second copy of part of it, typically a Ty LTR or delta in a
    flank, which also hits dozens of places in every other genome.
    """
    mask = {}
    for locus_id, locus in loci.items():
        g0, g1 = int(locus["genome_start"]), int(locus["genome_end"])
        intervals = []
        for a in reference_alignments:
            if a["tname"] != locus_id:
                continue
            if a["qname"] == locus["chrom"] and a["qs"] < g1 and a["qe"] > g0:
                continue
            intervals.append((a["ts"], a["te"]))
        mask[locus_id] = merge_intervals(intervals)
    return mask


def best_chain(alignments, locus_len, mask):
    """The blocks of the one strain contig and strand that carry the locus.

    The anchor is the block with the most sequence outside the repeat mask, so
    a flank repeat cannot decide where the gene is. Other blocks join it only if
    they sit where collinearity puts them, give or take an insertion: a block
    that is all repeat is held to 10 kb of its expected place, which keeps a
    Ty1 insertion and excludes a delta element from the next gene over.
    """
    if not alignments:
        return []
    anchor = max(alignments, key=lambda b: unique_bases(b, mask))
    if unique_bases(anchor, mask) < 50:
        return []

    def expected_offset(block):
        if anchor["strand"] == "+":
            return block["qs"] - (anchor["qs"] + block["ts"] - anchor["ts"])
        return block["qe"] - (anchor["qe"] - (block["ts"] - anchor["ts"]))

    blocks = []
    for block in alignments:
        if block["qname"] != anchor["qname"] or block["strand"] != anchor["strand"]:
            continue
        limit = 10000 if unique_bases(block, mask) < 50 else 3 * locus_len + 50000
        if abs(expected_offset(block)) <= limit:
            blocks.append(block)
    blocks.sort(key=lambda b: (b["ts"], -b["te"]))
    kept = []
    for block in blocks:
        if kept and block["te"] <= kept[-1]["te"]:
            continue  # contained in a block already kept
        kept.append(block)
    return kept


def allele_events(blocks, locus_len, contig_seq):
    """Small and structural events, plus the locus intervals with no data."""
    events, uncovered, structural = [], [], []
    if not blocks:
        return events, [(0, locus_len)], ["locus not found in the assembly"]
    if blocks[0]["ts"] > 0:
        uncovered.append((0, blocks[0]["ts"]))
    for i, block in enumerate(blocks):
        start_cut = blocks[i - 1]["te"] if i else block["ts"]
        for event in cs_events(block):
            if event[0] >= start_cut or i == 0:
                events.append(event)
        if i == 0:
            continue
        prev = blocks[i - 1]
        overlap = max(0, prev["te"] - block["ts"])
        locus_gap = block["ts"] - prev["te"]
        if block["strand"] == "+":
            strain_gap = block["qs"] - prev["qe"]
            q0, q1 = prev["qe"], block["qs"]
        else:
            strain_gap = prev["qs"] - block["qe"]
            q0, q1 = block["qe"], prev["qs"]
        strain_gap += overlap  # duplicated target bases, such as a Ty1 target site
        inserted = ""
        if strain_gap > 0 and contig_seq is not None and q1 > q0:
            inserted = contig_seq[q0:q1]
            if block["strand"] == "-":
                inserted = revcomp(inserted)
        if strain_gap < 0:
            strain_gap = 0
        if locus_gap > 0:
            events.append((prev["te"], block["ts"], inserted))
            structural.append("%s of %d bp at locus %d-%d%s"
                              % ("replacement" if inserted else "deletion", locus_gap,
                                 prev["te"], block["ts"],
                                 ", %d bp inserted" % len(inserted) if inserted else ""))
        elif strain_gap > 0:
            events.append((prev["te"], prev["te"], inserted or None))
            structural.append("insertion of %d bp at locus %d" % (strain_gap, prev["te"]))
    if blocks[-1]["te"] < locus_len:
        uncovered.append((blocks[-1]["te"], locus_len))
    return events, uncovered, structural


def kmer_set(seq, k):
    out = set()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        out.add(min(kmer, revcomp(kmer)))
    return out


def name_insertion(seq, tools, k=21):
    """Foreign tool elements contained in an inserted sequence, longest first."""
    if not seq or len(seq) < k:
        return ""
    kmers = kmer_set(seq, k)
    hits = []
    for name, element_kmers in tools.items():
        if element_kmers:
            fraction = len(kmers & element_kmers) / float(len(element_kmers))
            if fraction >= 0.5:
                hits.append((len(element_kmers), name.split("__")[0], fraction))
    seen, out = set(), []
    for _, name, fraction in sorted(hits, reverse=True):
        if name not in seen:
            seen.add(name)
            out.append("%s %.2f" % (name, fraction))
    return "contains " + ", ".join(out) if out else ""


def short_hash(text):
    return hashlib.sha1(text.encode()).hexdigest()[:8]


def encode_events(events):
    out = []
    for start, end, alt in events:
        if alt is None:
            alt = "?"
        elif len(alt) > 50:
            alt = "INS%d" % len(alt)
        out.append("%d:%d:%s" % (start, end, alt))
    return ";".join(out)


def documented_alleles(path, genes):
    """strain -> {gene: documented allele token}, wild type when not listed."""
    out, sources = {}, {}
    if not path:
        return out, sources
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            tokens = row["documented_genotype"].split()
            alleles = {}
            for gene in genes:
                hits = [t for t in tokens if t.lower().startswith(gene.lower())]
                mutant = [t for t in hits if not t.startswith(gene.upper())]
                alleles[gene] = mutant[0] if mutant else gene.upper()
            out[row["strain"]] = alleles
            sources[row["strain"]] = row["source"]
    return out, sources


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--loci", required=True, help="marker_loci.tsv")
    parser.add_argument("--loci-fasta", required=True, help="marker_loci.fasta")
    parser.add_argument("--assembly-dir", required=True)
    parser.add_argument("--reference", required=True,
                        help="reference genome FASTA, aligned to the loci to find their repeats")
    parser.add_argument("--repeats-out", help="marker_repeats.tsv, the repeat mask written out")
    parser.add_argument("--documented", help="documented_genotypes.tsv")
    parser.add_argument("--tools", help="tools.fasta, to name inserted sequences")
    parser.add_argument("--tools-manifest", help="tools.tsv, only foreign elements are used")
    parser.add_argument("--paf-dir", required=True, help="where the PAF files are cached")
    parser.add_argument("--minimap2", default="minimap2",
                        help="minimap2 command, e.g. 'apptainer exec img minimap2'")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--max-absent", type=int, default=4,
                        help="exclude an assembly in which this many loci are not found")
    parser.add_argument("--out", required=True, help="marker_alleles.tsv")
    parser.add_argument("--report", help="text summary")
    args = parser.parse_args()

    with open(args.loci, newline="") as handle:
        loci = {row["locus_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    locus_seqs = dict((name.split()[0], seq) for name, seq in read_fasta(args.loci_fasta))
    genes = [row["gene"] for row in loci.values()]
    documented, doc_sources = documented_alleles(args.documented, genes)

    tools = {}
    if args.tools:
        foreign = None
        if args.tools_manifest:
            with open(args.tools_manifest, newline="") as handle:
                foreign = set(r["element_id"] for r in csv.DictReader(handle, delimiter="\t")
                              if r.get("origin") == "foreign")
        for name, seq in read_fasta(args.tools):
            element = name.split()[0]
            if foreign is None or element in foreign:
                tools[element] = kmer_set(seq, 21)

    assemblies = select_assemblies(args.assembly_dir)
    os.makedirs(args.paf_dir, exist_ok=True)
    sys.stderr.write("%d assemblies selected\n" % len(assemblies))

    def job(item):
        label, path = item
        return label, path, run_minimap2(args.minimap2, args.loci_fasta, path,
                                         os.path.join(args.paf_dir, label + ".paf"))

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        done = 0
        results = []
        for result in pool.map(job, assemblies):
            results.append(result)
            done += 1
            if done % 20 == 0:
                sys.stderr.write("  %d of %d aligned\n" % (done, len(assemblies)))

    reference_paf = run_minimap2(args.minimap2, args.loci_fasta, args.reference,
                                 os.path.join(args.paf_dir, "_reference.paf"))
    mask = repeat_mask(parse_paf(reference_paf), loci)
    if args.repeats_out:
        with open(args.repeats_out, "w", newline="\n") as handle:
            handle.write("locus_id\tstart\tend\n")
            for locus_id in loci:
                for start, end in mask[locus_id]:
                    handle.write("%s\t%d\t%d\n" % (locus_id, start, end))

    rows = []
    for label, path, paf in results:
        alignments = parse_paf(paf)
        strain = label.split(".")[0]
        contigs = None
        for locus_id, locus in loci.items():
            locus_seq = locus_seqs[locus_id]
            locus_len = len(locus_seq)
            orf_start, orf_end = int(locus["orf_start"]), int(locus["orf_end"])
            blocks = best_chain([a for a in alignments if a["tname"] == locus_id], locus_len,
                                mask[locus_id])
            needs_seq = any(abs(blocks[i]["qs"] - blocks[i - 1]["qe"]) > 0
                            for i in range(1, len(blocks)))
            contig_seq = None
            if needs_seq:
                if contigs is None:
                    contigs = dict((n.split()[0], s) for n, s in read_fasta(path))
                contig_seq = contigs.get(blocks[0]["qname"])
            events, uncovered, structural = allele_events(blocks, locus_len, contig_seq)

            orf_uncovered = sum(max(0, min(e, orf_end) - max(s, orf_start))
                                for s, e in uncovered)
            if not blocks:
                effect = {"effect": "absent", "activity": "inactive", "protein": None,
                          "protein_change": "locus not found in the assembly",
                          "coding": [], "upstream": [], "downstream": []}
            elif orf_uncovered >= 0.5 * (orf_end - orf_start):
                effect = {"effect": "orf_deleted", "activity": "inactive", "protein": None,
                          "protein_change": "ORF missing, deletion extends past the locus flank",
                          "coding": [], "upstream": [], "downstream": []}
            else:
                effect = classify(locus_seq, orf_start, orf_end, events, uncovered)

            insert_names = []
            for start, end, alt in events:
                if alt and len(alt) >= 100:
                    named = name_insertion(alt, tools)
                    insert_names.append(named or "unidentified %d bp" % len(alt))

            doc = documented.get(strain, {}).get(locus["gene"], "")
            agreement = ""
            if doc:
                expected_inactive = doc != locus["gene"].upper()
                activity = effect["activity"]
                if activity == "uncertain":
                    agreement = "unconfirmed"
                elif (activity == "inactive") == expected_inactive:
                    agreement = "agree"
                else:
                    agreement = "DISAGREE"

            protein = effect.get("protein")
            rows.append({
                "assembly": label,
                "strain": strain,
                "gene": locus["gene"],
                "locus_id": locus_id,
                "allele_id": "%s-%s" % (locus["gene"],
                                        short_hash(encode_events(events)) if events else "ref"),
                "protein_id": ("%s-p%s" % (locus["gene"], short_hash(protein))) if protein else "",
                "effect": effect["effect"],
                "activity": effect["activity"],
                "protein_change": effect.get("protein_change", ""),
                "coding_events": ",".join(effect.get("coding", [])),
                "upstream_events": ",".join(effect.get("upstream", [])),
                "structural": "; ".join(structural),
                "inserted_element": "; ".join(insert_names),
                "uncovered_bp": sum(e - s for s, e in uncovered),
                "uncovered": ";".join("%d-%d" % (s, e) for s, e in uncovered),
                "events": encode_events(events),
                "documented_allele": doc,
                "documented_agreement": agreement,
                "contig": blocks[0]["qname"] if blocks else "",
                "contig_start": min(b["qs"] for b in blocks) if blocks else "",
                "contig_end": max(b["qe"] for b in blocks) if blocks else "",
                "assembly_file": os.path.basename(path),
            })

    # An assembly that misses most loci outright is not a divergent S. cerevisiae
    # strain, it is another species: ScRAPdb also holds S. paradoxus and others,
    # too far from the reference for asm5. Their "absent" calls would read as
    # deletions, so they are left out and listed.
    absent = {}
    for r in rows:
        if r["effect"] == "absent":
            absent[r["assembly"]] = absent.get(r["assembly"], 0) + 1
    excluded = sorted(a for a, n in absent.items() if n >= args.max_absent)
    rows = [r for r in rows if r["assembly"] not in excluded]

    fields = list(rows[0])
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    lines = ["%d assemblies, %d loci, %d alleles" % (len(results) - len(excluded), len(loci),
                                                      len(rows)),
             "%d excluded, absent at %d or more loci (another species): %s"
             % (len(excluded), args.max_absent, " ".join(excluded)), ""]
    for gene in genes:
        sub = [r for r in rows if r["gene"] == gene]
        tally = {}
        for r in sub:
            tally[r["effect"]] = tally.get(r["effect"], 0) + 1
        forms = set(r["protein_id"] for r in sub if r["protein_id"])
        lines.append("%-6s %d protein forms; %s" % (
            gene, len(forms), ", ".join("%s=%d" % kv for kv in sorted(tally.items()))))
    lines.append("")
    lines.append("repeat mask (locus bases also found elsewhere in the reference):")
    for locus_id in loci:
        lines.append("  %-16s %s" % (locus_id, ", ".join("%d-%d" % iv for iv in mask[locus_id])
                                     or "none"))
    lines.append("")
    lines.append("documented strains:")
    for r in rows:
        if r["documented_allele"]:
            lines.append("  %-12s %-6s documented %-14s found %-16s %-9s %s %s" % (
                r["assembly"], r["gene"], r["documented_allele"], r["effect"],
                r["documented_agreement"], r["protein_change"], r["structural"]))
    text = "\n".join(lines) + "\n"
    sys.stderr.write(text)
    if args.report:
        with open(args.report, "w", newline="\n") as handle:
            handle.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
