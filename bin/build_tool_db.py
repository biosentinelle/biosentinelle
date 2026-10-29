#!/usr/bin/env python3
"""
Build the biomolecular tool element database from a curated feature package.

Reads the feature manifest and multi-FASTA produced for the yeast integration
cassette package, deduplicates identical sequences, assigns a stable id and a
functional class to each element, and writes db/<organism>/tools/.

No sequence is created here: every record is copied from the input FASTA and
keeps a pointer back to its source plasmid and source URL.
"""

import argparse
import csv
import hashlib
import re
import sys
from collections import OrderedDict


# Functional classes. The first pattern that matches a feature label wins, so
# the more specific patterns are listed first.
CLASS_PATTERNS = [
    ("promoter", r"promoter$"),
    ("terminator", r"terminator$"),
    ("selection_marker_drug", r"^(kanmx|natmx|hphmx|ble|kanr|natr|hygr|nrsr|hph|bler)"),
    ("selection_marker_auxo", r"^(his3mx|his3|leu2|ura3|trp1|ade2|lys2|met15|klura3|kltrp1)"),
    ("reporter", r"(gfp|yegfp|mcherry|rfp|yfp|cfp|venus|vn173|vc155|luciferase|lacz)"),
    ("epitope_tag", r"(ha|myc|flag|tap|gst|pk9|v5|his6|6xhis|strep)"),
    ("degron", r"(aid|degron|tir1|dhfr)"),
    ("recombination_site", r"(loxp|frt|rox|att[bprl])"),
]

# Short elements cannot be scored by breadth of coverage: a 27 bp HA tag is
# shorter than a read, is repeated in tandem in real constructs, and matches by
# chance far more often than a kilobase cassette. They stay in the database but
# are flagged so the report can present them separately instead of pretending
# the same score means the same thing.
MIN_SCORABLE_LENGTH = 100


def classify(label, feature_type):
    """Assign a functional class from the feature label, falling back on type."""
    key = label.lower().replace("_", "").replace("-", "").replace(" ", "")
    for name, pattern in CLASS_PATTERNS:
        if re.search(pattern, key):
            return name
    if feature_type in ("promoter", "terminator"):
        return feature_type
    return "other"


def read_fasta(path):
    """Yield (header, sequence) pairs. Sequences are upper-cased and unwrapped."""
    header, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            else:
                chunks.append(line.upper())
    if header is not None:
        yield header, "".join(chunks)


def parse_header(header):
    """Pull the key=value fields out of a feature-package FASTA header."""
    parts = header.split("|")
    fields = {"source_plasmid": parts[0]}
    for part in parts[1:]:
        if "=" in part:
            key, value = part.split("=", 1)
            fields[key] = value
    return fields


def slugify(text):
    return re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="source feature manifest TSV")
    parser.add_argument("--fasta", required=True, help="source multi-FASTA")
    parser.add_argument("--outdir", required=True, help="output tools/ directory")
    parser.add_argument("--min-length", type=int, default=20,
                        help="skip elements shorter than this, too short to score reliably")
    args = parser.parse_args()

    # Source URLs are keyed by (plasmid, label) so each element keeps its citation.
    sources = {}
    with open(args.manifest, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["source_plasmid"], slugify(row["feature_label"]))
            sources[key] = row.get("source_page", "")

    # Deduplicate on sequence: the same kanMX appears in many cassettes, and we
    # want to report "kanMX is present", not the same hit once per plasmid.
    elements = OrderedDict()
    skipped_short = 0
    for header, sequence in read_fasta(args.fasta):
        fields = parse_header(header)
        label = fields.get("feature", "unknown")
        plasmid = fields.get("source_plasmid", "unknown")
        if len(sequence) < args.min_length:
            skipped_short += 1
            continue
        digest = hashlib.sha256(sequence.encode()).hexdigest()[:12]
        if digest in elements:
            elements[digest]["seen_in"].append(plasmid)
            continue
        elements[digest] = {
            "element_id": "%s__%s" % (slugify(label), digest[:6]),
            "name": label,
            "element_class": classify(label, fields.get("type", "")),
            "feature_type": fields.get("type", ""),
            "length_bp": len(sequence),
            "scorable": "yes" if len(sequence) >= MIN_SCORABLE_LENGTH else "short",
            "sequence": sequence,
            "sha256_12": digest,
            "source_plasmid": plasmid,
            "source_page": sources.get((plasmid, slugify(label)), ""),
            "seen_in": [plasmid],
        }

    import os
    os.makedirs(args.outdir, exist_ok=True)
    fasta_out = os.path.join(args.outdir, "tools.fasta")
    tsv_out = os.path.join(args.outdir, "tools.tsv")

    with open(fasta_out, "w") as handle:
        for element in elements.values():
            handle.write(">%s class=%s name=%s length=%d scorable=%s\n" % (
                element["element_id"], element["element_class"],
                element["name"], element["length_bp"], element["scorable"]))
            for i in range(0, len(element["sequence"]), 70):
                handle.write(element["sequence"][i:i + 70] + "\n")

    columns = ["element_id", "name", "element_class", "feature_type", "length_bp",
               "scorable", "sha256_12", "source_plasmid", "source_page", "seen_in"]
    with open(tsv_out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for element in elements.values():
            row = dict(element)
            row["seen_in"] = ";".join(sorted(set(row["seen_in"])))
            writer.writerow(row)

    by_class = {}
    for element in elements.values():
        by_class[element["element_class"]] = by_class.get(element["element_class"], 0) + 1

    print("wrote %d unique elements to %s" % (len(elements), fasta_out))
    if skipped_short:
        print("skipped %d elements shorter than %d bp" % (skipped_short, args.min_length))
    for name in sorted(by_class):
        print("  %-24s %d" % (name, by_class[name]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
