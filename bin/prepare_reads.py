#!/usr/bin/env python3
"""
Turn any nucleotide sequence input into the reads the pipeline works on.

Accepted, plain or compressed (.gz, .bz2, or inside a .zip):
  FASTQ      .fastq .fq
  FASTA      .fasta .fa .fna .fas .seq .txt   (reads or assembly contigs)
  GenBank    .gb .gbk .gbff                   (the ORIGIN section is read)
  EMBL       .embl .dat                       (the SQ section is read)
The format is recognised from the content, not from the extension. BAM, SAM and
CRAM are converted to FASTQ by samtools before reaching this script.

Every sequence is upper-cased, U becomes T (direct RNA reads) and any other
letter outside ACGT becomes N. A FASTA, GenBank or EMBL record has no base
quality, so it is given a uniform Q40: the aligners need one, and the pipeline
filters on mapping quality, not on base quality.

Two read sets are written, as file name prefixes in align/ and scan_only/:

  align/R1__*, align/R2__*  paired mates, aligned with bowtie2 -1/-2
  align/U__*                unpaired full-length reads, aligned with -U
  align/S__*                shredded reads, aligned with -U
  scan_only/F__*            the full-length sequences behind the S__ files

Bowtie2 is a short-read aligner, so a sequence longer than --max-length (a
long read, a contig) is cut into --frag-length pieces for the alignments only.
The k-mer scan reads every sequence whole: R1, R2, U and F files, never S.
Assembly contigs are cut with half-overlapping pieces, so that each base is in
two pieces and a contig end is not lost; reads are cut without overlap, so the
depth of the alignment stays the depth of the sequencing.

A FASTQ file that needs no change (no U, no long read, no assembly) is hard
linked, or copied when the file system refuses links, rather than rewritten.

--to-fasta mode only converts reference files (FASTA, GenBank, EMBL) to plain
FASTA, for the reference, feature and mating type folders.
"""

import argparse
import bz2
import gzip
import io
import os
import re
import shutil
import sys
import zipfile

SEQ_EXTENSIONS = ("fastq", "fq", "fasta", "fa", "fna", "fas", "seq", "txt",
                  "faa", "gb", "gbk", "gbff", "embl", "dat")
COMPRESSIONS = ("gz", "bz2")
# A mate number right before the extension: x_1.fq, x_R2.fastq.gz, x.1.fq
MATE = re.compile(r"^(.*?)[._]R?([12])$")
CLEAN = str.maketrans({"U": "T"})
NOT_ACGTN = re.compile(r"[^ACGTN]")
SAMPLE_RECORDS = 20000
QUAL = "I"  # Q40 in Phred+33


def stem_of(name):
    """File name without compression and sequence extensions, made path safe."""
    base = os.path.basename(name)
    parts = base.split(".")
    while len(parts) > 1 and parts[-1].lower() in COMPRESSIONS + SEQ_EXTENSIONS:
        parts.pop()
    return re.sub(r"[^A-Za-z0-9._-]", "_", ".".join(parts))


def compression_of(name):
    lower = name.lower()
    if lower.endswith(".gz"):
        return "gz"
    if lower.endswith(".bz2"):
        return "bz2"
    return ""


class Source(object):
    """One sequence file, on disk or inside a zip archive."""

    def __init__(self, label, opener, disk_path=None):
        self.label = label
        self.opener = opener
        self.disk_path = disk_path
        self.stem = stem_of(label)
        self.compression = compression_of(label)

    def open_text(self):
        raw = self.opener()
        if self.compression == "gz":
            raw = gzip.GzipFile(fileobj=raw)
        elif self.compression == "bz2":
            raw = bz2.BZ2File(raw)
        return io.TextIOWrapper(raw, encoding="ascii", errors="replace")


def list_sources(paths):
    sources = []
    for path in paths:
        if os.path.basename(path).startswith("."):
            continue
        if path.lower().endswith(".zip"):
            archive = zipfile.ZipFile(path)
            for member in archive.namelist():
                base = os.path.basename(member)
                if member.endswith("/") or not base or base.startswith(".") or "__MACOSX" in member:
                    continue
                if base.lower().split(".")[-1] in ("bam", "sam", "cram", "zip"):
                    sys.exit("ERROR: %s inside %s: unzip it first, BAM, SAM, CRAM and nested zip files"
                             " are only read directly" % (member, path))
                sources.append(Source(base, lambda a=archive, m=member: a.open(m)))
        else:
            sources.append(Source(os.path.basename(path), lambda p=path: open(p, "rb"), path))
    return sources


def detect_format(source):
    with source.open_text() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if re.match(r"^@(HD|SQ|RG|PG|CO)\t", line):
                return "sam"
            if line.startswith("@"):
                return "fastq"
            if line.startswith(">"):
                return "fasta"
            if line.startswith("LOCUS"):
                return "genbank"
            if line.startswith("ID   "):
                return "embl"
            return "unknown"
    return "empty"


def read_fastq(handle, label):
    while True:
        header = handle.readline()
        if not header:
            return
        if not header.strip():
            continue
        sequence = handle.readline().strip()
        plus = handle.readline()
        quality = handle.readline().strip()
        if not header.startswith("@") or not plus.startswith("+"):
            sys.exit("ERROR: %s is not a 4 line FASTQ file near read %s" % (label, header.strip()))
        yield header[1:].split()[0] if header[1:].split() else "read", sequence, quality


def read_fasta(handle):
    name, chunks = None, []
    for line in handle:
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(chunks), None
            name = line[1:].split()[0] if line[1:].split() else "seq"
            chunks = []
        elif name is not None:
            chunks.append(re.sub(r"[^A-Za-z]", "", line))
    if name is not None:
        yield name, "".join(chunks), None


def read_flatfile(handle, start, name_field):
    """GenBank (start ORIGIN, name from LOCUS) or EMBL (start SQ, name from ID)."""
    name, chunks, in_seq = None, [], False
    for line in handle:
        if line.startswith(name_field):
            fields = line[len(name_field):].replace(";", " ").split()
            name = fields[0] if fields else "seq"
        elif line.startswith(start):
            in_seq, chunks = True, []
        elif line.startswith("//"):
            if in_seq:
                yield name or "seq", "".join(chunks), None
            name, in_seq = None, False
        elif in_seq:
            chunks.append(re.sub(r"[^A-Za-z]", "", line))


def records(source, fmt):
    handle = source.open_text()
    if fmt == "fastq":
        return read_fastq(handle, source.label)
    if fmt == "fasta":
        return read_fasta(handle)
    if fmt == "genbank":
        return read_flatfile(handle, "ORIGIN", "LOCUS")
    if fmt == "embl":
        return read_flatfile(handle, "SQ ", "ID ")
    sys.exit("ERROR: %s is neither FASTQ, FASTA, GenBank nor EMBL (first line not recognised)" % source.label)


def clean(sequence):
    return NOT_ACGTN.sub("N", sequence.upper().translate(CLEAN))


def pieces(length, frag, step):
    if length <= frag:
        return [0]
    starts = list(range(0, length - frag + 1, step))
    if starts[-1] + frag < length:
        starts.append(length - frag)
    return starts


def sample(source, fmt, max_length):
    """Look at the first records: is a long sequence or a U present?"""
    long_seq, has_u = False, False
    for i, (_, sequence, _) in enumerate(records(source, fmt)):
        if len(sequence) > max_length:
            long_seq = True
        if "U" in sequence or "u" in sequence:
            has_u = True
        if i + 1 >= SAMPLE_RECORDS:
            break
    return long_seq, has_u


def link_or_copy(src, dst):
    try:
        os.link(os.path.realpath(src), dst)
    except OSError:
        shutil.copyfile(src, dst)


def write_fastq(handle, name, sequence, quality):
    handle.write("@%s\n%s\n+\n%s\n" % (name, sequence, quality))


def prepare(args):
    os.makedirs("align", exist_ok=True)
    os.makedirs("scan_only", exist_ok=True)
    sources = list_sources(args.inputs)
    if not sources:
        sys.exit("ERROR: no sequence file found in the input")

    # One decision per file, from a sample of its first records.
    plan = []
    for source in sources:
        fmt = detect_format(source)
        if fmt == "empty":
            sys.stderr.write("skipped, empty: %s\n" % source.label)
            continue
        if fmt == "sam":
            sys.exit("ERROR: %s is a SAM file: give it a .sam extension so that samtools converts it" % source.label)
        if fmt == "unknown":
            records(source, fmt)
        long_seq, has_u = sample(source, fmt, args.max_length)
        shred = args.kind == "assembly" or long_seq
        passthrough = (fmt == "fastq" and not shred and not has_u
                       and source.disk_path is not None and source.compression in ("", "gz"))
        plan.append({"source": source, "format": fmt, "shred": shred, "passthrough": passthrough})

    # Pair mates by stem: x_1 with x_2. Shredded files stay unpaired.
    by_stem = {}
    for item in plan:
        match = MATE.match(item["source"].stem)
        if match and not item["shred"]:
            by_stem.setdefault(match.group(1), {})[match.group(2)] = item
    for mates in by_stem.values():
        if "1" in mates and "2" in mates:
            mates["1"]["role"], mates["2"]["role"] = "R1", "R2"

    # Assembly pieces overlap by half, reads pieces do not (see the docstring).
    step = args.frag_length // 2 if args.kind == "assembly" else args.frag_length
    report = []
    seen = set()
    for item in plan:
        source = item["source"]
        role = item.get("role", "U")
        stem = source.stem
        while stem in seen:
            stem += "_"
        seen.add(stem)
        if item["passthrough"]:
            suffix = ".fastq.gz" if source.compression == "gz" else ".fastq"
            link_or_copy(source.disk_path, os.path.join("align", "%s__%s%s" % (role, stem, suffix)))
            report.append((source.label, item["format"], role, "used as is", "NA", "NA", "NA"))
            continue
        n_seq = n_bases = longest = n_pieces = 0
        if item["shred"]:
            full = gzip.open(os.path.join("scan_only", "F__%s.fastq.gz" % stem), "wt", compresslevel=1)
            out = gzip.open(os.path.join("align", "S__%s.fastq.gz" % stem), "wt", compresslevel=1)
        else:
            full = None
            out = gzip.open(os.path.join("align", "%s__%s.fastq.gz" % (role, stem)), "wt", compresslevel=1)
        for name, sequence, quality in records(source, item["format"]):
            sequence = clean(sequence)
            if not sequence:
                continue
            if quality is None or len(quality) != len(sequence):
                quality = QUAL * len(sequence)
            n_seq += 1
            n_bases += len(sequence)
            longest = max(longest, len(sequence))
            if full is None:
                write_fastq(out, name, sequence, quality)
                continue
            write_fastq(full, name, sequence, quality)
            for start in pieces(len(sequence), args.frag_length, step):
                end = start + args.frag_length
                write_fastq(out, "%s_%d" % (name, start + 1), sequence[start:end], quality[start:end])
                n_pieces += 1
        out.close()
        if full is not None:
            full.close()
        action = "cut into %d pieces of %d bp" % (n_pieces, args.frag_length) if item["shred"] else "converted"
        report.append((source.label, item["format"], "S" if item["shred"] else role, action,
                       n_seq, n_bases, longest))
        if args.kind == "reads" and longest > 100000:
            sys.stderr.write("WARNING: %s holds sequences of %d bp. If it is an assembly, set seq_kind"
                             " to \"assembly\" so that depths are interpreted as such.\n" % (source.label, longest))

    with open(args.report, "w") as handle:
        handle.write("file\tformat\trole\taction\tsequences\tbases\tlongest\n")
        for row in report:
            handle.write("\t".join(str(x) for x in row) + "\n")
    if not os.listdir("align"):
        sys.exit("ERROR: no sequence left after reading the input")


def to_fasta(args):
    """Reference mode: FASTA, GenBank or EMBL in, one FASTA per input file out."""
    os.makedirs(args.out_dir, exist_ok=True)
    for source in list_sources(args.inputs):
        fmt = detect_format(source)
        if fmt == "fastq" or fmt == "sam":
            sys.exit("ERROR: %s is a %s file, a reference must be FASTA, GenBank or EMBL" % (source.label, fmt.upper()))
        if fmt == "empty":
            continue
        with open(os.path.join(args.out_dir, "%s.fasta" % source.stem), "w") as out:
            if fmt == "fasta":
                # Copied as is: Split_fasta does its own header cleaning.
                with source.open_text() as handle:
                    shutil.copyfileobj(handle, out)
                continue
            for name, sequence, _ in records(source, fmt):
                out.write(">%s\n%s\n" % (name, sequence.upper()))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("inputs", nargs="+", help="sequence files, compressed or not, or .zip archives")
    parser.add_argument("--kind", choices=("reads", "assembly"), default="reads")
    parser.add_argument("--max-length", type=int, default=1000,
                        help="longer sequences are cut for the alignments (default 1000)")
    parser.add_argument("--frag-length", type=int, default=250,
                        help="length of the pieces (default 250)")
    parser.add_argument("--report", default="prepare_reads.tsv")
    parser.add_argument("--to-fasta", action="store_true", help="reference mode, see the docstring")
    parser.add_argument("--out-dir", default="fasta", help="output folder of --to-fasta")
    args = parser.parse_args()
    if args.to_fasta:
        to_fasta(args)
    else:
        prepare(args)


if __name__ == "__main__":
    main()
