#!/usr/bin/env python3
"""
Count panel k-mers in raw reads, alignment-free.

Reads FASTQ or FASTA (optionally gzipped), counts how many times each panel
k-mer occurs, and writes per-site reference and alternative allele counts.
Reads are canonicalised, so strand does not matter.

Throughput note
---------------
The k-mers of the reads are turned into 2-bit integer codes with numpy and
looked up among the codes of the panel k-mers, by --threads worker processes
that each take a block of reads: roughly 8 Mbp of reads per second and per
thread, so that a deep whole-genome sample (2.6 Gbp, 220x on yeast) is scanned
in about a minute on 12 threads. Without numpy, or for k above 32, the k-mers
are counted one by one in pure Python, in one process, at 1-3 Mbp per second.
The counts, hence the output files, are the same either way. Two ways to spend
less time:

  * --max-bases  Subsample to a number of bases. Coverage, not read count, is
    what decides whether a site gets called, and read length varies by an order
    of magnitude between platforms: 400000 reads is 60 Mbp at 150 bp and 14 Mbp
    at 35 bp, which is 5x and 1.2x on a 12 Mb genome. Capping by reads silently
    produces an unusable scan on short-read data, so the cap is in bases.
    --max-reads still exists and still works, but prefer --max-bases.

  * --jellyfish  Hand the counting to jellyfish when it is installed. The panel
    k-mers are then queried against its database instead of being counted here.

A site is called only when its total depth reaches --min-depth; below that the
site is reported as not covered rather than as reference, because absence of a
k-mer is also what low coverage looks like.
"""

import argparse
import collections
import csv
import gc
import gzip
import mmap
import multiprocessing
import os
import subprocess
import sys

try:
    import numpy as np
except ImportError:  # the k-mers are then counted in pure Python
    np = None

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def canonical(kmer):
    reverse = revcomp(kmer)
    return kmer if kmer <= reverse else reverse


def opener(path):
    if path == "-":
        return sys.stdin
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def iter_sequences(paths, max_reads, max_bases=0):
    """Yield read sequences from FASTQ or FASTA inputs, stopping at either cap."""
    emitted = 0
    emitted_bases = 0
    for path in paths:
        with opener(path) as handle:
            first = handle.readline()
            if not first:
                continue
            is_fastq = first.startswith("@")
            handle.seek(0) if path != "-" else None
            if is_fastq:
                while True:
                    header = handle.readline()
                    if not header:
                        break
                    sequence = handle.readline().strip().upper()
                    handle.readline()
                    handle.readline()
                    if not sequence:
                        break
                    yield sequence
                    emitted += 1
                    emitted_bases += len(sequence)
                    if max_reads and emitted >= max_reads:
                        return
                    if max_bases and emitted_bases >= max_bases:
                        return
            else:
                chunks = []
                for line in handle:
                    line = line.strip()
                    if line.startswith(">"):
                        if chunks:
                            sequence = "".join(chunks).upper()
                            yield sequence
                            emitted += 1
                            emitted_bases += len(sequence)
                            if max_reads and emitted >= max_reads:
                                return
                            if max_bases and emitted_bases >= max_bases:
                                return
                        chunks = []
                    elif line:
                        chunks.append(line)
                if chunks:  # the last sequence of the file counts towards the caps too
                    sequence = "".join(chunks).upper()
                    yield sequence
                    emitted += 1
                    emitted_bases += len(sequence)
                    if max_reads and emitted >= max_reads:
                        return
                    if max_bases and emitted_bases >= max_bases:
                        return


def count_python(paths, wanted, k, max_reads, progress_every, max_bases=0):
    counts = dict.fromkeys(wanted, 0)
    bases = 0
    for index, sequence in enumerate(iter_sequences(paths, max_reads, max_bases), start=1):
        bases += len(sequence)
        # Reads with ambiguous bases are split: an N can never be part of a
        # panel k-mer, and skipping the whole read would bias against sites in
        # lower-quality regions.
        for piece in sequence.split("N"):
            for i in range(len(piece) - k + 1):
                key = canonical(piece[i:i + k])
                if key in counts:
                    counts[key] += 1
        if progress_every and index % progress_every == 0:
            sys.stderr.write("  %d reads, %.1f Mbp\n" % (index, bases / 1e6))
            sys.stderr.flush()
    return counts, bases


def count_jellyfish(paths, wanted, k, threads, hash_size, workdir):
    """Count with jellyfish, then query only the panel k-mers."""
    database = os.path.join(workdir, "reads.jf")
    command = ["jellyfish", "count", "-C", "-m", str(k), "-s", hash_size,
               "-t", str(threads), "-o", database] + list(paths)
    sys.stderr.write("running: %s\n" % " ".join(command))
    subprocess.run(command, check=True)

    query_input = os.path.join(workdir, "panel.kmers")
    with open(query_input, "w") as handle:
        for kmer in wanted:
            handle.write(kmer + "\n")
    result = subprocess.run(["jellyfish", "query", "-s", query_input, database],
                            check=True, capture_output=True, text=True)
    counts = dict.fromkeys(wanted, 0)
    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) == 2:
            counts[canonical(parts[0])] = int(parts[1])
    return counts, 0


# numpy counter. Bases get 2-bit codes (A, C, G, T = 0, 1, 2, 3) and a k-mer
# the code of its bases, first base in the high bits, so that codes sort as
# the k-mers do. A read k-mer is canonicalised in count_python, which is the
# same as looking it up as it is among the panel k-mers and their reverse
# complements. The codes of the first 16 bases of every read window are
# hashed into a table that marks those of the panel k-mers: the few windows
# that pass get the code of all their k bases, which is looked up among the
# sorted panel codes. A window holding anything but A, C, G or T, of either
# case, is skipped, as the N that splits a read in count_python. Reads are
# joined with a newline, which no window may hold either.
READ_BLOCK = 1 << 22  # bytes of FASTQ records, or of FASTA sequences, sent to a worker at a time
SLOT_BYTES = 2 * READ_BLOCK  # a block is handed to a worker in a slot of shared memory, if it fits
BATCH_BASES = 1 << 21  # bases encoded at a time in a worker, which bounds its memory
HEAD = 16  # bases of a window hashed into the table
TABLE_BITS = 26  # the table holds 2^26 booleans (64 MB)
HASH_MULT = 0x9E3779B1
CODE_WIDTH = {1: "uint8", 2: "uint8", 4: "uint16", 8: "uint32"}  # holds the code of 2 x this many bases

# Set before the workers are forked, which share them: (k, table, sorted codes,
# panel index of each code), and the slots of shared memory, which carry the
# blocks 7 times faster than the pipe of the pool.
_lookup = None
_slots = None


def base_values(letters):
    """Base code of every byte value: 0 to 3 for the letters of A, C, G, T, 4 for anything else."""
    values = np.full(256, 4, dtype=np.uint8)
    for code, bases in enumerate(letters):
        for letter in bases:
            values[ord(letter)] = code
    return values


if np is not None:
    READ_VALUES = base_values(("Aa", "Cc", "Gg", "Tt"))  # reads are upper-cased in count_python
    PANEL_VALUES = base_values(("A", "C", "G", "T"))


def panel_codes(kmers, k):
    """Codes a read k-mer is looked up among, sorted, and the index in kmers of each.

    The panel k-mers a read can match, each with its reverse complement: as
    count_python finds a k-mer only once canonicalised and upper-cased, a panel
    k-mer that is not canonical, or that holds anything but A, C, G and T in
    upper case, never gets a count.
    """
    bases = PANEL_VALUES[np.frombuffer("".join(kmers).encode("ascii", "replace"), dtype=np.uint8)]
    bases = bases.reshape(len(kmers), k)
    usable = (bases < 4).all(axis=1)
    forward = np.zeros(len(kmers), dtype=np.uint64)
    reverse = np.zeros(len(kmers), dtype=np.uint64)
    for position in range(k):
        column = bases[:, position].astype(np.uint64)
        forward = (forward << np.uint64(2)) | column
        reverse |= (np.uint64(3) - column) << np.uint64(2 * position)
    owner = np.flatnonzero(usable & (forward <= reverse))
    codes = np.concatenate((forward[owner], reverse[owner]))
    owner = np.concatenate((owner, owner))
    order = np.argsort(codes)
    codes, owner = codes[order], owner[order]
    single = np.ones(len(codes), dtype=bool)  # a palindrome is its own reverse complement
    single[1:] = codes[1:] != codes[:-1]
    return codes[single], owner[single]


def hash_slots(codes):
    """Slot in the table of each code of HEAD bases (uint32)."""
    return (codes * np.uint32(HASH_MULT)) >> np.uint32(32 - TABLE_BITS)


def head_codes(bases, length):
    """Codes (uint32) of the windows of length bases (16 at most) of bases (0 to 3).

    The codes of 2, 4, 8 and 16 bases are built from those of half the length,
    then those of length bases from the powers of 2 that add up to it.
    """
    parts = [power for power in (16, 8, 4, 2, 1) if length & power]
    codes = {1: bases}
    power = 1
    while power < parts[0]:
        width = np.dtype(CODE_WIDTH[power]).type
        half = codes[power].astype(width, copy=False)
        codes[2 * power] = (half[:-power] << width(2 * power)) | half[power:]
        power *= 2
    size = len(bases) - length + 1
    head = codes[parts[0]][:size].astype(np.uint32)
    offset = parts[0]
    for power in parts[1:]:
        head = (head << np.uint32(2 * power)) | codes[power][offset:offset + size]
        offset += power
    return head


def panel_hits(data):
    """Index in the panel k-mers of each one found in data (reads joined by newlines)."""
    k, table, codes, owner = _lookup
    values = READ_VALUES[np.frombuffer(data, dtype=np.uint8)]
    found = [np.zeros(0, dtype=np.int64)]
    if len(codes):
        for start in range(0, len(values) - k + 1, BATCH_BASES):
            found.append(batch_hits(values[start:start + BATCH_BASES + k - 1], k, table, codes, owner))
    return np.concatenate(found)


def batch_hits(values, k, table, codes, owner):
    """panel_hits of the windows of a batch of read base codes (4: not a base)."""
    head = min(k, HEAD)
    bases = values & np.uint8(3)
    size = len(values) - k + 1
    prefix = head_codes(bases, head)[:size]
    start = np.flatnonzero(table[hash_slots(prefix)])  # windows whose first bases may be a panel k-mer's
    windows = prefix[start].astype(np.uint64)
    for position in range(head, k):
        windows = (windows << np.uint64(2)) | bases[start + position]
    missing = np.zeros(len(values) + 1, dtype=np.int32)  # running count of what is not a base
    np.cumsum(values > 3, dtype=np.int32, out=missing[1:])
    windows = windows[missing[start + k] == missing[start]]
    windows.sort()  # searchsorted is much faster on sorted keys
    index = np.searchsorted(codes, windows)
    np.minimum(index, len(codes) - 1, out=index)
    return owner[index[codes[index] == windows]]


def tally(sequences, first_empty):
    """What a worker sends back: the panel indices found and how many times, the read lengths."""
    lengths = np.fromiter(map(len, sequences), dtype=np.int64, count=len(sequences))
    index, times = np.unique(panel_hits(b"\n".join(sequences)), return_counts=True)
    return index.astype(np.int32), times.astype(np.int32), lengths, first_empty


def fastq_sequences(block, count=None):
    """Sequences of a block of FASTQ records, read as iter_sequences does.

    Four lines a record, the second one stripped. An empty sequence ends the
    file there: the sequences before it are returned, with its index (else -1).
    """
    lines = block.split(b"\n")
    if block.endswith(b"\n"):
        lines.pop()
    sequences = [line.strip() for line in lines[1:None if count is None else 4 * count:4]]
    first_empty = -1
    if not all(sequences):
        first_empty = sequences.index(b"")
        del sequences[first_empty:]
    return sequences, first_empty


def block_of(slot, data):
    """A block handed to a worker: data itself, or its length when it is in a slot."""
    if slot is None:
        return data
    return _slots[slot * SLOT_BYTES:slot * SLOT_BYTES + data]


def scan_fastq(slot, data):
    return tally(*fastq_sequences(block_of(slot, data)))


def scan_fasta(slot, data):
    return tally(block_of(slot, data).split(b"\n"), -1)


def record_end(buffer):
    """Length of the whole 4-line records at the start of buffer, 0 if there are none."""
    lines = buffer.count(b"\n")
    extra = lines % 4
    if lines == extra:
        return 0
    end = len(buffer)
    for _ in range(extra + 1):
        end = buffer.rfind(b"\n", 0, end)
    return end + 1


def open_binary(path):
    """Binary handle on a file, gunzipped by a gzip process when available (a core of its own)."""
    if not path.endswith(".gz"):
        return open(path, "rb"), None
    try:
        process = subprocess.Popen(["gzip", "-dc", path], stdout=subprocess.PIPE)
    except OSError:
        return gzip.open(path, "rb"), None
    return process.stdout, process


def read_blocks(paths, scan):
    """(file index, worker task, block) over the inputs, in order, until the scan stops.

    FASTQ: blocks of whole records. FASTA: blocks of sequences joined by
    newlines, read by iter_sequences itself.
    """
    for file_index, path in enumerate(paths):
        if scan.done:
            return
        handle, process = open_binary(path)
        complete = False  # read to its end
        try:
            first = handle.readline()
            fastq = first.startswith(b"@")
            complete = not first
            if fastq:
                buffer = first
                while not scan.done and file_index not in scan.ended:
                    chunk = handle.read(READ_BLOCK)
                    if not chunk:
                        complete = True
                        if buffer:
                            yield file_index, scan_fastq, buffer
                        break
                    buffer += chunk
                    end = record_end(buffer)
                    if end:
                        yield file_index, scan_fastq, buffer[:end]
                        buffer = buffer[end:]
        finally:
            if process is not None and not complete:
                process.kill()  # the rest of the file is not needed
            handle.close()
            if process is not None and process.wait() and complete:
                sys.exit("gzip -dc failed on %s" % path)
        if first and not fastq:
            batch, size = [], 0
            for sequence in iter_sequences([path], 0, 0):
                if scan.done:
                    return
                batch.append(sequence)
                size += len(sequence) + 1
                if size >= READ_BLOCK:
                    yield file_index, scan_fasta, "\n".join(batch).encode("ascii", "replace")
                    batch, size = [], 0
            if batch:
                yield file_index, scan_fasta, "\n".join(batch).encode("ascii", "replace")


class Scan:
    """Totals of the worker results, taken in file order, with the reads caps.

    The caps, the progress lines and the end of a file at an empty sequence
    follow iter_sequences and count_python read by read.
    """

    def __init__(self, size, max_reads, max_bases, progress_every):
        self.counts = np.zeros(size, dtype=np.int64)
        self.reads = 0
        self.bases = 0
        self.max_reads = max_reads
        self.max_bases = max_bases
        self.progress_every = progress_every
        self.done = False  # a cap is reached
        self.ended = set()  # files ended by an empty sequence

    def add(self, file_index, task, block, result):
        index, times, lengths, first_empty = result.get()  # also frees the slot of the block
        if self.done or file_index in self.ended:
            return
        stop = None
        if self.max_reads and self.reads + len(lengths) >= self.max_reads:
            stop = self.max_reads - self.reads
        if self.max_bases and self.bases + int(lengths.sum()) >= self.max_bases:
            reached = int(np.searchsorted(np.cumsum(lengths), self.max_bases - self.bases)) + 1
            stop = reached if stop is None else min(stop, reached)
        if stop is not None:
            self.done = True
            if stop < len(lengths):  # the reads of the block after the cap are not counted
                if task is scan_fastq:
                    sequences = fastq_sequences(block, stop)[0]
                else:
                    sequences = block.split(b"\n")[:stop]
                index, times = tally(sequences, -1)[:2]
                lengths = lengths[:stop]
        self.counts[index] += times
        if self.progress_every:
            cumulative = self.bases + np.cumsum(lengths)
            for read in range(self.reads // self.progress_every + 1,
                              (self.reads + len(lengths)) // self.progress_every + 1):
                read *= self.progress_every
                sys.stderr.write("  %d reads, %.1f Mbp\n" % (read, cumulative[read - self.reads - 1] / 1e6))
            sys.stderr.flush()
        self.reads += len(lengths)
        self.bases += int(lengths.sum())
        if first_empty >= 0:
            self.ended.add(file_index)


class Done:
    """Result of a task run in the main process, like a pool's AsyncResult."""

    def __init__(self, value):
        self.value = value

    def get(self):
        return self.value


def count_numpy(paths, wanted, k, max_reads, progress_every, max_bases, threads):
    """Same counts as count_python, with numpy and threads worker processes."""
    global _lookup, _slots
    kmers = list(wanted)  # all of length k
    codes, owner = panel_codes(kmers, k)
    table = np.zeros(1 << TABLE_BITS, dtype=bool)
    table[hash_slots((codes >> np.uint64(2 * (k - min(k, HEAD)))).astype(np.uint32))] = True
    _lookup = (k, table, codes, owner)
    scan = Scan(len(kmers), max_reads, max_bases, progress_every)
    pool = None
    free = []  # slots of shared memory not in use
    if threads > 1:
        # Forked workers share the parent's memory as long as nothing writes
        # to it: frozen, the panel objects are left alone by their garbage
        # collector. At most 2 x threads blocks are pending at a time.
        gc.freeze()
        free = list(range(2 * threads + 1))
        _slots = mmap.mmap(-1, len(free) * SLOT_BYTES)
        pool = multiprocessing.get_context("fork").Pool(threads)
    pending = collections.deque()
    try:
        for file_index, task, block in read_blocks(paths, scan):
            slot = None
            if pool is None:
                result = Done(task(None, block))
            elif len(block) <= SLOT_BYTES:
                slot = free.pop()
                _slots[slot * SLOT_BYTES:slot * SLOT_BYTES + len(block)] = block
                result = pool.apply_async(task, (slot, len(block)))
            else:  # a block larger than a slot (a long sequence) goes through the pipe
                result = pool.apply_async(task, (None, block))
            pending.append((file_index, task, block, result, slot))
            while len(pending) > 2 * threads:
                *taken, slot = pending.popleft()
                scan.add(*taken)
                if slot is not None:
                    free.append(slot)
        while pending:
            scan.add(*pending.popleft()[:4])
    finally:
        if pool:
            pool.terminate()
            pool.join()
    counts = dict.fromkeys(wanted, 0)
    for position in np.flatnonzero(scan.counts):
        counts[kmers[position]] = int(scan.counts[position])
    return counts, scan.bases


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--kmers", required=True, help="panel k-mer lookup TSV")
    parser.add_argument("--reads", required=True, nargs="+", help="FASTQ or FASTA inputs")
    parser.add_argument("--out", required=True, help="output per-site allele counts TSV")
    parser.add_argument("--min-depth", type=int, default=3,
                        help="total allele count needed to call a site (default 3)")
    parser.add_argument("--min-allele-fraction", type=float, default=0.8,
                        help="fraction of reads needed to call a homozygous allele")
    parser.add_argument("--max-bases", type=int, default=0,
                        help="stop after this many bases, 0 to use all (the default). "
                             "Prefer this over --max-reads: it is coverage that "
                             "decides whether a site is called, and bases are the "
                             "unit of coverage. For a 12 Mb genome, 300000000 is "
                             "about 25x, which is ample.")
    parser.add_argument("--max-reads", type=int, default=0,
                        help="stop after this many reads, 0 to use all (the default). "
                             "Kept for convenience, but a read cap means different "
                             "coverage on different platforms")
    parser.add_argument("--extra-kmers",
                        help="TSV kmer, class: further k-mers counted in the same pass, "
                             "reported per class (the mating-type k-mers)")
    parser.add_argument("--extra-out", help="output TSV of the --extra-kmers counts")
    parser.add_argument("--jellyfish", action="store_true",
                        help="count with jellyfish instead of in Python")
    parser.add_argument("--threads", type=int, default=4,
                        help="worker processes of the numpy counter, or jellyfish threads (default 4)")
    parser.add_argument("--hash-size", default="100M")
    parser.add_argument("--workdir", default=".")
    parser.add_argument("--progress-every", type=int, default=100000)
    args = parser.parse_args()

    # A site carries one k-mer per haplotype seen in the panel strains, each
    # labelled with the allele at its centre, so an allele is counted as the
    # sum over its k-mers.
    sites = {}
    wanted = set()
    with opener(args.kmers) as handle:
        rows = csv.reader(handle, delimiter="\t")  # 1.5 million rows: faster than a DictReader
        header = next(rows, None)
        if header:
            columns = [header.index(name) for name in ("kmer", "site_id", "allele")]
            for row in rows:
                if row:  # empty lines are skipped, as by a DictReader
                    kmer, site_id, allele = (row[column] for column in columns)
                    wanted.add(kmer)
                    sites.setdefault(site_id, {}).setdefault(allele, []).append(kmer)
    if not wanted:
        sys.exit("empty k-mer panel: %s" % args.kmers)
    extra = {}
    if args.extra_kmers:
        with opener(args.extra_kmers) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                extra.setdefault(row["class"], []).append(row["kmer"])
                wanted.add(row["kmer"])
    k = len(next(iter(wanted)))
    if any(len(kmer) != k for kmer in wanted):
        sys.exit("panel mixes k-mer lengths, rebuild it with one k")

    sys.stderr.write("panel: %d sites, %d k-mers, k=%d\n" % (len(sites), len(wanted), k))
    if args.jellyfish:
        counts, bases = count_jellyfish(args.reads, wanted, k, args.threads,
                                        args.hash_size, args.workdir)
    else:
        if np is not None and k <= 32 and "-" not in args.reads:
            counts, bases = count_numpy(args.reads, wanted, k, args.max_reads,
                                        args.progress_every, args.max_bases, args.threads)
        else:
            counts, bases = count_python(args.reads, wanted, k, args.max_reads,
                                         args.progress_every, args.max_bases)
        sys.stderr.write("scanned %.1f Mbp\n" % (bases / 1e6))

    rows = []
    for site_id, alleles in sites.items():
        ref_count = sum(counts.get(kmer, 0) for kmer in alleles.get("R", ()))
        alt_count = sum(counts.get(kmer, 0) for kmer in alleles.get("A", ()))
        depth = ref_count + alt_count
        if depth < args.min_depth:
            call, fraction = "?", 0.0
        else:
            fraction = alt_count / float(depth)
            if fraction >= args.min_allele_fraction:
                call = "A"
            elif fraction <= 1.0 - args.min_allele_fraction:
                call = "R"
            else:
                # Both alleles well represented: a heterozygous diploid, a mixed
                # culture, or a cross. Reported as its own state, never forced
                # into one of the two homozygous calls.
                call = "H"
        rows.append({"site_id": site_id, "ref_count": ref_count, "alt_count": alt_count,
                     "depth": depth, "alt_fraction": round(fraction, 4), "call": call})

    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["site_id", "ref_count", "alt_count", "depth",
                                "alt_fraction", "call"], delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: row["site_id"]))

    if extra and args.extra_out:
        with open(args.extra_out, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["class", "n_kmers", "kmers_seen", "total_count", "median_count"])
            for name in sorted(extra):
                values = sorted(counts.get(kmer, 0) for kmer in extra[name])
                writer.writerow([name, len(values), sum(1 for v in values if v),
                                 sum(values), values[len(values) // 2]])

    tally = {}
    for row in rows:
        tally[row["call"]] = tally.get(row["call"], 0) + 1
    sys.stderr.write("calls: %s\n" % ", ".join(
        "%s=%d" % (call, tally[call]) for call in sorted(tally)))
    if tally.get("?", 0) > len(rows) * 0.5:
        sys.stderr.write(
            "warning: over half the panel is below --min-depth after %.1f Mbp. "
            "Either the scan was capped too low (raise --max-bases, 0 scans "
            "everything), the run is too shallow, or k is too close to the read "
            "length, which leaves too few k-mers per read\n" % (bases / 1e6))
    return 0


if __name__ == "__main__":
    sys.exit(main())
