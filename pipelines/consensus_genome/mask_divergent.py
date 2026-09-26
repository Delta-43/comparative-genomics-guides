#!/usr/bin/env python3
"""Mask reference bases that differ from aligned sequences in a multiple alignment (MAF).

Given a MAF whose first sequence in every block is the reference (e.g. `multiz` output built
from UCSC hg38 pairwise alignments) and the reference FASTA, write a FASTA of *identical length
and coordinates* in which every reference base that is not supported by all aligned species is
replaced with N:

  * mismatch  - any aligned species has a different base (or N) at that column: mask the base,
                plus --mismatch-flank bases on each side (default 0);
  * indel     - any aligned species has a gap at a reference base (deletion in that species), or
                a base where the reference has a gap (insertion): mask the affected reference
                bases plus --indel-flank bases on each side (default 6).

Regions outside any alignment block are left unchanged. Because lengths are preserved, gene
annotations (GTF) for the reference still apply, and reads from every aligned species map to the
masked genome without favouring the reference species.

Usage:
  mask_divergent.py --maf aln.maf --fasta hg38.fa --ref-prefix hg38. --out consensus.fa \
                    [--mismatch-flank 0] [--indel-flank 6] [--bed masked.bed]
Memory: one chromosome at a time (~250 MB for human chr1). Only the Python standard library.
"""
import argparse, collections, gzip, sys

def opener(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode.replace("t", "") if "b" in mode else mode)

def index_maf(path, ref_prefix):
    """chrom -> list of byte offsets of 'a' lines whose reference row is on that chrom."""
    idx = collections.defaultdict(list)
    with open(path, "rb") as fh:
        off = fh.tell(); line = fh.readline(); pending = None
        while line:
            if line.startswith(b"a"):
                pending = off
            elif line.startswith(b"s") and pending is not None:
                src = line.split()[1].decode()
                if not src.startswith(ref_prefix):
                    sys.exit(f"first sequence of block at byte {pending} is {src}, not the reference ({ref_prefix}*)")
                idx[src[len(ref_prefix):]].append(pending)
                pending = None
            off = fh.tell(); line = fh.readline()
    return idx

def read_block(fh, offset):
    fh.seek(offset); fh.readline()          # skip 'a' line
    rows = []
    for line in iter(fh.readline, b""):
        if not line.strip():
            break
        if line.startswith(b"s"):
            f = line.split()
            rows.append((f[1].decode(), int(f[2]), int(f[3]), f[4].decode(), f[6].upper()))
    return rows

def mask_chrom(seq, blocks, fh, mm_flank, indel_flank):
    """Returns the number of newly masked bases; edits `seq` (bytearray) in place."""
    n = len(seq); marks = bytearray(n)       # 1 = mask this position
    def mark(a, b):                          # half-open [a, b)
        marks[max(a, 0):min(b, n)] = b"\x01" * (min(b, n) - max(a, 0))
    for off in blocks:
        rows = read_block(fh, off)
        src, start, size, strand, ref = rows[0]
        if strand != "+":
            sys.exit(f"reference row on '-' strand at {src}:{start}; expected '+'")
        queries = [r[4] for r in rows[1:]]
        pos = start - 1                      # reference coordinate of the last reference base seen
        for col in range(len(ref)):
            r = ref[col]
            if r != 45:                      # 45 = '-'
                pos += 1
                for q in queries:
                    c = q[col]
                    if c == 45:              # deletion in query
                        mark(pos - indel_flank, pos + 1 + indel_flank); break
                    if c != r:               # mismatch (incl. N in query)
                        mark(pos - mm_flank, pos + 1 + mm_flank)
            else:                            # reference gap: insertion in a query, between pos and pos+1
                if any(q[col] != 45 for q in queries):
                    mark(pos + 1 - indel_flank, pos + 1 + indel_flank)
        if pos != start + size - 1:
            sys.exit(f"block at {src}:{start} has {pos - start + 1} reference bases, MAF says {size}")
    newly = 0
    for i in range(n):
        if marks[i] and seq[i] not in (78, 110):   # N / n
            seq[i] = 78; newly += 1
    return newly, marks

def fasta_records(path):
    name, chunks = None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name: yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name: yield name, "".join(chunks)

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--maf", required=True); ap.add_argument("--fasta", required=True)
    ap.add_argument("--ref-prefix", required=True, help="MAF source prefix of the reference, e.g. 'hg38.'")
    ap.add_argument("--out", required=True)
    ap.add_argument("--mismatch-flank", type=int, default=0)
    ap.add_argument("--indel-flank", type=int, default=6)
    ap.add_argument("--bed", help="optional BED of masked intervals")
    a = ap.parse_args()
    idx = index_maf(a.maf, a.ref_prefix)
    out, bed = open(a.out, "w"), (open(a.bed, "w") if a.bed else None)
    total = 0
    with open(a.maf, "rb") as fh:
        for name, s in fasta_records(a.fasta):
            seq = bytearray(s, "ascii")
            if name in idx:
                newly, marks = mask_chrom(seq, idx[name], fh, a.mismatch_flank, a.indel_flank)
                total += newly
                if bed:
                    i = 0
                    while i < len(marks):
                        if marks[i]:
                            j = i
                            while j < len(marks) and marks[j]: j += 1
                            bed.write(f"{name}\t{i}\t{j}\n"); i = j
                        else: i += 1
                print(f"{name}\t{len(idx[name])} blocks\t{newly} bases masked", file=sys.stderr)
            out.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60].decode() + "\n")
    print(f"total newly masked bases: {total}", file=sys.stderr)

if __name__ == "__main__":
    main()
