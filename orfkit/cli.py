"""Command-line interface: python -m orfkit <command> FASTA [options]."""

from __future__ import annotations

import argparse
import csv
import sys

from .core import codon_counts, find_orfs, gc_content, gc_skew, parse_fasta, rscu


def _records(path: str):
    with open(path) as fh:
        yield from parse_fasta(fh)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="orfkit", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    o = sub.add_parser("orfs", help="find ORFs and write a TSV")
    o.add_argument("fasta")
    o.add_argument("--min-aa", type=int, default=30)
    o.add_argument("--forward-only", action="store_true")

    c = sub.add_parser("codons", help="RSCU table per record (assumes records are CDSs)")
    c.add_argument("fasta")

    s = sub.add_parser("skew", help="windowed GC skew per record")
    s.add_argument("fasta")
    s.add_argument("--window", type=int, default=1000)
    s.add_argument("--step", type=int, default=None)

    g = sub.add_parser("gc", help="GC content per record")
    g.add_argument("fasta")

    args = p.parse_args(argv)
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")

    if args.cmd == "orfs":
        w.writerow(["record", "strand", "frame", "start", "end", "length_aa", "protein"])
        for name, seq in _records(args.fasta):
            for orf in find_orfs(seq, args.min_aa, not args.forward_only):
                w.writerow([name.split()[0], orf.strand, orf.frame, orf.start,
                            orf.end, orf.length_aa, orf.protein])
    elif args.cmd == "codons":
        w.writerow(["record", "codon", "rscu"])
        for name, seq in _records(args.fasta):
            for codon, val in sorted(rscu(codon_counts(seq)).items()):
                w.writerow([name.split()[0], codon, f"{val:.3f}"])
    elif args.cmd == "skew":
        w.writerow(["record", "window_start", "gc_skew"])
        for name, seq in _records(args.fasta):
            for start, val in gc_skew(seq, args.window, args.step):
                w.writerow([name.split()[0], start, f"{val:.4f}"])
    elif args.cmd == "gc":
        w.writerow(["record", "length", "gc_content"])
        for name, seq in _records(args.fasta):
            w.writerow([name.split()[0], len(seq), f"{gc_content(seq):.4f}"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
