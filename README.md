# orf-codon-toolkit

A small, dependency-free Python toolkit for basic prokaryotic sequence analysis:
ORF finding on both strands, codon usage (RSCU), and GC content / GC skew.
Pure standard library, so it runs anywhere Python 3.9+ does.

## Install

```bash
pip install -e .
```

## Command line

```bash
python -m orfkit orfs  examples/example.fa --min-aa 30   # ORFs, both strands (TSV)
python -m orfkit codons cds.fa                           # RSCU per record
python -m orfkit skew  genome.fa --window 10000          # windowed GC skew
python -m orfkit gc    examples/example.fa               # GC content per record
```

Output is tab-separated and written to stdout.

## Python API

```python
from orfkit import find_orfs, codon_counts, rscu, gc_skew

orfs = find_orfs(seq, min_aa=30)          # longest first
usage = rscu(codon_counts(cds))           # relative synonymous codon usage
skew = gc_skew(genome, window=10_000)     # [(window_start, skew), ...]
```

## Method notes

- **ORFs** run from an `ATG` to the next in-frame stop codon (standard genetic
  code, table 1). Within a frame, each ORF starts at the first `ATG` after the
  previous stop, so nested ORFs are not reported separately. ORFs that run off
  the end of the sequence without a stop are ignored.
- **Coordinates** are 1-based, inclusive, on the forward strand, and include
  the stop codon. For `-` strand ORFs, `start <= end` still holds and `frame`
  is numbered on the reverse complement.
- **RSCU** is observed count divided by the mean count of the synonymous codons
  for that amino acid; amino acids never observed are omitted.
- **GC skew** is `(G - C) / (G + C)` per window (0 when a window has no G or C).
  A switch in the sign of cumulative skew is often used to locate replication
  origin and terminus in bacterial genomes.

## Limitations

- Standard genetic code only; no alternative start codons (`GTG`, `TTG`) and
  no alternative translation tables.
- This is a simple ORF scanner, not a gene predictor: short or spurious ORFs
  are common, especially on the opposite strand. Use a tool such as Prodigal
  for annotation-grade gene calling.

## Tests

```bash
pip install pytest
pytest
```
