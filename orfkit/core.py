"""Core sequence-analysis functions: FASTA parsing, ORF finding, codon usage, GC skew."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Tuple

_BASES = "TCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE: Dict[str, str] = {
    a + b + c: _AMINO[i]
    for i, (a, b, c) in enumerate(
        (a, b, c) for a in _BASES for b in _BASES for c in _BASES
    )
}
STOP_CODONS = frozenset(c for c, aa in CODON_TABLE.items() if aa == "*")
_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def parse_fasta(lines: Iterable[str]) -> Iterator[Tuple[str, str]]:
    """Yield (header, sequence) pairs from an iterable of FASTA lines."""
    header, chunks = None, []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks).upper()
            header, chunks = line[1:].strip(), []
        elif header is None:
            raise ValueError("FASTA data must begin with a '>' header line")
        else:
            chunks.append(line)
    if header is not None:
        yield header, "".join(chunks).upper()


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def translate(seq: str) -> str:
    """Translate a DNA sequence; codons containing non-ACGT become 'X'."""
    seq = seq.upper()
    return "".join(
        CODON_TABLE.get(seq[i : i + 3], "X") for i in range(0, len(seq) - 2, 3)
    )


def gc_content(seq: str) -> float:
    """Fraction of G/C among A/C/G/T bases (0.0 for an empty sequence)."""
    seq = seq.upper()
    acgt = sum(seq.count(b) for b in "ACGT")
    return (seq.count("G") + seq.count("C")) / acgt if acgt else 0.0


def gc_skew(seq: str, window: int, step: int | None = None) -> List[Tuple[int, float]]:
    """Windowed GC skew, (G - C) / (G + C). Returns (window_start_1based, skew)."""
    if window <= 0:
        raise ValueError("window must be positive")
    step = step or window
    seq = seq.upper()
    out = []
    for start in range(0, max(len(seq) - window + 1, 0), step):
        w = seq[start : start + window]
        g, c = w.count("G"), w.count("C")
        out.append((start + 1, (g - c) / (g + c) if g + c else 0.0))
    return out


@dataclass(frozen=True)
class ORF:
    strand: str  # '+' or '-'
    frame: int  # 1-3 within the strand
    start: int  # 1-based, forward-strand coordinates, start <= end
    end: int
    length_aa: int
    protein: str


def _scan_strand(seq: str, strand: str, min_aa: int) -> Iterator[ORF]:
    n = len(seq)
    for frame in range(3):
        start_pos = None
        for i in range(frame, n - 2, 3):
            codon = seq[i : i + 3]
            if start_pos is None and codon == "ATG":
                start_pos = i
            elif start_pos is not None and codon in STOP_CODONS:
                protein = translate(seq[start_pos:i])
                if len(protein) >= min_aa:
                    s, e = start_pos, i + 3  # e is exclusive, includes stop codon
                    if strand == "-":
                        s, e = n - e, n - s
                    yield ORF(strand, frame + 1, s + 1, e, len(protein), protein)
                start_pos = None


def find_orfs(seq: str, min_aa: int = 30, both_strands: bool = True) -> List[ORF]:
    """Find stop-terminated ORFs beginning at ATG, longest-first.

    Within a frame, each ORF starts at the first ATG after the previous stop
    (so nested ORFs are not reported separately). ORFs running off the end of
    the sequence without a stop codon are ignored.
    """
    seq = seq.upper()
    orfs = list(_scan_strand(seq, "+", min_aa))
    if both_strands:
        orfs += list(_scan_strand(reverse_complement(seq), "-", min_aa))
    return sorted(orfs, key=lambda o: (-o.length_aa, o.start))


def codon_counts(cds: str) -> Counter:
    """Count in-frame ACGT codons of a coding sequence (trailing bases ignored)."""
    cds = cds.upper()
    return Counter(
        cds[i : i + 3]
        for i in range(0, len(cds) - 2, 3)
        if set(cds[i : i + 3]) <= set("ACGT")
    )


def rscu(counts: Counter) -> Dict[str, float]:
    """Relative synonymous codon usage: observed / mean count of synonymous codons.

    Amino acids with no observed codons are omitted. Stop codons are included.
    """
    by_aa: Dict[str, List[str]] = {}
    for codon, aa in CODON_TABLE.items():
        by_aa.setdefault(aa, []).append(codon)
    out = {}
    for codons in by_aa.values():
        total = sum(counts.get(c, 0) for c in codons)
        if total == 0:
            continue
        expected = total / len(codons)
        for c in codons:
            out[c] = counts.get(c, 0) / expected
    return out
