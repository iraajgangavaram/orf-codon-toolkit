from collections import Counter

import pytest

from orfkit import (
    codon_counts,
    find_orfs,
    gc_content,
    gc_skew,
    parse_fasta,
    reverse_complement,
    rscu,
    translate,
)

ORF_SEQ = "ATG" + "AAA" * 3 + "TAA"  # M K K K stop, 15 nt


def test_translate_and_reverse_complement():
    assert translate("ATGAAATAA") == "MK*"
    assert translate("ATGNNNTAA") == "MX*"
    assert reverse_complement("AACG") == "CGTT"


def test_gc_content_and_skew():
    assert gc_content("GGCCAATT") == 0.5
    assert gc_content("") == 0.0
    assert gc_skew("GGGGCC", window=6) == [(1, pytest.approx(1 / 3))]
    assert gc_skew("AAAA", window=4) == [(1, 0.0)]  # no G/C -> 0, not a crash
    with pytest.raises(ValueError):
        gc_skew("ACGT", window=0)


def test_forward_orf_coordinates():
    (orf,) = find_orfs("GG" + ORF_SEQ, min_aa=4, both_strands=False)
    assert (orf.strand, orf.frame, orf.start, orf.end) == ("+", 3, 3, 17)
    assert orf.protein == "MKKK" and orf.length_aa == 4


def test_frame_numbering():
    (orf,) = find_orfs("C" + ORF_SEQ, min_aa=4, both_strands=False)
    assert orf.frame == 2 and orf.start == 2


def test_reverse_strand_orf_maps_to_forward_coordinates():
    seq = reverse_complement("GG" + ORF_SEQ)  # ORF on minus strand, 2 nt of padding
    orfs = [o for o in find_orfs(seq, min_aa=4) if o.strand == "-"]
    assert len(orfs) == 1
    assert (orfs[0].start, orfs[0].end) == (1, 15)
    assert orfs[0].protein == "MKKK"


def test_nested_atg_reports_single_orf_from_first_start():
    (orf,) = find_orfs("ATGATGAAATAA", min_aa=1, both_strands=False)
    assert orf.protein == "MMK"


def test_unterminated_orf_ignored_and_min_length_respected():
    assert find_orfs("ATG" + "AAA" * 40, min_aa=1, both_strands=False) == []
    assert find_orfs(ORF_SEQ, min_aa=5, both_strands=False) == []


def test_orfs_sorted_longest_first():
    short = "ATGAAATAA"
    long_ = "ATG" + "GCT" * 10 + "TGA"
    orfs = find_orfs(short + "CC" + long_, min_aa=1, both_strands=False)
    assert [o.length_aa for o in orfs] == sorted((o.length_aa for o in orfs), reverse=True)
    assert orfs[0].length_aa == 11


def test_codon_counts_ignores_ambiguous_and_trailing_bases():
    assert codon_counts("ATGNNNAAAA") == Counter({"ATG": 1, "AAA": 1})


def test_rscu_values():
    values = rscu(Counter({"CTG": 3}))  # only one of six Leu codons used
    assert values["CTG"] == pytest.approx(6.0)
    assert values["TTA"] == 0.0
    values = rscu(Counter({"GCT": 1, "GCC": 1}))  # two of four Ala codons, equal
    assert values["GCT"] == values["GCC"] == pytest.approx(2.0)
    assert "ATG" not in values  # amino acids never observed are omitted


def test_parse_fasta_multiline_and_errors():
    lines = [">a desc", "atg", "AAA", "", ">b", "CCC"]
    assert list(parse_fasta(lines)) == [("a desc", "ATGAAA"), ("b", "CCC")]
    with pytest.raises(ValueError):
        list(parse_fasta(["ATG"]))
