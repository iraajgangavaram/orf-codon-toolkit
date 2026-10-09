"""Validation of orfkit against simulated sequences with known truth.

A. ORF scanner vs analytic expectation on random sequence (null model).
B. Recovery of planted genes (both strands) in random background.
C. Location of a planted GC-skew switch (replication origin/terminus analogue).

Run: python scripts/validate_orfs.py   (about a minute; all runs seeded)
"""

import csv
import json
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from orfkit import find_orfs, gc_skew  # noqa: E402
from orfkit.core import CODON_TABLE, STOP_CODONS, reverse_complement  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "docs")
SENSE = [c for c in CODON_TABLE if c not in STOP_CODONS]


def random_seq(rng, n, gc):
    return "".join(rng.choices("ACGT", [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2], k=n))


def expected_orfs_per_kb(gc, min_aa):
    """Renewal-process expectation for the forward strand, 3 frames, per kb.

    a = P(ATG) and s = P(stop) per codon; a cycle is a wait for ATG (mean 1/a
    codons) then a run to the next stop (mean 1/s codons). An ORF reaches
    min_aa amino acids if the min_aa - 1 codons after ATG are all non-stop.
    """
    pa = pt = (1 - gc) / 2
    pg = gc / 2
    a = pa * pt * pg
    s = pt * (pa * pa + pa * pg + pg * pa)
    rate = 1.0 / (1.0 / a + 1.0 / s)  # ORFs per codon per frame
    return rate * (1 - s) ** (min_aa - 1) * 3 * 1000 / 3


def experiment_a(min_aa=30, length=300_000, reps=5):
    rows = []
    for gc in (0.30, 0.40, 0.50, 0.60, 0.70):
        obs = []
        for r in range(reps):
            rng = random.Random(1000 * r + int(gc * 100))
            seq = random_seq(rng, length, gc)
            obs.append(len(find_orfs(seq, min_aa, both_strands=False)) / (length / 1000))
        mean = sum(obs) / reps
        sd = (sum((x - mean) ** 2 for x in obs) / (reps - 1)) ** 0.5
        rows.append({"gc": gc, "observed_per_kb": round(mean, 4), "sd": round(sd, 4),
                     "expected_per_kb": round(expected_orfs_per_kb(gc, min_aa), 4)})
    return rows


def plant(rng, n_genes=200, gc=0.5):
    """Random background with planted genes on both strands; returns seq, truth."""
    parts, truth, pos = [], [], 0
    for _ in range(n_genes):
        gap = random_seq(rng, rng.randint(200, 600), gc)
        n_codons = rng.randint(50, 500)
        body = "".join(rng.choices(SENSE, k=n_codons - 1))
        gene = "ATG" + body + rng.choice(sorted(STOP_CODONS))
        strand = rng.choice("+-")
        piece = gene if strand == "+" else reverse_complement(gene)
        parts += [gap, piece]
        start = pos + len(gap) + 1
        truth.append((strand, start, start + len(gene) - 1))
        pos += len(gap) + len(gene)
    return "".join(parts), truth


def experiment_b(reps=5):
    rows = []
    for r in range(reps):
        rng = random.Random(7000 + r)
        seq, truth = plant(rng)
        found = find_orfs(seq, min_aa=40)
        index = {(o.strand, o.end if o.strand == "+" else o.start): o for o in found}
        both = stop_only = missed = 0
        for strand, s, e in truth:
            key = (strand, e if strand == "+" else s)  # stop-codon end of the gene
            o = index.get(key)
            if o is None:
                missed += 1
            elif (o.start, o.end) == (s, e):
                both += 1
            else:
                stop_only += 1
        rows.append({"run": r, "genes": len(truth), "exact": both,
                     "stop_correct_start_upstream": stop_only, "missed": missed,
                     "orfs_reported": len(found)})
    return rows


def experiment_c(length=1_000_000, window=5000):
    rows = []
    curves = {}
    for delta in (0.01, 0.02, 0.05, 0.10):
        errs = []
        for r in range(10):
            rng = random.Random(500 + r)
            half = length // 2
            pg1, pc1 = 0.25 + delta / 2, 0.25 - delta / 2  # first half G-rich
            seq = "".join(rng.choices("ACGT", [0.25, pc1, pg1, 0.25], k=half)
                          + rng.choices("ACGT", [0.25, pg1, pc1, 0.25], k=length - half))
            cum, total = [], 0.0
            for _, sk in gc_skew(seq, window):
                total += sk
                cum.append(total)
            peak = max(range(len(cum)), key=cum.__getitem__) * window + window // 2
            errs.append(abs(peak - half))
            if r == 0:
                curves[delta] = cum
        rows.append({"skew_delta": delta, "true_switch": length // 2,
                     "median_abs_error_bp": sorted(errs)[len(errs) // 2],
                     "max_abs_error_bp": max(errs)})
    return rows, curves


def write_csv(name, rows):
    with open(os.path.join(OUT, "results", name), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def figure(a, curves, window=5000):
    from _style import BLUE, ORANGE, INK_2, apply
    import matplotlib.pyplot as plt

    apply()
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    gc = [r["gc"] for r in a]
    ax[0].plot(gc, [r["expected_per_kb"] for r in a], "--", color=INK_2, marker="s", label="Expected (analytic)")
    ax[0].errorbar(gc, [r["observed_per_kb"] for r in a], yerr=[r["sd"] for r in a],
                   color=BLUE, marker="o", capsize=3, label="Observed (mean +/- SD, 5 runs)")
    ax[0].set_ylim(0, max(r["expected_per_kb"] for r in a) * 1.25)
    ax[0].set_xlabel("GC content of random sequence")
    ax[0].set_ylabel("Spurious ORFs >= 30 aa per kb (one strand)")
    ax[0].set_title("Chance ORFs follow the null model")
    ax[0].legend(loc="upper left")
    for (delta, cum), col, mk in zip(curves.items(), (BLUE, ORANGE, "#1baf7a", INK_2), "osv^"):
        x = [(i + 0.5) * window / 1e6 for i in range(len(cum))]
        ax[1].plot(x, cum, color=col, label=f"skew difference {delta:g}", marker=mk, markevery=40, markersize=5)
    ax[1].axvline(0.5, color=INK_2, lw=0.8, ls=":")
    ax[1].text(0.51, ax[1].get_ylim()[0] * 0.9, "planted switch", color=INK_2, fontsize=9)
    ax[1].set_xlabel("Genome position (Mb)")
    ax[1].set_ylabel("Cumulative GC skew")
    ax[1].set_title("Cumulative skew peaks at the planted switch")
    ax[1].set_ylim(top=ax[1].get_ylim()[1] * 1.3)
    ax[1].legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "figures", "orf_validation.png"), dpi=200)


if __name__ == "__main__":
    a = experiment_a()
    b = experiment_b()
    c, curves = experiment_c()
    write_csv("validation_null_orfs.csv", a)
    write_csv("validation_planted_genes.csv", b)
    write_csv("validation_gc_skew.csv", c)
    figure(a, curves)
    tot = {k: sum(r[k] for r in b) for k in ("genes", "exact", "stop_correct_start_upstream", "missed")}
    with open(os.path.join(OUT, "results", "validation_summary.json"), "w") as fh:
        json.dump({"planted_genes_total": tot}, fh, indent=2)
    for name, rows in (("A", a), ("B", b), ("C", c)):
        print(name)
        for r in rows:
            print("  ", r)
    print(tot)
