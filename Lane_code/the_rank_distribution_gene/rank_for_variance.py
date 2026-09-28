# -*- coding: utf-8 -*-
"""Per-gene rank r that explains a given fraction (default 99%) of variance.

The rank r of the Lowrank Wu pipelines truncates K_g = Z_g Z_g' per gene, where
Z_g is the gene's block of the column-standardized additive design (the same
additive_design the estimator uses).  The eigenvalues of K_g are the squared
singular values s_i^2 of Z_g, so the variance explained by the top k is

    sum_{i<=k} s_i^2 / sum_i s_i^2        (denominator = tr(K_g) = n m_g),

and r_g is the smallest k where this reaches the threshold.

Genes can be given two ways (combinable):
    python rank_for_variance.py gene1.csv gene2.csv ...
        each CSV is one gene: n rows (individuals) by m_g columns (SNP dosages)
    python rank_for_variance.py --genotype geno.csv --G 10
        one n-by-m genotype split into G contiguous genes, exactly as the
        pipeline's split_into_genes does

Prints one r per line (one per gene, in gene order); --out writes the same.
"""
import argparse

import numpy as np
import pandas as pd


# Copied from the pipelines' Function_MCREML.py so this pre-step is standalone;
# keep them identical so the genes and their scaling match the main step.
def additive_design(real_data, stability_std=1e-12):
    """Column-standardize to mean 0, variance 1 (ddof=0).

    A near-constant column (std < stability_std) is left un-scaled.
    """
    M = np.asarray(real_data, dtype=float)
    mu = M.mean(axis=0)
    sd = M.std(axis=0)
    sd = np.where(sd < stability_std, 1.0, sd)
    return (M - mu) / sd


def split_into_genes(Z, G):
    """Split the m SNP columns of Z into G contiguous gene blocks.

    np.array_split makes the first (m mod G) genes one SNP larger.
    """
    m = Z.shape[1]
    return [Z[:, cols] for cols in np.array_split(np.arange(m), G)]


def rank_for_variance(Zg, threshold=0.99):
    """Smallest k with (sum of top-k eigenvalues of Z_g Z_g') / total >= threshold.

    Zg is a standardized n-by-m_g block.  Returns (r, cumulative fractions).
    """
    s = np.linalg.svd(np.asarray(Zg, dtype=float), compute_uv=False)
    lam = s ** 2                                 # eigenvalues of K_g, descending
    total = lam.sum()
    if total <= 0:                               # all columns monomorphic
        return 0, np.zeros_like(lam)
    frac = np.cumsum(lam) / total
    # small tolerance so a threshold hit exactly (e.g. 1.0) is not missed to rounding
    r = int(np.searchsorted(frac, threshold - 1e-12) + 1)
    return min(r, len(lam)), frac


def ranks_for_genes(gene_blocks, threshold=0.99):
    """One r per standardized gene block."""
    return [rank_for_variance(Zg, threshold)[0] for Zg in gene_blocks]


def _load(path):
    return pd.read_csv(path, header=None).to_numpy()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genes", nargs="*", help="one CSV per gene (n x m_g dosages)")
    ap.add_argument("--genotype", help="n x m genotype CSV to split into --G genes")
    ap.add_argument("--G", type=int, help="number of contiguous genes for --genotype")
    ap.add_argument("--threshold", type=float, default=0.99)
    ap.add_argument("--out", help="optional file to write the r values, one per line")
    args = ap.parse_args()

    blocks = []
    for path in args.genes:
        # each gene file is standardized on its own, as a gene block would be
        blocks.append(additive_design(_load(path)))
    if args.genotype:
        if not args.G:
            ap.error("--genotype needs --G")
        # split the RAW dosages, then standardize per block: standardization is
        # per column so this equals splitting additive_design(X), but never holds
        # a second full n-by-m float copy (n=32000, m=16000 is ~4 GB each)
        X = _load(args.genotype)
        for Xg in split_into_genes(X, args.G):
            blocks.append(additive_design(Xg))
    if not blocks:
        ap.error("give gene CSVs and/or --genotype with --G")

    # one number per gene, one per line, in gene order
    ranks = ranks_for_genes(blocks, args.threshold)
    text = "\n".join(str(r) for r in ranks) + "\n"
    print(text, end="")
    if args.out:
        with open(args.out, "w") as f:
            f.write(text)


if __name__ == "__main__":
    main()
