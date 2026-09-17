# -*- coding: utf-8 -*-
# Pooled within-gene AxA extension of the single-gene AI-REML pipeline.
# The pooled kernel W = (1/G) sum_g K_g averages the AxA kernels of G
# subsampled contiguous SNP regions ("genes"), i.e. the estimand of eq. (3) in
# gxg_reml_within_gene.pdf. Shared building blocks (build_W_batched, reml, the
# single-gene simulators) are reused from Function_AIREML.
import numpy as np
import pandas as pd
from scipy.linalg import cholesky

from Function_AIREML import (
    build_W_batched,
    simulate_remove_sampling_err,
    reml,
    AI_REML,
)


def build_gene_kernel(Zg, pair_batch_size=5000):
    """Per-gene (within-region) AxA kernel  K_g = (1/P_g) H_g H_g^T.

    Zg : (n, m_g) standardized genotype columns for a single gene/region.
    H_g holds the standardized within-region pairwise products Z_a o Z_b for
    every pair (a, b) of columns in Zg; P_g = C(m_g, 2). This is exactly the
    per-gene kernel K_g of eq. (3) in gxg_reml_within_gene.pdf, and is
    identical to build_W_batched applied to the region.
    """
    return build_W_batched(Zg, pair_batch_size=pair_batch_size)


def build_pooled_W(Z, G, region_size, rng=None, pair_batch_size=5000):
    """Pooled within-gene AxA kernel from G subsampled contiguous regions.

    Treats the genotype matrix as a chromosome and subsamples G contiguous
    SNP windows ("genes") of `region_size` markers each. Windows are placed at
    random start positions and MAY OVERLAP. Each window g contributes a
    per-gene kernel K_g (build_gene_kernel), and the pooled kernel is their
    average

        W = (1/G) sum_{g=1..G} K_g,          K_g = (1/P_g) H_g H_g^T,

    i.e. the estimand W of eq. (3) in gxg_reml_within_gene.pdf, with the G
    genes here defined by subsampled contiguous regions rather than real gene
    boundaries. Pooling over features (regions) keeps all ~n^2/2 individual
    pairs in every K_g while lowering the pooled kernel's effective-marker
    count by ~sqrt(G) (eq. 9).

    Parameters
    ----------
    Z            : (n, m) standardized genotype matrix.
    G            : number of contiguous regions (genes) to subsample.
    region_size  : number of contiguous markers per region (m_g).
    rng          : optional np.random.Generator / RandomState for reproducible
                   region placement; if None a fresh default_rng() is used.
    pair_batch_size : column batch size for the H_g H_g^T accumulation.

    Returns
    -------
    W      : (n, n) pooled within-gene AxA kernel.
    starts : (G,) array of the region start indices used.
    """
    n, m = Z.shape
    if region_size < 2:
        raise ValueError("region_size must be >= 2 to contain at least one pair")
    if region_size > m:
        raise ValueError(f"region_size ({region_size}) exceeds number of markers ({m})")
    if rng is None:
        rng = np.random.default_rng()

    # Random contiguous start positions in [0, m - region_size]; overlap allowed.
    starts = rng.integers(0, m - region_size + 1, size=G)

    W = np.zeros((n, n))
    for start in starts:
        Zg = Z[:, start:start + region_size]
        W += build_gene_kernel(Zg, pair_batch_size=pair_batch_size)
    W /= G

    return W, starts


def simulate_Cholesky_pooled_withadd(real_data, G, region_size, seed=None,
                                     s2a=0.2, s2gxg=0.7, s2e=0.1, stability=1e-10):
    """Pooled-kernel counterpart of simulate_Cholesky_from_std_withadd.

    Identical to the single-gene builder except the epistatic GRM is the
    pooled within-gene kernel W_pooled = (1/G) sum_g K_g over G subsampled
    contiguous regions (build_pooled_W, eq. 3) instead of the all-pairs W.

    The pooled kernel W is also returned so it can be stored once and reused
    directly by AI_REML_pooled -- build_pooled_W then runs exactly once and the
    simulated y and the fitted kernel share the same W by construction (no need
    to re-draw regions from a seed at fit time).

    Returns (Lgxg, La, W): the lower-Cholesky factors of s2gxg*W and s2a*K, and
    the pooled kernel W itself.
    """
    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0)
    n, m = Z.shape

    # Pooled epistatic GRM (built once)
    rng = np.random.default_rng(seed)
    W, _ = build_pooled_W(Z, G, region_size, rng=rng)
    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)

    K = (Z @ Z.T) / m
    La = cholesky(s2a * K + stability * np.eye(n), lower=True)

    return Lgxg, La, W


def AI_REML_pooled(Z, y, W, iters=12):
    """Exact AI-REML for the pooled within-gene AxA model.

    Same interface / return ordering as AI_REML, but the epistatic component
    is a *precomputed* pooled kernel W instead of the all-pairs kernel:

        V = s2a*K + s2gxg*W + s2e*I,
        K = Z Z' / m                                     (additive, full Z)
        W = (1/G) sum_g K_g   (built once by build_pooled_W, eq. 3, and stored)

    W is passed in rather than rebuilt so build_pooled_W runs exactly once (in
    the setup / Cholesky step). To keep the estimand aligned, y must have been
    simulated from this *same* W (see simulate_Cholesky_pooled_withadd).

    Returns (s2a_hat, s2gxg_hat, s2e_hat, AI) -- same ordering as AI_REML.
    """
    n, m = Z.shape

    K = (Z @ Z.T) / m
    Ks = [K, W, np.eye(n)]

    s, AI = reml(y, Ks, iters=iters)

    s2a_hat, s2gxg_hat, s2e_hat = s
    return s2a_hat, s2gxg_hat, s2e_hat, AI
