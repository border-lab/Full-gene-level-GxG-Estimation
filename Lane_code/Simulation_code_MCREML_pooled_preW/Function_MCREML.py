# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.linalg import cholesky
import time

####################################################################
# POOLED within-gene pairwise-epistasis phenotype simulation + MC AI-REML.
#
# Model (no fixed effects; y is mean-centred) -- the "Pooled Model" of
# generative_model.typ:
#     y = g_gxg + e,
#     V = Var(y) = s2gxg * W + s2e * I ,
#     W = (1/P) sum_{g=1}^G H_g H_g' = (1/P) sum_g sum_{a<b in g} h_ab h_ab' ,
#     h_ab = std(Z_a . Z_b),   P = sum_{g=1}^G C(m_g, 2)  (total within-gene pairs).
#
# The kernel takes SEVERAL Z -- a list [Z_1, ..., Z_G], one column-standardized
# genotype block per gene -- and pools ONLY within-gene SNP pairs; cross-gene
# (trans) pairs are excluded because the genes are disjoint SNP sets, so the
# columns of H partition gene-by-gene and H H' = sum_g H_g H_g'.
#
# NORMALIZATION.  The single 1/P division weights every within-gene PAIR
# equally, so a larger gene carries proportionally more heritability (prop m_g^2).
# This is the Pooled Model, distinct from the equal-weight-per-GENE convention
# W = (1/G) sum_g K_g of Simulation_code_MCREML_gxg_1_over_G (= "The Model in
# Notes").  The two coincide EXACTLY when all genes are equal-sized and diverge
# otherwise.  tr(W) = N for both.
#
# SIMULATION builds the pooled W once (build_W_pooled) and forms the Cholesky
# factor Lgxg (Lgxg Lgxg' = s2gxg W) to draw a correctly-correlated epistasis
# effect.  W is DETERMINISTIC per (genotype, gene split), so it is PRE-COMPUTED
# and cached to disk once, then REUSED by every replicate / variance setting.
#
# ESTIMATION (mc_reml / MC_REML) loads that pre-computed dense W and applies it
# as a plain dense mat-vec  W @ B  inside conjugate gradient (O(n^2 c) per CG
# iteration) -- identical to Simulation_code_MCREML_gxg_1_over_G; only the
# kernel W differs.
####################################################################


# ------------------------------------------------------------------ designs
def _standardize_cols(M, stability_std=1e-12):
    """Column-standardize to mean 0, variance 1 (ddof=0).

    A near-constant column (std < stability_std) is left un-scaled to avoid a
    divide-by-zero; such columns contribute ~0 to the GRM anyway.
    """
    M = np.asarray(M, dtype=float)
    mu = M.mean(axis=0)
    sd = M.std(axis=0)
    sd = np.where(sd < stability_std, 1.0, sd)
    return (M - mu) / sd


def additive_design(real_data):
    """Additive design Z_a: column-standardized allele dosages."""
    return _standardize_cols(real_data)


def parse_ratio(ratio, G):
    """Normalize a gene-size ratio spec into a length-G array summing to 1.

    Accepts a comma/space separated string ("0.1,0.2,0.3,0.2,0.2"), a sequence,
    or None / "" / "even" (meaning the equal contiguous split -- returns None,
    which every downstream helper reads as "use np.array_split").  The values
    need not sum to 1; they are rescaled.
    """
    if ratio is None:
        return None
    if isinstance(ratio, str):
        s = ratio.strip()
        if s == "" or s.lower() in ("even", "equal", "none"):
            return None
        parts = [p for p in s.replace(",", " ").split() if p]
        r = np.asarray([float(p) for p in parts], dtype=float)
    else:
        r = np.asarray(ratio, dtype=float)

    if r.size != G:
        raise ValueError(f"--ratio has {r.size} entries but G={G}.")
    if np.any(r <= 0):
        raise ValueError(f"--ratio entries must all be positive; got {r}.")
    return r / r.sum()


def ratio_to_counts(m, ratio):
    """Turn a normalized ratio into G integer SNP counts summing exactly to m.

    Largest-remainder (Hamilton) apportionment: floor every share, then hand the
    leftover SNPs to the genes with the biggest fractional parts.  This keeps
    each gene as close as possible to its requested share of the m SNPs while
    guaranteeing sum(counts) == m and no gene left empty.
    """
    raw = np.asarray(ratio, dtype=float) * m
    counts = np.floor(raw).astype(int)
    rem = m - int(counts.sum())
    if rem > 0:
        order = np.argsort(-(raw - counts))
        counts[order[:rem]] += 1

    if np.any(counts < 1):
        raise ValueError(
            f"ratio {np.asarray(ratio)} leaves an empty gene at m={m} "
            f"(counts={counts}); raise m or the smallest ratio.")
    if np.any(counts < 2):
        print(f"WARNING: gene sizes {counts} include a 1-SNP gene, which "
              f"contributes no within-gene pair and is skipped in W.")
    return counts


def split_into_genes(Z, G, ratio=None):
    """Split the m SNP columns of Z into G contiguous gene blocks.

    Returns a list [Z_1, ..., Z_G] of column sub-matrices -- the SEVERAL Z, one
    per gene, that build_W_pooled consumes.  The blocks are always CONTIGUOUS
    and in SNP order: gene 1 takes the leading columns, gene 2 the next, etc.

    ratio=None gives the equal split (np.array_split: the first (m mod G) genes
    get one SNP extra).  Otherwise ratio is a length-G share vector (see
    parse_ratio) and gene g takes ratio[g] of the m SNPs -- e.g. G=5 with
    ratio=0.1,0.2,0.3,0.2,0.2 and m=1000 cuts the genotype as
    [0:100], [100:300], [300:600], [600:800], [800:1000].  Unequal genes are
    exactly the case the Pooled Model (per-PAIR weight) is built for.
    """
    m = Z.shape[1]
    if ratio is None:
        return [Z[:, cols] for cols in np.array_split(np.arange(m), G)]

    counts = ratio_to_counts(m, ratio)
    bounds = np.concatenate([[0], np.cumsum(counts)])
    return [Z[:, bounds[g]:bounds[g + 1]] for g in range(G)]


def gene_split_tag(G, ratio_str=None):
    """Filename tag identifying a gene split: "G5" or "G5_r0.1-0.2-0.3-0.2-0.2".

    W and the Cholesky factor are deterministic per (genotype, gene split), so
    the split -- not just G -- has to key the caches, or two ratio settings
    would collide on one file.  The tag is built from the RAW --ratio string
    (separators normalized to "-"), never from the rescaled floats, so the
    pipeline shell script can reproduce it with a plain `tr` and the four SLURM
    steps agree on every path.
    """
    if ratio_str is None:
        return f"G{G}"
    s = str(ratio_str).strip()
    if s == "" or s.lower() in ("even", "equal", "none"):
        return f"G{G}"
    return f"G{G}_r" + "-".join(p for p in s.replace(",", " ").split() if p)


# ------------------------------------------------- epistasis GRM (within-gene)
def _within_gene_sum(Zg, pair_batch_size=5000):
    """UN-normalized within-gene epistasis GRM  S = sum_{a<b} h_ab h_ab'.

    h_ab is the column-standardized element-wise product of within-gene SNP
    columns a, b.  Built in pair-batches to bound memory.  Returns (S, p_g):
    the n-by-n sum over the C(m_g, 2) within-gene pairs and the pair count p_g.
    There is NO 1/p_g division here -- the Pooled Model divides ONCE by the
    global pair total P (see build_W_pooled), giving every pair equal weight.
    """
    n, mg = Zg.shape
    pg = mg * (mg - 1) // 2
    idx_i, idx_j = np.triu_indices(mg, k=1)

    S = np.zeros((n, n))
    for start in range(0, pg, pair_batch_size):
        end = min(start + pair_batch_size, pg)
        H = Zg[:, idx_i[start:end]] * Zg[:, idx_j[start:end]]
        mu = H.mean(axis=0)
        sig = H.std(axis=0, ddof=0)
        mask = sig > 1e-10
        H[:, mask] = (H[:, mask] - mu[mask]) / sig[mask]
        H[:, ~mask] = 0.0
        S += H @ H.T
    return S, pg


def build_W_pooled(Z_list, pair_batch_size=5000):
    """Pooled WITHIN-gene pairwise-epistasis GRM (Pooled Model, per-PAIR weight).

    Takes SEVERAL Z -- a list [Z_1, ..., Z_G], one column-standardized genotype
    block per gene -- and pools only within-gene SNP pairs:

        W = (1/P) sum_{g=1}^G H_g H_g' = (1/P) sum_g sum_{a<b in g} h_ab h_ab',
        P = sum_{g=1}^G C(m_g, 2) .

    A 1-SNP gene contributes no pair and is skipped.  The single 1/P division
    weights every within-gene pair EQUALLY (so a larger gene carries prop m_g^2
    more heritability) -- the Pooled Model of generative_model.typ, distinct
    from the equal-weight-per-gene W = (1/G) sum_g K_g of the _1_over_G pipeline;
    the two coincide only for equal-sized genes.  Deterministic per (genotype,
    gene split); built once, cached, applied densely as W @ b at estimation.
    tr(W) = N.
    """
    n = Z_list[0].shape[0]
    W = np.zeros((n, n))
    P = 0
    for Zg in Z_list:
        if Zg.shape[1] < 2:                  # a 1-SNP gene has no within-gene pair
            continue
        S, pg = _within_gene_sum(Zg, pair_batch_size=pair_batch_size)
        W += S
        P += pg
    if P == 0:
        raise ValueError("No gene block has >= 2 SNPs; increase m/G.")
    W /= P
    return W


# ------------------------------------------------------------- simulation
def simulate_Cholesky_gxg(real_data, G, s2gxg=0.5, s2e=0.5, stability=1e-10,
                          ratio=None):
    """Cholesky factor of the pooled epistasis covariance, plus the GRM W itself.

    The m SNPs of real_data are column-standardized and split into G contiguous
    genes (split_into_genes -> several Z), and the pooled within-gene kernel is
    built by build_W_pooled.  ratio (see parse_ratio / split_into_genes) sets
    the per-gene share of SNPs; None keeps the equal split.
    Returns (Lgxg, W, w_build_time) with

        Lgxg Lgxg' = s2gxg W .

    W (deterministic per (genotype, G)) is returned too so the caller can cache
    it once for the estimation step's dense mat-vec.  w_build_time is the
    wall-clock seconds spent building W only (tracked separately from the
    downstream estimation time).
    """
    Za = additive_design(real_data)
    n, m = Za.shape

    genes = split_into_genes(Za, G, ratio=ratio)   # several Z, one per gene
    print(f"gene sizes: {[g.shape[1] for g in genes]}")

    t_start = time.perf_counter()
    W = build_W_pooled(genes)
    w_build_time = time.perf_counter() - t_start

    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)
    return Lgxg, W, w_build_time


def simulate_remove_sampling_err(Lgxg, n, s2gxg=0.5, s2e=0.5):
    """Phenotype y = g_gxg + e with per-component sampling-error removal.

    g_gxg = Lgxg u1  (~ N(0, s2gxg W)),  e white noise.  Each component is
    rescaled to its exact target variance and y is mean-centred, so the fitted
    model carries no fixed effect.  Returns y (n,).
    """
    u1 = np.random.randn(n)
    u2 = np.random.randn(n)

    gxg = Lgxg @ u1                       # epistasis effect ~ N(0, s2gxg W)
    e = np.sqrt(s2e) * u2                 # residual noise

    # Eliminate sampling variances: rescale each component to its exact target.
    gxg *= np.sqrt(s2gxg / np.var(gxg, ddof=0))
    e *= np.sqrt(s2e / np.var(e, ddof=0))

    y = gxg + e
    y -= y.mean()
    return y


# ------------------------------------------------ matrix-free linear algebra
def _v_matvec(W, s2gxg, s2e, B):
    """Apply V = s2gxg W + s2e I to B using the pre-computed dense W.

        V B = s2gxg (W B) + s2e B .

    B may be (n,) or (n, c); the result matches its shape.  The dense product
    W @ B is O(n^2 c) and dominates; nothing is rebuilt from the genotype.
    """
    return s2gxg * (W @ B) + s2e * B


def _cg_batched(matvec, Bmat, x0=None, tol=1e-6, maxiter=1000):
    """Conjugate gradient for the SPD system V X = Bmat.

    matvec : callable  X -> V X   (accepts / returns (n, c) arrays).
    Bmat   : (n, c) right-hand sides -- all c columns solved together with
             per-column CG scalars, so one V-pass advances every column.
    x0     : (n, c) warm start (e.g. the previous REML iteration's solution).
    """
    n, c = Bmat.shape
    X = np.zeros((n, c)) if x0 is None else x0.copy()
    R = Bmat - matvec(X)
    P = R.copy()
    rs_old = np.sum(R * R, axis=0)
    b_norm = np.sqrt(np.sum(Bmat * Bmat, axis=0))
    b_norm[b_norm == 0.0] = 1.0

    for _ in range(maxiter):
        VP = matvec(P)
        alpha = rs_old / np.sum(P * VP, axis=0)
        X += alpha * P
        R -= alpha * VP
        rs_new = np.sum(R * R, axis=0)
        if np.max(np.sqrt(rs_new) / b_norm) < tol:
            break
        beta = rs_new / rs_old
        P = R + beta * P
        rs_old = rs_new
    return X


# ----------------------------------------------------------- MC AI-REML (k=2)
def mc_reml(W, y, iters=30, Nmc=50, cg_tol=1e-6, cg_maxiter=1000,
            jitter=1e-8, tol=1e-8, lm=1e-3, step_frac=0.5, upper_mult=5.0,
            seed=None, verbose=False):
    """Monte-Carlo AI-REML for V = s2gxg W + s2e I (dense pre-computed W).

    Per iteration: 1 CG solve for x = V^{-1} y, Nmc solves for the Hutchinson
    probes, and 2 solves for V^{-1}(K_j x).  W enters only as the dense product
    W @ B; W U is fixed and formed once.

    Returns
    -------
    s  : (2,) estimated (s2gxg, s2e).
    AI : (2, 2) final average-information matrix.
    """
    y = np.asarray(y, dtype=float).flatten()
    W = np.asarray(W, dtype=float)
    n = y.shape[0]
    k = 2

    vary = y.var()
    s_upper = upper_mult * vary          # no component can exceed ~total var

    rng = np.random.default_rng(seed)
    U = rng.choice([-1.0, 1.0], size=(n, Nmc))      # Rademacher probes in {+-1}

    s = np.full(k, vary / k)
    AI = np.eye(k)

    yc = y.reshape(n, 1)
    xbuf = None       # warm-start buffers for the three CG solve groups
    P = None
    G = None

    WU = W @ U                            # W U for the fixed probes: once

    for it in range(iters):
        s2gxg, s2e = s
        matvec = lambda B: _v_matvec(W, s2gxg, s2e, B)

        # --- x = V^{-1} y ---
        xbuf = _cg_batched(matvec, yc, x0=xbuf, tol=cg_tol, maxiter=cg_maxiter)
        x = xbuf[:, 0]

        # data quadratics x'K_i x
        Wx = W @ x
        xWx = x @ Wx
        xIx = x @ x

        # --- Hutchinson trace: P = V^{-1} U ---
        P = _cg_batched(matvec, U, x0=P, tol=cg_tol, maxiter=cg_maxiter)
        trV1W = np.mean(np.sum(P * WU, axis=0))
        trV1I = np.mean(np.sum(P * U, axis=0))

        score = np.array([0.5 * (xWx - trV1W),
                          0.5 * (xIx - trV1I)])

        # --- average information: A_ij = 0.5 (K_i x)' V^{-1}(K_j x) ---
        KX = np.column_stack([Wx, x])
        G = _cg_batched(matvec, KX, x0=G, tol=cg_tol, maxiter=cg_maxiter)
        AI = 0.5 * (KX.T @ G)
        AI = 0.5 * (AI + AI.T)

        # --- damped, bounded AI-Newton step ---------------------------------
        # W is often nearly collinear with I (interactions of standardized SNPs
        # are ~independent), so the AI matrix can be near-singular and an
        # undamped step explodes.  Levenberg-Marquardt ridge (scaled to AI) +
        # trust region on the step + a box clamp keep the path stable without
        # perturbing well-identified cases.
        dA = np.abs(np.diag(AI))
        ridge = lm * (dA.mean() + 1e-12)
        step = np.linalg.solve(AI + (ridge + jitter) * np.eye(k), score)

        mx = np.abs(step).max()
        max_step = step_frac * vary                  # trust region
        if mx > max_step:
            step *= max_step / mx
        s = np.clip(s + step, 1e-9, s_upper)         # box clamp

        if verbose:
            print(f"iter {it:2d}  s={s}  max|step|={np.abs(step).max():.3e}")
        if np.abs(step).max() < tol:
            break

    return s, AI


def MC_REML(W, y, iters=30, Nmc=50, cg_tol=1e-6, cg_maxiter=1000, seed=None):
    """Wrapper: returns (s2gxg_hat, s2e_hat, AI)."""
    s, AI = mc_reml(W, y, iters=iters, Nmc=Nmc, cg_tol=cg_tol,
                    cg_maxiter=cg_maxiter, seed=seed)
    s2gxg_hat, s2e_hat = s
    return s2gxg_hat, s2e_hat, AI
