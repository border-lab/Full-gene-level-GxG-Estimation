# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.linalg import cholesky
import time

####################################################################
# POOLED within-gene pairwise-epistasis phenotype simulation + MC AI-REML
# (pre-computed dense W variant).
#
# Model (no fixed effects; y is mean-centred):
#     y = g_gxg + e,
#     V = Var(y) = s2gxg * W + s2e * I ,
#     W = (1/G) sum_{g=1}^G K_g   (POOLED within-gene GRM),
#     K_g = (1/p_g) sum_{a<b in gene g} h_ab h_ab' ,  h_ab = std(Z_a . Z_b),
#     p_g = C(m/G, 2)  within-gene pairs per gene.
#
# The m SNPs are split into G contiguous gene blocks of m/G SNPs each (e.g.
# G = 10, m = 1000 -> ten 100-SNP genes).  Only WITHIN-gene SNP pairs enter W;
# cross-gene pairs are excluded.  For G = 1 this reduces to the genome-wide
# all-pairs kernel of Simulation_code_MCREML_gxg.
#
# Z_a : column-standardized allele dosages.
#
# SIMULATION builds the pooled W once (build_W_pooled_within_gene) and forms the
# Cholesky factor Lgxg (Lgxg Lgxg' = s2gxg W) to draw a correctly-correlated
# epistasis effect.  W is DETERMINISTIC per (genotype, G), so it is PRE-COMPUTED
# and cached to disk once, then REUSED by every replicate / variance setting.
#
# ESTIMATION (mc_reml / MC_REML) loads that pre-computed dense W and applies it
# as a plain dense mat-vec  W @ B  inside conjugate gradient (O(n^2 c) per CG
# iteration) -- identical to Simulation_code_MCREML_gxg; only the kernel W
# differs.  See MCREML.typ for the derivation.
#
# SIMULATE WHOLE / ESTIMATE PART.  The phenotype is ALWAYS simulated from the
# full G-gene kernel W_full = (1/G) sum_g K_g, but the REML fit can use a kernel
# pooled over a SUBSET of genes, W_est = (1/G_est) sum_{g in subset} K_g (see
# build_W_pooled_within_gene_with_subset / parse_gene_subset).  That is a
# deliberately MISSPECIFIED fit; because each gene enters W with the SAME 1/G
# weight, the recovered s2gxg is attenuated towards s2gxg * (G_est / G_full) --
# the subset's SHARE OF GENES (count-based), independent of gene sizes.  This is
# the equal-weight-per-gene analogue of the per-PAIR share P_est/P_full of
# Simulation_code_MCREML_pooled_preW.
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


# ----------------------------------------------------- gene split / cache tags
def parse_ratio(ratio, G):
    """Normalize a gene-size ratio spec into a length-G array summing to 1.

    Accepts a comma/space separated string ("0.1,0.2,0.3,0.2,0.2"), a sequence,
    or None / "" / "even" (meaning the equal contiguous split -- returns None,
    which every downstream helper reads as "use np.array_split").  The values
    need not sum to 1; they are rescaled.  The ratio only sets which SNPs land in
    which gene; under the equal-weight-per-gene (1/G) kernel each gene still
    contributes the SAME weight regardless of its size.
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
    per gene.  The blocks are always CONTIGUOUS and in SNP order: gene 1 takes
    the leading columns, gene 2 the next, etc.

    ratio=None gives the equal split (np.array_split: the first (m mod G) genes
    get one SNP extra) -- exactly the block layout of the original
    build_W_pooled_within_gene.  Otherwise ratio is a length-G share vector (see
    parse_ratio) and gene g takes ratio[g] of the m SNPs -- e.g. G=5 with
    ratio=0.1,0.2,0.3,0.2,0.2 and m=1000 cuts the genotype as
    [0:100], [100:300], [300:600], [600:800], [800:1000].
    """
    m = Z.shape[1]
    if ratio is None:
        return [Z[:, cols] for cols in np.array_split(np.arange(m), G)]

    counts = ratio_to_counts(m, ratio)
    bounds = np.concatenate([[0], np.cumsum(counts)])
    return [Z[:, bounds[g]:bounds[g + 1]] for g in range(G)]


def ratio_tag(ratio_str):
    """Filename suffix identifying the gene-SIZE split: "" or "_r0.1-0.2-...".

    W and the Cholesky factor are deterministic per (genotype, gene split), so an
    unequal split has to key the caches or two ratio settings would collide.  The
    tag is built from the RAW --ratio string (separators normalized to "-"),
    never from the rescaled floats, so the pipeline shell script can reproduce it
    with a plain `tr` and the four SLURM steps agree on every path.  Empty for
    the even split -- keeping every even-split path byte-identical to the
    pipeline's previous behaviour.  (The gene COUNT G already sits in the base
    tag, so only the ratio part is added here.)
    """
    if ratio_str is None:
        return ""
    s = str(ratio_str).strip()
    if s == "" or s.lower() in ("even", "equal", "none"):
        return ""
    return "_r" + "-".join(p for p in s.replace(",", " ").split() if p)


def parse_gene_subset(spec, G):
    """1-based gene indices to USE AT ESTIMATION -> sorted 0-based list.

    "1,2" means the estimation kernel is built from the front two genes only,
    while the phenotype is still simulated from ALL G genes -- a deliberately
    MISSPECIFIED fit that asks how much of s2gxg a partial gene panel recovers.
    None / "" / "all" returns None, meaning "use every gene" (the correctly
    specified fit, i.e. the pipeline's previous behaviour).

    An EXPLICIT full list ("1,2,3,4,5" at G=5) is NOT collapsed to None: it
    still returns every index, so W_est is built and cached under the matching
    _est1-2-3-4-5 tag.  Collapsing it would make this function disagree with the
    pipeline's shell-side tag (which cannot tell "all genes" from "some genes"),
    and the estimation step would then write results to the UN-suffixed
    directory while the combine step read the suffixed one -- yielding an empty
    combined .txt next to a full result dir.  Only an EMPTY spec means "all".
    """
    if spec is None:
        return None
    if isinstance(spec, str):
        s = spec.strip()
        if s == "" or s.lower() in ("all", "none"):
            return None
        parts = [p for p in s.replace(",", " ").split() if p]
        idx = [int(p) for p in parts]
    else:
        idx = [int(i) for i in spec]

    if len(idx) == 0:
        return None
    bad = [i for i in idx if i < 1 or i > G]
    if bad:
        raise ValueError(f"--estimate indices {bad} out of range 1..{G}.")
    return [i - 1 for i in sorted(set(idx))]     # to 0-based


def gene_subset_tag(spec, G):
    """Filename suffix for the estimation subset: "" (all genes) or "_est1-2".

    Built from the RAW --estimate string so the pipeline shell script can
    reproduce it with `tr`, exactly like ratio_tag.  Empty ONLY when the spec
    itself is empty -- an explicit "1,2,3,4,5" still tags as _est1-2-3-4-5,
    because the shell side cannot detect "this list happens to be every gene"
    and the two must agree on every path.  G is accepted for signature symmetry
    with parse_gene_subset but deliberately unused.
    """
    if spec is None:
        return ""
    s = str(spec).strip()
    if s == "" or s.lower() in ("all", "none"):
        return ""
    parts = [p for p in s.replace(",", " ").split() if p]
    return "_est" + "-".join(parts)


# ------------------------------------------------- epistasis GRM (from onlyW)
def build_W_batched(Z, pair_batch_size=5000):
    """Explicit pairwise-epistasis GRM  W = (1/p) sum_{a<b} h_ab h_ab'.

    h_ab is the column-standardized element-wise product of SNP columns a, b.
    Built in pair-batches to bound memory.  O(n^2 p) time, n-by-n storage.
    Computed ONCE per genotype and cached to disk; reused everywhere else.
    """
    n, m = Z.shape
    p = m * (m - 1) // 2
    idx_i, idx_j = np.triu_indices(m, k=1)

    W = np.zeros((n, n))
    for start in range(0, p, pair_batch_size):
        end = min(start + pair_batch_size, p)
        H = Z[:, idx_i[start:end]] * Z[:, idx_j[start:end]]
        mu = H.mean(axis=0)
        sig = H.std(axis=0, ddof=0)
        mask = sig > 1e-10
        H[:, mask] = (H[:, mask] - mu[mask]) / sig[mask]
        H[:, ~mask] = 0.0
        W += H @ H.T
    W /= p
    return W


def build_W_pooled_within_gene(Z, G, pair_batch_size=5000):
    """Pooled WITHIN-GENE pairwise-epistasis GRM.

    The m SNP columns are split into G contiguous gene blocks (e.g. G = 10 and
    m = 1000 -> ten 100-SNP genes).  Each block g contributes only its OWN
    within-gene epistasis GRM K_g = build_W_batched(Z_g) (standardized products
    of within-gene SNP pairs), and the pooled kernel is the uniform average

        W = (1/G) sum_{g=1}^G K_g .

    Cross-gene SNP pairs are excluded, so W uses G * C(m/G, 2) within-gene pairs
    instead of the genome-wide C(m, 2).  Deterministic per (genotype, G); built
    once and cached, then applied densely as W @ b at estimation time.  For
    G = 1 this is exactly build_W_batched (the genome-wide all-pairs kernel).
    """
    n, m = Z.shape
    blocks = np.array_split(np.arange(m), G)     # G contiguous gene blocks
    W = np.zeros((n, n))
    n_genes = 0
    for cols in blocks:
        if cols.size < 2:                        # a 1-SNP gene has no pair
            continue
        W += build_W_batched(Z[:, cols], pair_batch_size=pair_batch_size)
        n_genes += 1
    if n_genes == 0:
        raise ValueError("No gene block has >= 2 SNPs; increase m/G.")
    W /= n_genes
    return W


def build_W_pooled_within_gene_with_subset(Z_list, subset, pair_batch_size=5000):
    """Full and SUBSET pooled within-gene GRMs, equal-weight-per-gene (1/G).

    W_full = (1/G_full) sum_{g: m_g>=2} K_g            -- simulate the phenotype from this
    W_est  = (1/G_est ) sum_{g in subset, m_g>=2} K_g  -- fit the REML model with this

    K_g = build_W_batched(Z_g) is each gene's OWN pair-normalized epistasis GRM
    (tr(K_g) = N), and both pooled kernels AVERAGE their K_g -- every gene enters
    with the same 1/G weight regardless of size (tr = N for each; this is "The
    Model in Notes", distinct from the per-PAIR Pooled Model of
    Simulation_code_MCREML_pooled_preW).

    Because the truth is s2gxg * W_full and the fit spans only the subset's
    genes, the fitted s2gxg is attenuated: if the omitted genes' K_g are
    near-orthogonal to the retained ones, E[s2gxg_hat] ~ s2gxg * (G_est / G_full)
    -- the subset's SHARE OF GENES (count-based, NOT the pair share of the preW
    model, because here every gene weighs the same).  Multiply the estimate by
    G_full / G_est to put it back on the full-panel scale.

    subset is a 0-based index list (parse_gene_subset); None gives W_est = None.
    Returns (W_full, W_est, G_full, G_est).  Each K_g is formed once and
    accumulated into both averages, so the subset costs one extra n-by-n array
    and no extra genotype work.
    """
    n = Z_list[0].shape[0]
    sel = None if subset is None else set(subset)

    W_full = np.zeros((n, n))
    W_est = None if sel is None else np.zeros((n, n))
    G_full = 0
    G_est = 0

    for g, Zg in enumerate(Z_list):
        if Zg.shape[1] < 2:                  # a 1-SNP gene has no within-gene pair
            continue
        Kg = build_W_batched(Zg, pair_batch_size=pair_batch_size)
        W_full += Kg
        G_full += 1
        if sel is not None and g in sel:
            W_est += Kg
            G_est += 1

    if G_full == 0:
        raise ValueError("No gene block has >= 2 SNPs; increase m/G.")
    W_full /= G_full
    if sel is not None:
        if G_est == 0:
            raise ValueError(
                "No gene in --estimate has >= 2 SNPs, so the estimation kernel "
                "would be empty; pick genes with more SNPs.")
        W_est /= G_est
    return W_full, W_est, G_full, G_est


# ------------------------------------------------------------- simulation
def simulate_Cholesky_gxg(real_data, s2gxg=0.5, s2e=0.5, stability=1e-10):
    """Cholesky factor of the epistasis covariance, plus the GRM W itself.

    Returns (Lgxg, W) with
        Lgxg Lgxg' = s2gxg W .
    W (the epistasis GRM, deterministic per genotype) is returned too so the
    caller can cache it once for the estimation step's dense mat-vec.  Built
    once per (genotype, s2gxg) and reused across reps.
    """
    Za = additive_design(real_data)
    n, m = Za.shape

    W = build_W_batched(Za)
    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)
    return Lgxg, W


def simulate_Cholesky_pooled(real_data, G, s2gxg=0.5, s2e=0.5, stability=1e-10,
                             ratio=None, subset=None):
    """Cholesky factor of the POOLED within-gene epistasis covariance, plus GRMs.

    The m SNPs of real_data are column-standardized and split into G contiguous
    genes (split_into_genes -> several Z), and the equal-weight-per-gene kernels
    are built by build_W_pooled_within_gene_with_subset.  ratio (see parse_ratio)
    sets the per-gene share of SNPs; None keeps the equal split -- identical to
    the previous build_W_pooled_within_gene layout.

    The phenotype is ALWAYS simulated from the full-panel kernel:

        Lgxg Lgxg' = s2gxg W_full ,   W_full = (1/G_full) sum_g K_g .

    subset (0-based gene indices, see parse_gene_subset) additionally builds the
    ESTIMATION kernel W_est = (1/G_est) sum_{g in subset} K_g from those genes
    only -- the misspecified fit.  It never touches Lgxg, so simulation is
    unchanged by it.

    Returns (Lgxg, W_full, W_est, info) where W_est is None when subset is None
    and info is a dict with the gene sizes, contributing-gene counts G_full /
    G_est, the expected attenuation gene_share = G_est/G_full, and the W build
    time in seconds (tracked separately from the downstream estimation time).
    """
    Za = additive_design(real_data)
    n, m = Za.shape

    genes = split_into_genes(Za, G, ratio=ratio)   # several Z, one per gene
    sizes = [g.shape[1] for g in genes]

    t_start = time.perf_counter()
    W_full, W_est, G_full, G_est = build_W_pooled_within_gene_with_subset(
        genes, subset)
    w_build_time = time.perf_counter() - t_start

    info = {
        "gene_sizes": sizes,
        "subset": None if subset is None else [g + 1 for g in subset],
        "G_full": G_full,
        "G_est": G_est if subset is not None else G_full,
        "gene_share": 1.0 if subset is None else G_est / G_full,
        "w_build_time": w_build_time,
    }

    Lgxg = cholesky(s2gxg * W_full + stability * np.eye(n), lower=True)
    return Lgxg, W_full, W_est, info


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
