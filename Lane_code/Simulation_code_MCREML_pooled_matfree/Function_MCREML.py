# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.linalg import cholesky
import time


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


def split_into_genes(Z, G):
    """Split the m SNP columns of Z into G contiguous gene blocks.

    Returns a list [Z_1, ..., Z_G] of column sub-matrices -- the SEVERAL Z, one
    per gene, that the pooled kernel consumes.  np.array_split handles m not
    divisible by G by making the first (m mod G) genes one SNP larger, so the
    genes may be unequal -- exactly the case the Pooled Model is built for.
    """
    m = Z.shape[1]
    return [Z[:, cols] for cols in np.array_split(np.arange(m), G)]


# ------------------------------------- epistasis GRM (dense, SIMULATION only)
def _within_gene_sum(Zg, pair_batch_size=5000):
    """UN-normalized within-gene epistasis GRM  S = sum_{a<b} h_ab h_ab'.

    h_ab is the RAW element-wise product of within-gene SNP columns a, b:

        h_ab = Z_a .* Z_b        (NO mean-centering, NO 1/sigma_ab scaling).

    Built in pair-batches to bound memory.  Returns (S, p_g).  There is NO
    1/p_g division here -- the Pooled Model divides ONCE by the global pair
    total P (see build_W_pooled), giving every pair equal weight.

    The degenerate-pair mask that the standardized version needed is GONE, and
    must be: it existed only to avoid dividing by a near-zero sigma_ab.  With no
    division a near-constant product column is a perfectly good rank-one
    contribution, so every pair now enters unconditionally.  The matrix-free
    apply has no mask either -- it sums over all a != b (see compute_WU_pooled)
    -- which is what keeps the two sides equal.
    """
    n, mg = Zg.shape
    pg = mg * (mg - 1) // 2
    idx_i, idx_j = np.triu_indices(mg, k=1)

    S = np.zeros((n, n))
    for start in range(0, pg, pair_batch_size):
        end = min(start + pair_batch_size, pg)
        H = Zg[:, idx_i[start:end]] * Zg[:, idx_j[start:end]]
        S += H @ H.T
    return S, pg


def build_W_pooled(Z_list, pair_batch_size=5000):
    """Pooled WITHIN-gene pairwise-epistasis GRM (Pooled Model, per-PAIR weight).

        W = (1/P) sum_{g=1}^G H_g H_g' ,   P = sum_{g=1}^G C(m_g, 2) ,

    with H_g the RAW within-gene interaction design, columns h_ab = Z_a .* Z_b.

    A 1-SNP gene contributes no pair and is skipped.  O(n^2 P) time and n-by-n
    storage.  USED ONLY BY THE SIMULATION (the Cholesky factor needs an explicit
    matrix); estimation applies W matrix-free via compute_WU_pooled.

    tr(W) = n and W 1 = 0 were both consequences of standardizing h_ab and
    NEITHER survives: tr(W)/n is now 1 only in expectation under independent
    SNPs, and 1 is no longer a null vector, so lambda_min(W) > 0 in general
    (W is still PSD, being a Gram matrix).  Nothing downstream in this pipeline
    relied on either identity.
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


# ------------------------------------------- matrix-free pooled W apply
def setup_pooled(Z_list):
    """SETUP phase: everything that depends on the genotype alone.

    Returns (D_list, P) with

        D_list : [D_1, ..., D_G]  (None for a 1-SNP gene), D_g = Z_g .* Z_g,
                 n-by-m_g each -- the u-INDEPENDENT half of the operator
        P      : sum_g C(m_g, 2)                          -- global pair total

    The RAW interaction column removes the pair weights entirely: there is no
    1/sigma_ab^2 scaling to invert and no mu_ab centering to subtract, so the
    m_g-by-m_g weight matrices (V_g, R_g, T_g) and the pooled centering pair
    (v_R, s_T) that the standardized kernel needed are GONE -- not zeroed and
    carried, gone.  What is left to precompute is only D_g, which is what the
    D D' term of the operator needs on every apply.

    O(n m) time and storage -- no m-by-m object is built at any point, in setup
    or in the apply.
    """
    D_list = []
    P = 0

    for Zg in Z_list:
        mg = Zg.shape[1]
        if mg < 2:                               # no within-gene pair
            D_list.append(None)
            continue
        P += mg * (mg - 1) // 2
        D_list.append(Zg * Zg)                   # D_g = Z_g .* Z_g, u-independent

    if P == 0:
        raise ValueError("No gene block has >= 2 SNPs; increase m/G.")
    return D_list, P


def compute_WU_pooled(Z_list, D_list, P, U):
    """Matrix-free product  W @ U  for the pooled within-gene epistasis GRM.

    U may be (n,) or (n, c); the result matches its shape.  W (n-by-n) is never
    formed, nor any interaction block H_g (n-by-p_g), nor -- and this is the
    change from the standardized version -- anything m-by-m.  Per gene, with
    K_w = Z_g Z_g' the additive GRM (also never formed) and D_g = Z_g .* Z_g:

        D_g D_g' u    = (Z_g .* Z_g)((Z_g .* Z_g)' u)              O(n m_g)
        (K_w .* K_w)u = (Z_g (Z_g'(u .* Z_g)) .* Z_g) 1            O(n m_g^2)

        W u = 1/(2P) sum_g [ (K_w .* K_w) u  -  D_g D_g' u ]       O(n m_g^2)

    WHY THIS IS THE WHOLE OPERATOR.  Summing the raw pair columns over UNORDERED
    pairs and doubling,

        2 sum_{a<b} (Z_a .* Z_b)(Z_a .* Z_b)' = sum_{a,b} (...) - sum_{a=b} (...)
                                              = (K_w .* K_w) - D_g D_g' ,

    since sum_{a,b} Z_ta Z_tb Z_sa Z_sb = (sum_a Z_ta Z_sa)^2 = (K_w)_ts^2.  The
    D_g D_g' term is exactly the a = b diagonal the pair sum excludes.  There is
    no weight matrix and no centering term left: dropping the standardization
    removed both, so the two lines above ARE W u.

    IMPLEMENTATION.  M_g = Z_g'(u .* Z_g) is the one u-dependent contraction and
    is column-specific, so the columns are looped; each iteration is two real
    gemms of n m_g^2 flops, and only the m_g-by-m_g M_g is allocated.  The
    standardized version instead built the n-by-m_g^2 outer-product design
    D_g[i, a*m_g+b] = Z_g[i,a] Z_g[i,b] TWICE per apply, to batch the same
    contraction over columns in one gemm.  That was worth it there because the
    elementwise V_g .* M_g weight had to be applied in the m_g^2 basis anyway.
    Here there is no weight, so the m_g^2 basis buys nothing while costing a
    buf_elems-sized buffer and two passes over it -- and, crucially, that cost
    is FIXED in the number of columns, so it is ruinous at the narrow widths the
    REML loop actually spends its time on.  Measured against the old route,
    identical output to 2e-16:

        n=1000, m=1000, G=10        n=2000, m=1000, G=10
        width  1   0.004 vs 0.370s  (93x)     0.017 vs 0.747s  (45x)
        width  2   0.009 vs 0.379s  (42x)     0.041 vs 0.780s  (19x)
        width 50   0.198 vs 0.485s  (2.5x)    0.842 vs 0.967s  (1.1x)

    The REML iteration issues width-1 applies inside every CG step for x, plus
    width-2 for the AI matrix, so the left column is what dominates a run.

    Equals build_W_pooled(Z_list) @ U to float precision.
    """
    n = Z_list[0].shape[0]
    U = np.asarray(U, dtype=float)
    single = (U.ndim == 1)
    if single:
        U = U.reshape(n, 1)
    c = U.shape[1]

    T1 = np.zeros((n, c))
    for Zg, Dg in zip(Z_list, D_list):
        if Dg is None:                           # 1-SNP gene: no pair, no work
            continue
        for k in range(c):
            M = Zg.T @ (U[:, k:k + 1] * Zg)      # m_g-by-m_g,  Z_g'(u .* Z_g)
            T1[:, k] += np.einsum('na,na->n', Zg @ M, Zg)     # (Z_g M .* Z_g) 1
        T1 -= Dg @ (Dg.T @ U)                    # the a = b terms, O(n m_g c)

    out = T1 / (2.0 * P)
    return out[:, 0] if single else out


# ------------------------------------------------------------- simulation
def simulate_Cholesky_gxg(real_data, G, s2gxg=0.5, s2e=0.5, stability=1e-10):
    """Cholesky factor of the pooled epistasis covariance.

    The m SNPs of real_data are column-standardized and split into G contiguous
    genes (split_into_genes -> several Z), and the pooled within-gene kernel is
    built densely by build_W_pooled.  Returns (Lgxg, w_build_time) with

        Lgxg Lgxg' = s2gxg W .

    The dense W is used here and then discarded: unlike the _preW pipeline it is
    NOT cached, because estimation rebuilds its action matrix-free from the
    genotype.  w_build_time is the wall-clock seconds spent building W only.
    """
    Za = additive_design(real_data)
    n, m = Za.shape

    genes = split_into_genes(Za, G)          # several Z, one per gene

    t_start = time.perf_counter()
    W = build_W_pooled(genes)
    w_build_time = time.perf_counter() - t_start

    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)
    return Lgxg, w_build_time


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
def _v_matvec(Z_list, D_list, P, s2gxg, s2e, B):
    """Apply V = s2gxg W + s2e I to B, with W applied matrix-free.

        V B = s2gxg (W B) + s2e B .

    B may be (n,) or (n, c); the result matches its shape.  The epistasis term
    is rebuilt from the genotype on every call -- no n-by-n W exists.
    """
    return s2gxg * compute_WU_pooled(Z_list, D_list, P, B) + s2e * B


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
def mc_reml(Z, y, G, iters=30, Nmc=50, cg_tol=1e-6, cg_maxiter=1000,
            jitter=1e-8, tol=1e-8, lm=1e-3, step_frac=0.5, upper_mult=5.0,
            seed=None, verbose=False):
    """Monte-Carlo AI-REML for V = s2gxg W + s2e I (pooled W, matrix-free).

    Z is the column-standardized genotype; it is split into G contiguous genes
    and the pooled setup (D_g = Z_g .* Z_g per gene, and the pair total P) is
    built ONCE here, outside the iteration -- that is the u-independent
    precomputation of Wu_complexity.typ.  With the RAW interaction column there
    are no pair-weight matrices left to precompute, so setup is O(n m).

    Per iteration: 1 CG solve for x = V^{-1} y, Nmc solves for the Hutchinson
    probes, and 2 solves for V^{-1}(K_j x).  W U for the FIXED probes is formed
    once up front, exactly as in the dense pipeline.

    Returns
    -------
    s  : (2,) estimated (s2gxg, s2e).
    AI : (2, 2) final average-information matrix.
    """
    y = np.asarray(y, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float)
    n = y.shape[0]
    k = 2

    # --- setup: genotype-only, done once (D_g = Z_g .* Z_g per gene) ---
    genes = split_into_genes(Z, G)
    D_list, P = setup_pooled(genes)
    Wapply = lambda B: compute_WU_pooled(genes, D_list, P, B)

    vary = y.var()
    s_upper = upper_mult * vary          # no component can exceed ~total var

    rng = np.random.default_rng(seed)
    U = rng.choice([-1.0, 1.0], size=(n, Nmc))      # Rademacher probes in {+-1}

    s = np.full(k, vary / k)
    AI = np.eye(k)

    yc = y.reshape(n, 1)
    xbuf = None       # warm-start buffers for the three CG solve groups
    Pbuf = None
    Gbuf = None

    WU = Wapply(U)                        # W U for the fixed probes: once

    for it in range(iters):
        s2gxg, s2e = s
        matvec = lambda B: _v_matvec(genes, D_list, P, s2gxg, s2e, B)

        # --- x = V^{-1} y ---
        xbuf = _cg_batched(matvec, yc, x0=xbuf, tol=cg_tol, maxiter=cg_maxiter)
        x = xbuf[:, 0]

        # data quadratics x'K_i x
        Wx = Wapply(x)
        xWx = x @ Wx
        xIx = x @ x

        # --- Hutchinson trace: Pmat = V^{-1} U ---
        Pbuf = _cg_batched(matvec, U, x0=Pbuf, tol=cg_tol, maxiter=cg_maxiter)
        trV1W = np.mean(np.sum(Pbuf * WU, axis=0))
        trV1I = np.mean(np.sum(Pbuf * U, axis=0))

        score = np.array([0.5 * (xWx - trV1W),
                          0.5 * (xIx - trV1I)])

        # --- average information: A_ij = 0.5 (K_i x)' V^{-1}(K_j x) ---
        KX = np.column_stack([Wx, x])
        Gbuf = _cg_batched(matvec, KX, x0=Gbuf, tol=cg_tol, maxiter=cg_maxiter)
        AI = 0.5 * (KX.T @ Gbuf)
        AI = 0.5 * (AI + AI.T)

        # --- damped, bounded AI-Newton step ---------------------------------
        # W can be nearly collinear with I, making the AI matrix near-singular
        # so that an undamped step explodes.  Levenberg-Marquardt ridge (scaled
        # to AI) + trust region on the step + a box clamp keep the path stable
        # without perturbing well-identified cases.  The RAW interaction column
        # weights each pair by its own sigma_ab instead of normalising it away,
        # which spreads the spectrum of W (measured [0.12, 3.97] at n=200, m=60,
        # G=3, against tr(W)/n = 1.0) and moves it away from I -- so this
        # damping should bind less often than it did with the standardized
        # kernel.  It is unchanged regardless.
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


def MC_REML(Z, y, G, iters=30, Nmc=50, cg_tol=1e-6, cg_maxiter=1000, seed=None):
    """Wrapper: returns (s2gxg_hat, s2e_hat, AI)."""
    s, AI = mc_reml(Z, y, G, iters=iters, Nmc=Nmc, cg_tol=cg_tol,
                    cg_maxiter=cg_maxiter, seed=seed)
    s2gxg_hat, s2e_hat = s
    return s2gxg_hat, s2e_hat, AI
