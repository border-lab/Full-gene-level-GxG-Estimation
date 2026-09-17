# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.linalg import cholesky
import time

####################################################################
# Pairwise-epistasis-ONLY phenotype simulation + Monte-Carlo AI-REML.
# MEAN-CENTRED-KERNEL variant of Simulation_code_MCREML_gxg.
#
# Model (no fixed effects; y is mean-centred):
#     y = g_gxg + e,
#     V = Var(y) = s2gxg * W + s2e * I ,
#     W = (1/p) sum_{a<b} h_ab h_ab' ,  h_ab = Z_a.Z_b - mean(Z_a.Z_b) 1 ,
#     p = m(m-1)/2               (pairwise-epistasis GRM).
#
# Z_a : column-standardized allele dosages.
# W   : GRM of the MEAN-CENTRED element-wise products of every SNP pair.
#
# THE ONLY DIFFERENCE FROM MCREML_gxg is the kernel: the interaction column
# h_ab is mean-centred but NOT variance-scaled -- only the 1/sigma_ab scaling
# is dropped -- where MCREML_gxg uses h_ab = std(Z_a . Z_b).  The SNPs
# themselves are still column-standardized (additive_design), exactly as
# before; it is only the pair column that changes.  Everything downstream --
# the dense cache, the Cholesky draw, mc_reml, the 4-step SLURM chain -- is
# unchanged, so a run here differs from an MCREML_gxg run in the kernel ALONE.
#
# WHAT THE CENTRED KERNEL GIVES UP.  Of the two identities of the standardized
# kernel, one goes away and nothing in this pipeline depended on it:
#     tr(W) = n     ->  now n only in expectation under independent SNPs
#     W 1 = 0       ->  STILL HOLDS: every pair column is mean-centred
# W stays PSD (a Gram matrix), so the Cholesky draw is unaffected.  In exchange
# the spectrum spreads out -- pairs are weighted by their own sigma_ab instead
# of having it normalised away -- which moves W AWAY from I.  That is the
# helpful direction: the near-collinearity of W and I is exactly what makes the
# 2-component AI matrix near-singular and h^2 weakly identified.
#
# SIMULATION builds W once (build_W_batched) and forms the Cholesky factor
# Lgxg (Lgxg Lgxg' = s2gxg W) to draw a correctly-correlated epistasis effect.
# W is DETERMINISTIC per genotype, so it is PRE-COMPUTED and cached to disk once,
# then REUSED by every replicate / variance setting.
#
# ESTIMATION (mc_reml / MC_REML) loads that pre-computed dense W and applies it
# as a plain dense mat-vec  W @ B  inside conjugate gradient (O(n^2 c) per CG
# iteration).  Reusing the cached dense W avoids the O(n m^2) matrix-free
# rebuild every CG iteration and needs no genotype at estimation time.
# See MCREML_gxg.typ for the derivation.
#
# CACHE NAMING -- READ THIS BEFORE CHANGING IT.  MCREML_gxg caches its W as
# stored_genotype/W_<mode>_n<n>_m<m>.npy, keyed by the GENOTYPE ONLY, and writes
# it only `if not os.path.exists`.  The kernel is not in that name.  This
# pipeline therefore caches to W_ctr_<mode>_n<n>_m<m>.npy (previously
# W_raw_... when the kernel was the raw product): sharing a name would mean
# that whichever pipeline ran first wins, and the second one would silently
# draw phenotypes from one kernel and estimate with the other.  Any stale
# W_raw_*.npy caches belong to the OLD raw kernel and must not be reused.
####################################################################


# ------------------------------------------------------------------ designs
def _standardize_cols(M, stability_std=1e-12):
    """Column-standardize to mean 0, variance 1 (ddof=0).

    A near-constant column (std < stability_std) is left un-scaled to avoid a
    divide-by-zero; such columns contribute ~0 to the GRM anyway.

    This standardizes the SNPs, which is unchanged from MCREML_gxg.  It is the
    1/sigma_ab scaling of the PAIR column h_ab that this variant drops -- the
    pair column keeps its mean-centring.
    """
    M = np.asarray(M, dtype=float)
    mu = M.mean(axis=0)
    sd = M.std(axis=0)
    sd = np.where(sd < stability_std, 1.0, sd)
    return (M - mu) / sd


def additive_design(real_data):
    """Additive design Z_a: column-standardized allele dosages."""
    return _standardize_cols(real_data)


# ------------------------------------------------- epistasis GRM (from onlyW)
def build_W_batched(Z, method='hadamard', pair_batch_size=5000):
    """Explicit pairwise-epistasis GRM  W = (1/p) sum_{a<b} h_ab h_ab',
    with the MEAN-CENTRED interaction column  h_ab = Z_a.*Z_b - mean(Z_a.*Z_b).

    Mean-centred but NOT variance-scaled -- only the 1/sigma_ab weight of
    MCREML_gxg's build_W_batched is dropped.  Computed ONCE per genotype and
    cached to disk; reused everywhere else.  n-by-n storage.

    method
    ------
    'hadamard' (default)
        The closed form.  Summing over UNORDERED pairs and doubling,

            2 sum_{a<b} (Z_a.*Z_b)(Z_a.*Z_b)' = sum_{a,b}(...) - sum_{a=b}(...)
                                              = (K .* K) - D D' ,

        with K = Z Z' the additive GRM and D = Z .* Z, since
        sum_{a,b} Z_ta Z_tb Z_sa Z_sb = (sum_a Z_ta Z_sa)^2 = K_ts^2.
        Centring every pair column is the projection M = I - 11'/n applied to
        each h_ab, and it factors straight through the Gram sum:

            sum_{a<b} (M h_ab)(M h_ab)' = M [ sum_{a<b} h_ab h_ab' ] M ,

        i.e. double-centre the raw closed form.  So

            W = M [ (K .* K) - D D' ] M / (2 p) ,     O(n^2 m) time,

        where the double-centring itself is O(n^2).  The standardized kernel
        has NO such shortcut -- its 1/sigma_ab weight does not factor through
        K, which is why MCREML_gxg must loop over all p = m(m-1)/2 pairs at
        O(n^2 m^2).  Dropping the scaling keeps the factor ~m saving on the
        build.
    'pairs'
        The literal pair sum, batched to bound memory, at O(n^2 p) -- exactly
        MCREML_gxg's loop with the 1/sigma_ab lines removed but the centring
        kept.  Retained as the independent check that 'hadamard' really is the
        pair sum and not merely something close to it.  Not used by the
        pipeline.

    The degenerate-pair mask of the standardized version is GONE from both
    routes, and must be: it existed only to avoid dividing by a near-zero
    sigma_ab.  With no division a near-constant product column centres to ~0
    and drops out on its own, so every pair now enters unconditionally.
    """
    n, m = Z.shape
    p = m * (m - 1) // 2
    if p == 0:
        raise ValueError("Need m >= 2 SNPs to form a pair.")

    if method == 'hadamard':
        K = Z @ Z.T
        D = Z * Z
        W = (K * K - D @ D.T) / (2.0 * p)
        W -= W.mean(axis=0, keepdims=True)   # double-centre: W <- M W M, i.e.
        W -= W.mean(axis=1, keepdims=True)   # centre every pair column h_ab
        return W

    if method != 'pairs':
        raise ValueError(f"method must be 'hadamard' or 'pairs'; got {method!r}.")

    idx_i, idx_j = np.triu_indices(m, k=1)
    W = np.zeros((n, n))
    for start in range(0, p, pair_batch_size):
        end = min(start + pair_batch_size, p)
        H = Z[:, idx_i[start:end]] * Z[:, idx_j[start:end]]
        H -= H.mean(axis=0)                                   # centred product
        W += H @ H.T
    W /= p
    return W


# ------------------------------------------------------------- simulation
def simulate_Cholesky_gxg(real_data, s2gxg=0.5, s2e=0.5, stability=1e-10):
    """Cholesky factor of the epistasis covariance, plus the GRM W itself.

    Returns (Lgxg, W, w_build_time) with
        Lgxg Lgxg' = s2gxg W .
    W (the epistasis GRM, deterministic per genotype) is returned too so the
    caller can cache it once for the estimation step's dense mat-vec.  Built
    once per (genotype, s2gxg) and reused across reps.  w_build_time is the
    wall-clock seconds spent in build_W_batched only (the kernel-build cost,
    tracked separately from the downstream estimation time).

    W is PSD but singular in general, so the `stability` jitter is what makes
    the factorisation succeed; at 1e-10 against a kernel whose diagonal is O(1)
    it is numerically invisible.
    """
    Za = additive_design(real_data)
    n, m = Za.shape

    t_start = time.perf_counter()
    W = build_W_batched(Za)
    w_build_time = time.perf_counter() - t_start

    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)
    return Lgxg, W, w_build_time


def simulate_remove_sampling_err(Lgxg, n, s2gxg=0.5, s2e=0.5):
    """Phenotype y = g_gxg + e with per-component sampling-error removal.

    g_gxg = Lgxg u1  (~ N(0, s2gxg W)),  e white noise.  Each component is
    rescaled to its exact target variance and y is mean-centred, so the fitted
    model carries no fixed effect.  Returns y (n,).

    NOTE the mean-centring is free here, exactly as in MCREML_gxg: the pair
    columns are mean-centred, so W 1 = 0 and 1 is a null direction of the
    epistasis component -- removing y's mean removes none of the signal.
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

    Identical to MCREML_gxg's mc_reml -- the centred kernel changes what W IS,
    not how it is used.

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
        # W can be nearly collinear with I, making the AI matrix near-singular
        # so that an undamped step explodes.  Levenberg-Marquardt ridge (scaled
        # to AI) + trust region on the step + a box clamp keep the path stable
        # without perturbing well-identified cases.  The un-scaled pair column
        # spreads the spectrum of W and moves it away from I, so this damping
        # should bind less often than with the standardized kernel.  It is
        # unchanged regardless, so the two pipelines remain comparable.
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
