# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.linalg import cholesky
import time

####################################################################
# Pairwise-epistasis-ONLY phenotype simulation + Monte-Carlo AI-REML.
#
# Model (no fixed effects; y is mean-centred):
#     y = g_gxg + e,
#     V = Var(y) = s2gxg * W + s2e * I ,
#     W = (1/p) sum_{a<b} h_ab h_ab' ,  h_ab = std(Z_a . Z_b),  p = m(m-1)/2
#                                (pairwise-epistasis GRM).
#
# Z_a : column-standardized allele dosages.
# W   : GRM of the standardized element-wise products of every SNP pair (the
#       epistasis relationship used by the MoM `onlyW` pipeline).
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
