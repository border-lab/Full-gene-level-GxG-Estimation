# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.linalg import cholesky
import time

####################################################################
# Pairwise-epistasis-ONLY phenotype simulation + Monte-Carlo AI-REML.
# MATRIX-FREE ESTIMATION variant (simulation still uses a dense W).
#
# This is the epistasis-only (k = 2) reduction of the additive+dominance+gxg
# estimator in `Simulation_code_MCREML_three_var`; the matrix-free W apply is
# the SAME `compute_WU` used there and in the MoM `onlyW` pipeline.
# See MCREML_three_var.typ for the full derivation.
#
# Model (no fixed effects; y is mean-centred):
#     y = g_gxg + e,
#     V = Var(y) = s2gxg * W + s2e * I ,
#     W = (1/p) sum_{a<b} h_ab h_ab' ,  h_ab = std(Z_a . Z_b),  p = m(m-1)/2
#                                (pairwise-epistasis GRM).
#
# Z_a : column-standardized allele dosages.
# W   : GRM of the standardized element-wise products of every SNP pair.
#
# SIMULATION (build_W_batched + simulate_Cholesky_gxg): builds the dense W ONCE
# and forms the Cholesky factor Lgxg (Lgxg Lgxg' = s2gxg W) to draw a
# correctly-correlated epistasis effect -- exactly as the dense (pre-computed-W)
# pipeline and `onlyW` do.  Using the O(n^2) W here is intentional and fine: it
# is a one-off per genotype.  W is NOT cached for estimation.
#
# ESTIMATION (mc_reml / MC_REML): fully MATRIX-FREE.  W enters only through the
# storage-free product  compute_WU(Z, U, S, R, T)  -- O(n m^2 c) per apply,
# never forming the n-by-n W nor the n-by-p interaction matrix H.  Its weight
# matrices (S, R, T) are the m-by-m contraction of the pair-standardization
# (compute_weight_matrices), built once from the genotype.  V^{-1} is applied by
# conjugate gradient and every tr(V^{-1} K_i) is a Hutchinson estimate.
#
# compute_WU(Z, U, S, R, T) == build_W_batched(Z) @ U to floating-point
# precision (verified), so this matrix-free estimator matches the dense-W
# estimator on the same phenotype; it trades arithmetic for memory so it scales
# to the large-m, memory-bound regime where a dense W will not fit.
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


# ------------------------------------------- epistasis GRM (dense, for sim)
def build_W_batched(Z, pair_batch_size=5000):
    """Explicit pairwise-epistasis GRM  W = (1/p) sum_{a<b} h_ab h_ab'.

    h_ab is the column-standardized element-wise product of SNP columns a, b.
    Built in pair-batches to bound memory; used ONLY in SIMULATION (for the
    Cholesky factor).  O(n^2 p) time, n-by-n storage.  A near-constant product
    column (std <= 1e-10) is masked to 0.  Estimation never calls this -- it
    applies W matrix-free via compute_WU.
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


# ---------------------------------- matrix-free W apply (from three_var/onlyW)
def compute_weight_matrices(Z):
    """Weight matrices (S, R, T) encoding the pair-standardization of W.

    With E[Z_a Z_b] = (Z'Z)/n and Var(Z_a Z_b) = (D'D)/n - E[Z_a Z_b]^2
    (D = Z.Z), define S = 1/V, R = E/V, T = E^2/V (diagonals zeroed, a != b).
    These let `compute_WU` apply W without forming it.  O(n m^2) once, stored
    m-by-m.

    NB: like the reference (three_var / onlyW) this does not mask degenerate
    pairs; if a SNP is constant its products have V = 0 and S = 1/V is +-inf.
    Real (MAF-filtered) genotypes have V > 0, so W entries stay finite -- match
    the simulation's build_W_batched, which masks such pairs to 0.
    """
    n, m = Z.shape
    E_ZaZb = (Z.T @ Z) / n
    D = Z * Z
    term1 = (D.T @ D) / n
    term2 = E_ZaZb * E_ZaZb
    V = term1 - term2

    S = 1.0 / V
    R = E_ZaZb / V
    T = term2 / V
    np.fill_diagonal(S, 0.0)
    np.fill_diagonal(R, 0.0)
    np.fill_diagonal(T, 0.0)
    return S, R, T


def compute_WU(Z, U, S, R, T):
    """Matrix-free product  W @ U  for the epistasis GRM.

    U : (n, c).  Returns (n, c) = W U without ever forming W or the m^2/2
    standardized interaction columns.  Cost O(c n m^2): the per-column term
    M = Z'(u . Z) is the m-by-m bottleneck.  (Identical to the three_var / onlyW
    routine; equals build_W_batched(Z) @ U to float precision.)
    """
    n, m = Z.shape
    c = U.shape[1]
    p = m * (m - 1) // 2

    ZR = Z @ R
    term2_base = 0.5 * np.sum(Z * ZR, axis=1, keepdims=True)     # (n,1)
    sum_T = 0.5 * np.sum(T)
    sum_U = np.sum(U, axis=0, keepdims=True)                     # (1,c)

    term2 = term2_base @ sum_U
    term4 = sum_T * np.ones((n, 1)) @ sum_U

    term1 = np.zeros((n, c))
    term3 = np.zeros((n, c))
    for k in range(c):
        u = U[:, k]
        M = Z.T @ (u[:, None] * Z)
        term1[:, k] = 0.5 * np.sum(Z * (Z @ (S * M)), axis=1)
        term3[:, k] = 0.5 * np.sum(R * M)

    return (term1 - term2 - term3 + term4) / p


# ------------------------------------------------------------- simulation
def simulate_Cholesky_gxg(real_data, s2gxg=0.5, s2e=0.5, stability=1e-10):
    """Cholesky factor of the epistasis covariance.

    Returns Lgxg with  Lgxg Lgxg' = s2gxg W .  Builds the dense W once (via
    build_W_batched) and factorises s2gxg W + stability I.  Built once per
    (genotype, s2gxg) and reused across reps.  W itself is not returned/cached
    because estimation applies it matrix-free from the genotype.
    """
    Za = additive_design(real_data)
    n, m = Za.shape

    W = build_W_batched(Za)
    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)
    return Lgxg


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
def _v_matvec(Z, S, R, T, s2gxg, s2e, U):
    """Apply V = s2gxg W + s2e I to U, forming no GRM.

        V U = s2gxg (W U) + s2e U ,

    with W U via compute_WU (O(n m^2 c), the bottleneck).  Nothing n-by-n or
    n-by-p is ever formed.  U may be (n, c).
    """
    return s2gxg * compute_WU(Z, U, S, R, T) + s2e * U


def _cg_batched(matvec, Bmat, x0=None, tol=1e-6, maxiter=1000):
    """Batched conjugate gradient for the SPD system V X = Bmat.

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
def mc_reml(Z, y, iters=30, Nmc=50, cg_tol=1e-6, cg_maxiter=1000,
            jitter=1e-8, tol=1e-8, lm=1e-3, step_frac=0.5, upper_mult=5.0,
            seed=None, verbose=False):
    """Matrix-free Monte-Carlo AI-REML for V = s2gxg W + s2e I.

    Epistasis-only (k = 2) reduction of the three_var estimator.  W is applied
    MATRIX-FREE via compute_WU (never formed); its weight matrices (S, R, T) are
    built once.  Per iteration: 1 CG solve for x = V^{-1} y, Nmc solves for the
    Hutchinson probes, and 2 solves for V^{-1}(K_j x).  W U is fixed and formed
    once.

    Parameters
    ----------
    Z : (n, m) column-standardized genotype (see additive_design).
    y : (n,) phenotype.

    Returns
    -------
    s  : (2,) estimated (s2gxg, s2e).
    AI : (2, 2) final average-information matrix.
    """
    y = np.asarray(y, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float)
    n, m = Z.shape
    k = 2

    vary = y.var()
    s_upper = upper_mult * vary          # no component can exceed ~total var

    # epistasis weight matrices: built once, reused every iteration
    S, R, T = compute_weight_matrices(Z)

    rng = np.random.default_rng(seed)
    U = rng.choice([-1.0, 1.0], size=(n, Nmc))      # Rademacher probes in {+-1}

    s = np.full(k, vary / k)
    AI = np.eye(k)

    yc = y.reshape(n, 1)
    xbuf = None       # warm-start buffers for the three CG solve groups
    P = None
    G = None

    WU = compute_WU(Z, U, S, R, T)       # W U (matrix-free), fixed -> once

    for it in range(iters):
        s2gxg, s2e = s
        matvec = lambda B: _v_matvec(Z, S, R, T, s2gxg, s2e, B)

        # --- x = V^{-1} y ---
        xbuf = _cg_batched(matvec, yc, x0=xbuf, tol=cg_tol, maxiter=cg_maxiter)
        x = xbuf[:, 0]

        # data quadratics x'K_i x
        Wx = compute_WU(Z, x.reshape(n, 1), S, R, T)[:, 0]
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


def MC_REML(Z, y, iters=30, Nmc=50, cg_tol=1e-6, cg_maxiter=1000, seed=None):
    """Wrapper: returns (s2gxg_hat, s2e_hat, AI).

    Z is the column-standardized genotype (additive_design output).
    """
    s, AI = mc_reml(Z, y, iters=iters, Nmc=Nmc, cg_tol=cg_tol,
                    cg_maxiter=cg_maxiter, seed=seed)
    s2gxg_hat, s2e_hat = s
    return s2gxg_hat, s2e_hat, AI
