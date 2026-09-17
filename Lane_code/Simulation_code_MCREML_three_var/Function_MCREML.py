# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.linalg import cholesky
import time

####################################################################
# Additive + dominance + pairwise-epistasis phenotype simulation + MC AI-REML.
#
# Model (no fixed effects; y is mean-centred):
#     y = g_a + g_d + g_gxg + e,
#     V = Var(y) = s2a K_a + s2d K_d + s2gxg W + s2e I,
#     K_a = Z_a Z_a' / m         (additive GRM),
#     K_d = Z_d Z_d' / m         (dominance GRM),
#     W   = (1/p) sum_{a<b} h_ab h_ab' ,  h_ab = std(Z_a . Z_b),  p = m(m-1)/2
#                                (pairwise-epistasis GRM).
#
# Z_a : column-standardized allele dosages.
# Z_d : column-standardized GCTA dominance coding  {0,1,2} -> {-p/q, 1, -q/p}.
# W   : GRM of the standardized element-wise products of every SNP pair (the
#       epistasis relationship used by the MoM `onlyW` pipeline).
#
# SIMULATION forms the Cholesky factors to draw correctly correlated effects.
# ESTIMATION (mc_reml / MC_REML) is fully MATRIX-FREE: additive / dominance GRMs
# enter through thin mat-vecs  K_i b = Z_i(Z_i'b)/m, and W enters through the
# storage-free product  compute_WU  (O(n m^2 c) per apply; never forms the
# n-by-n W nor the n-by-p interaction matrix H).  V^{-1} is applied by conjugate
# gradient and every tr(V^{-1} K_i) is a Hutchinson estimate.
# See MCREML_three_var.typ for the derivation.
####################################################################


# ------------------------------------------------------------------ designs
def _standardize_cols(M, stability_std=1e-12):
    """Column-standardize to mean 0, variance 1 (ddof=0)."""
    M = np.asarray(M, dtype=float)
    mu = M.mean(axis=0)
    sd = M.std(axis=0)
    sd = np.where(sd < stability_std, 1.0, sd)
    return (M - mu) / sd


def additive_design(real_data):
    """Additive design Z_a: column-standardized allele dosages."""
    return _standardize_cols(real_data)


def dominance_design(real_data, maf_floor=1e-6):
    """Dominance design Z_d: column-standardized GCTA dominance coding.

    Genotypes expected as allele counts {0,1,2}; with p = column mean / 2,
    q = 1 - p, the orthogonal dominance deviation maps 0->-p/q, 1->1, 2->-q/p.
    """
    X = np.rint(np.asarray(real_data, dtype=float))
    p = X.mean(axis=0) / 2.0
    p = np.clip(p, maf_floor, 1.0 - maf_floor)
    q = 1.0 - p
    W = (X == 0) * (-p / q) + (X == 1) * 1.0 + (X == 2) * (-q / p)
    return _standardize_cols(W)


# ------------------------------------------------- epistasis GRM (from onlyW)
def build_W_batched(Z, pair_batch_size=5000):
    """Explicit pairwise-epistasis GRM  W = (1/p) sum_{a<b} h_ab h_ab'.

    h_ab is the column-standardized element-wise product of SNP columns a, b.
    Built in pair-batches to bound memory; used only in SIMULATION (for the
    Cholesky factor).  O(n^2 p) time, n-by-n storage.
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


def compute_weight_matrices(Z):
    """Weight matrices (S, R, T) encoding the pair-standardization of W.

    With E[Z_a Z_b] = (Z'Z)/n and Var(Z_a Z_b) = (D'D)/n - E[Z_a Z_b]^2
    (D = Z.Z), define S = 1/V, R = E/V, T = E^2/V (diagonals zeroed, a != b).
    These let `compute_WU` apply W without forming it.  O(n m^2) once.
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
    M = Z'(u . Z) is the m-by-m bottleneck.  (Identical to the onlyW routine.)
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
def simulate_Cholesky_three_var(real_data, s2a=0.4, s2d=0.2, s2gxg=0.2, s2e=0.2,
                                stability=1e-10):
    """Cholesky factors of the additive, dominance, and epistasis covariances.

    Returns (La, Ld, Lgxg, W) with
        La La'   = s2a  K_a ,
        Ld Ld'   = s2d  K_d ,
        Lgxg Lgxg' = s2gxg W .
    W itself (the epistasis GRM, deterministic per genotype) is returned too, so
    the caller can save it once for the estimation step's dense mat-vec.  Built
    once per (genotype, variance targets) and reused across reps.
    """
    Za = additive_design(real_data)
    Zd = dominance_design(real_data)
    n, m = Za.shape

    Ka = (Za @ Za.T) / m
    Kd = (Zd @ Zd.T) / m
    W = build_W_batched(Za)

    La = cholesky(s2a * Ka + stability * np.eye(n), lower=True)
    Ld = cholesky(s2d * Kd + stability * np.eye(n), lower=True)
    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)
    return La, Ld, Lgxg, W


def simulate_remove_sampling_err(real_data, La, Ld, Lgxg,
                                 s2a=0.4, s2d=0.2, s2gxg=0.2, s2e=0.2):
    """Phenotype y = g_a + g_d + g_gxg + e with per-component sampling-error removal.

    Each genetic / residual component is rescaled to its exact target variance
    and y is mean-centred.  Returns (Za, Zd, y).
    """
    Za = additive_design(real_data)
    Zd = dominance_design(real_data)
    n, m = Za.shape

    u1 = np.random.randn(n)
    u2 = np.random.randn(n)
    u3 = np.random.randn(n)
    u4 = np.random.randn(n)

    a = La @ u1
    d = Ld @ u2
    gxg = Lgxg @ u3
    e = np.sqrt(s2e) * u4

    a *= np.sqrt(s2a / np.var(a, ddof=0))
    d *= np.sqrt(s2d / np.var(d, ddof=0))
    gxg *= np.sqrt(s2gxg / np.var(gxg, ddof=0))
    e *= np.sqrt(s2e / np.var(e, ddof=0))

    y = a + d + gxg + e
    y -= y.mean()
    return Za, Zd, y


# ------------------------------------------------ matrix-free linear algebra
def _v_matvec(Za, Zd, S, R, T, s2a, s2d, s2gxg, s2e, U):
    """Apply V = s2a K_a + s2d K_d + s2gxg W + s2e I to U, forming no GRM.

        V U = (s2a/m) Z_a(Z_a'U) + (s2d/m) Z_d(Z_d'U) + s2gxg (W U) + s2e U .

    The additive / dominance terms are O(n m c); the epistasis term (W U via
    compute_WU) is O(n m^2 c) and dominates for large m.  Nothing n-by-n or
    n-by-p is ever formed.
    """
    m = Za.shape[1]
    return ((s2a / m) * (Za @ (Za.T @ U))
            + (s2d / m) * (Zd @ (Zd.T @ U))
            + s2gxg * compute_WU(Za, U, S, R, T)
            + s2e * U)


def _cg_batched(matvec, Bmat, x0=None, tol=1e-6, maxiter=1000):
    """Batched conjugate gradient for the SPD system V X = Bmat."""
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


# ----------------------------------------------------------- MC AI-REML (k=4)
def mc_reml(Za, Zd, y, iters=30, Nmc=50, cg_tol=1e-6, cg_maxiter=1000,
            jitter=1e-8, tol=1e-8, lm=1e-3, step_frac=0.5, upper_mult=5.0,
            seed=None, verbose=False):
    """Monte-Carlo AI-REML for V = s2a K_a + s2d K_d + s2gxg W + s2e I.

    Per iteration: 1 CG solve for x = V^{-1} y, Nmc solves for the Hutchinson
    probes, and 4 solves for V^{-1}(K_j x).  W is applied MATRIX-FREE via
    compute_WU (never formed); its weight matrices (S, R, T) are built once.

    Returns
    -------
    s  : (4,) estimated (s2a, s2d, s2gxg, s2e).
    AI : (4, 4) final average-information matrix.
    """
    y = np.asarray(y, dtype=float).flatten()
    n, m = Za.shape
    k = 4

    vary = y.var()
    s_upper = upper_mult * vary          # no component can exceed ~total var

    # epistasis weight matrices: built once, reused every iteration
    S, R, T = compute_weight_matrices(Za)

    rng = np.random.default_rng(seed)
    U = rng.choice([-1.0, 1.0], size=(n, Nmc))      # Rademacher probes U in {+-1}^(n x Nmc)

    s = np.full(k, vary / k)
    AI = np.eye(k)

    yc = y.reshape(n, 1)
    xbuf = None       # warm-start buffers for the three CG solve groups
    P = None
    G = None

    KaU = (Za @ (Za.T @ U)) / m         # K_a U for the fixed probes: once
    KdU = (Zd @ (Zd.T @ U)) / m         # K_d U
    WU = compute_WU(Za, U, S, R, T)     # W U (matrix-free), fixed -> once

    for it in range(iters):
        s2a, s2d, s2gxg, s2e = s
        matvec = lambda U: _v_matvec(Za, Zd, S, R, T, s2a, s2d, s2gxg, s2e, U)

        # --- x = V^{-1} y ---
        xbuf = _cg_batched(matvec, yc, x0=xbuf, tol=cg_tol, maxiter=cg_maxiter)
        x = xbuf[:, 0]

        # data quadratics x'K_i x
        Zatx = Za.T @ x
        Zdtx = Zd.T @ x
        Wx = compute_WU(Za, x.reshape(n, 1), S, R, T)[:, 0]
        xKax = (Zatx @ Zatx) / m
        xKdx = (Zdtx @ Zdtx) / m
        xWx = x @ Wx
        xIx = x @ x

        # --- Hutchinson trace: P = V^{-1} U ---
        P = _cg_batched(matvec, U, x0=P, tol=cg_tol, maxiter=cg_maxiter)
        trV1Ka = np.mean(np.sum(P * KaU, axis=0))
        trV1Kd = np.mean(np.sum(P * KdU, axis=0))
        trV1W = np.mean(np.sum(P * WU, axis=0))
        trV1I = np.mean(np.sum(P * U, axis=0))

        score = np.array([0.5 * (xKax - trV1Ka),
                          0.5 * (xKdx - trV1Kd),
                          0.5 * (xWx - trV1W),
                          0.5 * (xIx - trV1I)])

        # --- average information: A_ij = 0.5 (K_i x)' V^{-1}(K_j x) ---
        Kax = (Za @ Zatx) / m
        Kdx = (Zd @ Zdtx) / m
        KX = np.column_stack([Kax, Kdx, Wx, x])
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


def MC_REML(Za, Zd, y, iters=30, Nmc=50, cg_tol=1e-6, cg_maxiter=1000, seed=None):
    """Wrapper: returns (s2a_hat, s2d_hat, s2gxg_hat, s2e_hat, AI)."""
    s, AI = mc_reml(Za, Zd, y, iters=iters, Nmc=Nmc, cg_tol=cg_tol,
                    cg_maxiter=cg_maxiter, seed=seed)
    s2a_hat, s2d_hat, s2gxg_hat, s2e_hat = s
    return s2a_hat, s2d_hat, s2gxg_hat, s2e_hat, AI
