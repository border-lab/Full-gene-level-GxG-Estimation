# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.linalg import cholesky
import time

####################################################################
# Additive-only phenotype simulation + Monte-Carlo AI-REML.
#
# Model (no fixed effects; y is mean-centred):
#     y = g_a + e,   V = Var(y) = s2a * K + s2e * I,   K = Z Z' / m.
#
# SIMULATION forms the additive GRM K and its Cholesky factor to draw a
# correctly-correlated additive effect g_a ~ N(0, s2a K), exactly as the
# MoM / AI-REML pipelines do.  ESTIMATION (mc_reml / MC_REML) is entirely
# MATRIX-FREE: K is never formed there; it enters only through the fast
# mat-vec  K b = Z (Z' b) / m  (O(n m)), V^{-1} is applied by conjugate
# gradient, and tr(V^{-1} K_i) is a Hutchinson stochastic estimate over
# Rademacher probes.  See MCREML_additive.typ for the full derivation.
####################################################################


def simulate_Cholesky_from_std_withadd(real_data, s2a=0.5, s2e=0.5, stability=1e-10):
    """Cholesky factor of the additive covariance.

    Returns La with
        La La' = s2a * K ,   K = Z Z' / m   (additive GRM),
    so that La @ u  (u ~ N(0, I))  is a draw of the additive effect
    g_a ~ N(0, s2a K).  Built once per (genotype, s2a) and reused across reps,
    mirroring the AI-REML pipeline.
    """
    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0)
    n, m = Z.shape

    K = (Z @ Z.T) / m
    La = cholesky(s2a * K + stability * np.eye(n), lower=True)
    return La


def simulate_remove_sampling_err(real_data, La, s2a=0.5, s2e=0.5):
    """Additive-only phenotype y = g_a + e with sampling-error removal.

    g_a = La u1  (u1 ~ N(0, I))  has covariance s2a K; e is white noise.  Each
    component is rescaled to its exact target variance and y is mean-centred, so
    the fitted model carries no fixed effect.

    Returns (Z, y) with Z the standardized genotype matrix.
    """
    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0)
    n, m = Z.shape

    u1 = np.random.randn(n)
    u2 = np.random.randn(n)

    # Components
    a = La @ u1                        # additive effect ~ N(0, s2a K)
    e = np.sqrt(s2e) * u2              # residual noise

    # Eliminate sampling variances: rescale each component to its exact target.
    cur_var_a = np.var(a, ddof=0)
    scale_a = np.sqrt(s2a / cur_var_a)
    a = a * scale_a

    cur_var_e = np.var(e, ddof=0)
    scale_e = np.sqrt(s2e / cur_var_e)
    e = e * scale_e

    # Phenotype
    y = a + e
    y -= y.mean()

    return Z, y


####################################################################
# Matrix-free linear algebra: the V mat-vec and a batched CG solver.
####################################################################
def _v_matvec(Z, s2a, s2e, B):
    """Apply V = s2a K + s2e I to B without forming K or V.

        V B = (s2a / m) Z (Z' B) + s2e B .

    B may be (n,) or (n, c); the result matches its shape.  Two thin passes
    over Z, so O(n m c) time and no n-by-n storage.
    """
    m = Z.shape[1]
    return (s2a / m) * (Z @ (Z.T @ B)) + s2e * B


def _cg_batched(matvec, Bmat, x0=None, tol=1e-6, maxiter=1000):
    """Conjugate gradient for the SPD system V X = Bmat.

    matvec : callable  X -> V X   (accepts / returns (n, c) arrays).
    Bmat   : (n, c) right-hand sides -- all c columns are solved together with
             per-column CG scalars, so one Z-pass advances every column.
    x0     : (n, c) warm start (e.g. the previous REML iteration's solution).

    Each iteration is a single V mat-vec, so a solve costs O(t_cg * n m c).
    """
    n, c = Bmat.shape
    X = np.zeros((n, c)) if x0 is None else x0.copy()
    R = Bmat - matvec(X)
    P = R.copy()
    rs_old = np.sum(R * R, axis=0)                       # (c,)
    b_norm = np.sqrt(np.sum(Bmat * Bmat, axis=0))
    b_norm[b_norm == 0.0] = 1.0                          # guard zero RHS

    for _ in range(maxiter):
        VP = matvec(P)
        alpha = rs_old / np.sum(P * VP, axis=0)          # (c,)
        X += alpha * P
        R -= alpha * VP
        rs_new = np.sum(R * R, axis=0)
        if np.max(np.sqrt(rs_new) / b_norm) < tol:
            break
        beta = rs_new / rs_old
        P = R + beta * P
        rs_old = rs_new

    return X


####################################################################
# Monte-Carlo AI-REML for the additive model V = s2a K + s2e I.
####################################################################
def mc_reml(Z, y, iters=30, nmc=50, cg_tol=1e-6, cg_maxiter=1000,
            jitter=1e-8, tol=1e-8, seed=None, verbose=False):
    """Monte-Carlo average-information REML (additive-only, matrix-free).

    Each REML iteration takes the same average-information Newton step as exact
    AI-REML, but evaluates it through CG solves and a stochastic trace:
      * 1   CG solve   for  u = V^{-1} y            (data quadratics),
      * nmc CG solves  for  V^{-1} r_b   (probes)   (Hutchinson trace),
      * 2   CG solves  for  V^{-1}(K_j u)           (average information).
    Cost O(iters * (nmc + 3) * t_cg * n m); no n^2 / n^3 term, no GRM stored.

    Parameters
    ----------
    y          : (n,) phenotype.
    iters      : max REML iterations.
    nmc        : number B of Rademacher probes for the trace estimate.
    cg_tol     : relative-residual tolerance for every CG solve.
    cg_maxiter : CG iteration cap per solve.
    jitter     : ridge on the k-by-k AI matrix for the Newton solve.
    tol        : stop early once the largest |update| drops below tol.
    seed       : RNG seed for the probes (fixed -> common random numbers).

    Returns
    -------
    s  : (2,) estimated (s2a, s2e).
    AI : (2, 2) final average-information matrix.
    """
    y = np.asarray(y, dtype=float).flatten()
    n, m = Z.shape
    k = 2

    rng = np.random.default_rng(seed)
    # Rademacher probes, drawn ONCE and reused across iterations
    # (common random numbers -> smooth score, stable Newton path).
    Rp = rng.choice([-1.0, 1.0], size=(n, nmc))

    s = np.full(k, y.var() / k)          # equal-share initialization
    AI = np.eye(k)

    yc = y.reshape(n, 1)
    u = None          # warm-start buffers for the three CG solve groups
    W = None
    G = None

    for it in range(iters):
        s2a, s2e = s
        matvec = lambda B: _v_matvec(Z, s2a, s2e, B)

        # --- u = V^{-1} y  (one solve) ---
        u = _cg_batched(matvec, yc, x0=u, tol=cg_tol, maxiter=cg_maxiter)
        uu = u[:, 0]

        # Data quadratics (exact given u): u'K u = ||Z'u||^2 / m,  u'I u = ||u||^2
        Ztu = Z.T @ uu
        uKu = (Ztu @ Ztu) / m
        uIu = uu @ uu

        # --- Hutchinson trace: W = V^{-1} Rp  (nmc solves, shared) ---
        W = _cg_batched(matvec, Rp, x0=W, tol=cg_tol, maxiter=cg_maxiter)
        KRp = (Z @ (Z.T @ Rp)) / m                       # K r_b   (n, nmc)
        trV1K = np.mean(np.sum(W * KRp, axis=0))         # ~ tr(V^{-1} K)
        trV1I = np.mean(np.sum(W * Rp,  axis=0))         # ~ tr(V^{-1} I)

        # Score: 0.5 * ( u'K_i u - tr(V^{-1} K_i) )
        score = np.array([0.5 * (uKu - trV1K),
                          0.5 * (uIu - trV1I)])

        # --- Average information: A_ij = 0.5 (K_i u)' V^{-1}(K_j u) ---
        Ku = (Z @ Ztu) / m                               # K_1 u = K u
        KU = np.column_stack([Ku, uu])                   # [K_1 u, K_2 u]
        G = _cg_batched(matvec, KU, x0=G, tol=cg_tol, maxiter=cg_maxiter)
        AI = 0.5 * (KU.T @ G)                             # (K_i u)' V^{-1}(K_j u)
        AI = 0.5 * (AI + AI.T)                            # symmetrize CG noise

        # --- Newton / Fisher-scoring update with non-negativity clamp ---
        step = np.linalg.solve(AI + jitter * np.eye(k), score)
        s = np.clip(s + step, 1e-9, None)

        if verbose:
            print(f"iter {it:2d}  s={s}  max|step|={np.abs(step).max():.3e}")
        if np.abs(step).max() < tol:
            break

    return s, AI


def MC_REML(Z, y, iters=30, nmc=50, cg_tol=1e-6, cg_maxiter=1000, seed=None):
    """Wrapper mirroring AI_REML's interface for the additive-only model.

    V = s2a K + s2e I with K = ZZ'/m built implicitly.  Returns
    (s2a_hat, s2e_hat, AI).
    """
    s, AI = mc_reml(Z, y, iters=iters, nmc=nmc, cg_tol=cg_tol,
                    cg_maxiter=cg_maxiter, seed=seed)
    s2a_hat, s2e_hat = s
    return s2a_hat, s2e_hat, AI
