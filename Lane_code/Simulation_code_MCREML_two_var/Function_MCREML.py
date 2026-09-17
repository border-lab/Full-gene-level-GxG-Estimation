# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.linalg import cholesky
import time

####################################################################
# Additive + dominance phenotype simulation + Monte-Carlo AI-REML.
#
# Model (no fixed effects; y is mean-centred):
#     y = g_a + g_d + e,
#     V = Var(y) = s2a * K_a + s2d * K_d + s2e * I,
#     K_a = Z_a Z_a' / m   (additive GRM),
#     K_d = Z_d Z_d' / m   (dominance GRM).
#
# Z_a is the column-standardized allele-dosage matrix (as in the additive-only
# pipeline).  Z_d is the column-standardized *dominance* design: the standard
# GCTA / Zhu-2015 orthogonal dominance coding maps genotype {0,1,2} to
# {-p/q, 1, -q/p} (p = allele freq, q = 1-p), which is orthogonal to the
# additive coding under HWE.  Both designs are standardized so tr(K_.)/n = 1,
# keeping the sampling-error-removal step unbiased.
#
# SIMULATION forms K_a, K_d and their Cholesky factors to draw correctly
# correlated additive / dominance effects.  ESTIMATION (mc_reml / MC_REML) is
# entirely MATRIX-FREE: no GRM is stored; each K_i enters only through the fast
# mat-vec  K_i b = Z_i (Z_i' b) / m  (O(n m)), V^{-1} is applied by conjugate
# gradient, and tr(V^{-1} K_i) is a Hutchinson stochastic estimate over
# Rademacher probes.  See MCREML_two_var.typ for the full derivation.
####################################################################


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


def dominance_design(real_data, maf_floor=1e-6):
    """Dominance design Z_d: column-standardized GCTA dominance coding.

    Genotypes are expected as allele counts in {0, 1, 2}.  With p = allele
    frequency (= column mean / 2) and q = 1 - p, the orthogonal dominance
    deviation maps
        genotype 0 -> -p/q ,   1 -> 1 ,   2 -> -q/p ,
    which is uncorrelated with the additive coding under Hardy-Weinberg.  The
    columns are then standardized (mean 0, variance 1) so that tr(K_d)/n = 1,
    mirroring the additive design.
    """
    X = np.rint(np.asarray(real_data, dtype=float))       # snap to {0,1,2}
    p = X.mean(axis=0) / 2.0                               # allele frequency
    p = np.clip(p, maf_floor, 1.0 - maf_floor)
    q = 1.0 - p

    is0 = (X == 0)
    is1 = (X == 1)
    is2 = (X == 2)
    # per-column values broadcast across rows
    W = is0 * (-p / q) + is1 * 1.0 + is2 * (-q / p)
    return _standardize_cols(W)


def simulate_Cholesky_two_var(real_data, s2a=0.5, s2d=0.25, s2e=0.25,
                              stability=1e-10):
    """Cholesky factors of the additive and dominance covariances.

    Returns (La, Ld) with
        La La' = s2a * K_a ,   K_a = Z_a Z_a' / m ,
        Ld Ld' = s2d * K_d ,   K_d = Z_d Z_d' / m ,
    so La @ u  and  Ld @ v  (u, v ~ N(0, I)) draw the additive and dominance
    effects g_a ~ N(0, s2a K_a) and g_d ~ N(0, s2d K_d).  Built once per
    (genotype, s2a, s2d) and reused across reps.
    """
    Za = additive_design(real_data)
    Zd = dominance_design(real_data)
    n, m = Za.shape

    Ka = (Za @ Za.T) / m
    Kd = (Zd @ Zd.T) / m
    La = cholesky(s2a * Ka + stability * np.eye(n), lower=True)
    Ld = cholesky(s2d * Kd + stability * np.eye(n), lower=True)
    return La, Ld


def simulate_remove_sampling_err(real_data, La, Ld, s2a=0.5, s2d=0.25, s2e=0.25):
    """Additive + dominance phenotype y = g_a + g_d + e with sampling-error removal.

    g_a = La u1  (~ N(0, s2a K_a)),  g_d = Ld u2  (~ N(0, s2d K_d)),  e white
    noise.  Each component is rescaled to its exact target variance and y is
    mean-centred, so the fitted model carries no fixed effect.

    Returns (Za, Zd, y) with Za, Zd the additive / dominance design matrices.
    """
    Za = additive_design(real_data)
    Zd = dominance_design(real_data)
    n, m = Za.shape

    u1 = np.random.randn(n)
    u2 = np.random.randn(n)
    u3 = np.random.randn(n)

    # Components
    a = La @ u1                        # additive effect  ~ N(0, s2a K_a)
    d = Ld @ u2                        # dominance effect ~ N(0, s2d K_d)
    e = np.sqrt(s2e) * u3              # residual noise

    # Eliminate sampling variances: rescale each component to its exact target.
    a *= np.sqrt(s2a / np.var(a, ddof=0))
    d *= np.sqrt(s2d / np.var(d, ddof=0))
    e *= np.sqrt(s2e / np.var(e, ddof=0))

    # Phenotype
    y = a + d + e
    y -= y.mean()

    return Za, Zd, y


####################################################################
# Matrix-free linear algebra: the V mat-vec and a batched CG solver.
####################################################################
def _v_matvec(Za, Zd, s2a, s2d, s2e, B):
    """Apply V = s2a K_a + s2d K_d + s2e I to B without forming any GRM.

        V B = (s2a / m) Z_a (Z_a' B) + (s2d / m) Z_d (Z_d' B) + s2e B .

    B may be (n,) or (n, c); the result matches its shape.  Four thin passes
    over the design matrices, so O(n m c) time and no n-by-n storage.
    """
    m = Za.shape[1]
    return ((s2a / m) * (Za @ (Za.T @ B))
            + (s2d / m) * (Zd @ (Zd.T @ B))
            + s2e * B)


def _cg_batched(matvec, Bmat, x0=None, tol=1e-6, maxiter=1000):
    """Conjugate gradient for the SPD system V X = Bmat.

    matvec : callable  X -> V X   (accepts / returns (n, c) arrays).
    Bmat   : (n, c) right-hand sides -- all c columns are solved together with
             per-column CG scalars, so one V-pass advances every column.
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
# Monte-Carlo AI-REML for V = s2a K_a + s2d K_d + s2e I.
####################################################################
def mc_reml(Za, Zd, y, iters=30, nmc=50, cg_tol=1e-6, cg_maxiter=1000,
            jitter=1e-8, tol=1e-8, seed=None, verbose=False):
    """Monte-Carlo average-information REML (additive + dominance, matrix-free).

    Each REML iteration takes the same average-information Newton step as exact
    AI-REML, but evaluates it through CG solves and a stochastic trace:
      * 1   CG solve   for  u = V^{-1} y            (data quadratics),
      * nmc CG solves  for  V^{-1} r_b   (probes)   (Hutchinson trace),
      * 3   CG solves  for  V^{-1}(K_j u)           (average information).
    Cost O(iters * (nmc + 4) * t_cg * n m); no n^2 / n^3 term, no GRM stored.

    Parameters
    ----------
    Za, Zd     : (n, m) additive / dominance design matrices.
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
    s  : (3,) estimated (s2a, s2d, s2e).
    AI : (3, 3) final average-information matrix.
    """
    y = np.asarray(y, dtype=float).flatten()
    n, m = Za.shape
    k = 3

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
        s2a, s2d, s2e = s
        matvec = lambda B: _v_matvec(Za, Zd, s2a, s2d, s2e, B)

        # --- u = V^{-1} y  (one solve) ---
        u = _cg_batched(matvec, yc, x0=u, tol=cg_tol, maxiter=cg_maxiter)
        uu = u[:, 0]

        # Data quadratics (exact given u):
        #   u'K_a u = ||Z_a'u||^2 / m,  u'K_d u = ||Z_d'u||^2 / m,  u'I u = ||u||^2
        Zatu = Za.T @ uu
        Zdtu = Zd.T @ uu
        uKau = (Zatu @ Zatu) / m
        uKdu = (Zdtu @ Zdtu) / m
        uIu = uu @ uu

        # --- Hutchinson trace: W = V^{-1} Rp  (nmc solves, shared) ---
        W = _cg_batched(matvec, Rp, x0=W, tol=cg_tol, maxiter=cg_maxiter)
        KaRp = (Za @ (Za.T @ Rp)) / m                    # K_a r_b   (n, nmc)
        KdRp = (Zd @ (Zd.T @ Rp)) / m                    # K_d r_b   (n, nmc)
        trV1Ka = np.mean(np.sum(W * KaRp, axis=0))       # ~ tr(V^{-1} K_a)
        trV1Kd = np.mean(np.sum(W * KdRp, axis=0))       # ~ tr(V^{-1} K_d)
        trV1I = np.mean(np.sum(W * Rp,  axis=0))         # ~ tr(V^{-1} I)

        # Score: 0.5 * ( u'K_i u - tr(V^{-1} K_i) )
        score = np.array([0.5 * (uKau - trV1Ka),
                          0.5 * (uKdu - trV1Kd),
                          0.5 * (uIu - trV1I)])

        # --- Average information: A_ij = 0.5 (K_i u)' V^{-1}(K_j u) ---
        Kau = (Za @ Zatu) / m                            # K_1 u = K_a u
        Kdu = (Zd @ Zdtu) / m                            # K_2 u = K_d u
        KU = np.column_stack([Kau, Kdu, uu])             # [K_1 u, K_2 u, K_3 u]
        G = _cg_batched(matvec, KU, x0=G, tol=cg_tol, maxiter=cg_maxiter)
        AI = 0.5 * (KU.T @ G)                            # (K_i u)' V^{-1}(K_j u)
        AI = 0.5 * (AI + AI.T)                           # symmetrize CG noise

        # --- Newton / Fisher-scoring update with non-negativity clamp ---
        step = np.linalg.solve(AI + jitter * np.eye(k), score)
        s = np.clip(s + step, 1e-9, None)

        if verbose:
            print(f"iter {it:2d}  s={s}  max|step|={np.abs(step).max():.3e}")
        if np.abs(step).max() < tol:
            break

    return s, AI


def MC_REML(Za, Zd, y, iters=30, nmc=50, cg_tol=1e-6, cg_maxiter=1000, seed=None):
    """Wrapper mirroring AI_REML's interface for the additive + dominance model.

    V = s2a K_a + s2d K_d + s2e I with K_a = Z_a Z_a'/m, K_d = Z_d Z_d'/m built
    implicitly.  Returns (s2a_hat, s2d_hat, s2e_hat, AI).
    """
    s, AI = mc_reml(Za, Zd, y, iters=iters, nmc=nmc, cg_tol=cg_tol,
                    cg_maxiter=cg_maxiter, seed=seed)
    s2a_hat, s2d_hat, s2e_hat = s
    return s2a_hat, s2d_hat, s2e_hat, AI
