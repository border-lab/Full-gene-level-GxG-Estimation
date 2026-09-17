# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.linalg import cholesky
import time

####################################################################
# Phenotype simulation (identical to the additive MoM pipeline so
# that AI-REML can be run on the *same* simulated phenotypes).
####################################################################
def build_W_batched(Z, pair_batch_size=5000):
    n, m = Z.shape
    p = m * (m - 1) // 2

    idx_i, idx_j = np.triu_indices(m, k=1)  # All pairs globally

    W = np.zeros((n, n))

    for start in range(0, p, pair_batch_size):
        end = min(start + pair_batch_size, p)

        # Pairwise products for this batch
        H_batch = Z[:, idx_i[start:end]] * Z[:, idx_j[start:end]]

        # Standardize
        mu = H_batch.mean(axis=0)
        sig = H_batch.std(axis=0, ddof=0)
        mask = sig > 1e-10
        H_batch[:, mask] = (H_batch[:, mask] - mu[mask]) / sig[mask]
        H_batch[:, ~mask] = 0.0

        W += H_batch @ H_batch.T

    W /= p
    return W


def build_dominance_grm(SNP):
    """Dominance GRM following Hivert et al. 2021 (Zhu et al. 2015 coding).

    SNP : (n, m) array of allele dosages, values in {0, 1, 2}.

    Each locus i (allele frequency p_i of the counted allele) is given the
    dominance coding
        xD = 0 ,  2 p_i ,  4 p_i - 2      for genotypes with 0, 1, 2 copies,
    which ensures orthogonality with the additive coding under HWE.  It is then
    standardized with the theoretical HWE moments (mean 2 p_i^2, sd 2 p_i q_i),

        wD(i) = ( xD(i) - 2 p_i^2 ) / ( 2 p_i (1 - p_i) ) ,

    and the GRM is the Hivert normalization

        Theta_D = W_D W_D' / L .

    Fitting Theta_D jointly with K and W lets AI-REML absorb the dominance-like
    signal that leaks out of the additive/interaction terms when SNPs are in LD
    (Hivert et al. 2021: the additive and dominance GRMs "should always be
    fitted when estimating additive-by-additive variance").
    """
    SNP = np.asarray(SNP, dtype=float)
    p = SNP.mean(axis=0) / 2.0          # allele frequency of the counted allele
    q = 1.0 - p

    # Zhu et al. dominance coding, broadcast over individuals (rows).
    xD = np.where(SNP < 0.5, 0.0,
         np.where(SNP < 1.5, 2.0 * p,
                             4.0 * p - 2.0))

    # Standardize with theoretical HWE moments: centre 2p^2, scale 2pq.
    sig = 2.0 * p * q
    mask = sig > 1e-10
    Wd = np.zeros_like(xD)
    Wd[:, mask] = (xD[:, mask] - 2.0 * p[mask] ** 2) / sig[mask]

    L = xD.shape[1]
    return (Wd @ Wd.T) / L


def simulate_Cholesky_from_std_withadd(real_data, s2a=0.2, s2d=0.0, s2gxg=0.7, s2e=0.1, stability=1e-10):
    """Cholesky factors of the additive, dominance and epistatic covariances.

    Returns (Lgxg, La, Ld) with
        Lgxg Lgxg' = s2gxg * W   (epistatic / GxG GRM)
        La   La'   = s2a   * K   (additive GRM)
        Ld   Ld'   = s2d   * D   (dominance GRM, Hivert/Zhu coding)
    """
    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0)
    n, m = Z.shape

    # Epistatic GRM
    W = build_W_batched(Z)
    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)

    # Additive GRM
    K = (Z @ Z.T) / m
    La = cholesky(s2a * K + stability * np.eye(n), lower=True)

    # Dominance GRM (built from the raw dosages)
    D = build_dominance_grm(real_data)
    Ld = cholesky(s2d * D + stability * np.eye(n), lower=True)

    return Lgxg, La, Ld


def simulate_remove_sampling_err(real_data, Lgxg, La, Ld, s2a=0.2, s2d=0.0, s2gxg=0.7, s2e=0.1):
    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0)
    n, m = Z.shape

    u1 = np.random.randn(n)
    u2 = np.random.randn(n)
    u3 = np.random.randn(n)
    u4 = np.random.randn(n)

    # Components
    gxg = Lgxg @ u1
    a = La @ u2
    e = np.sqrt(s2e) * u3
    d = Ld @ u4

    # Eliminate sampling variances
    cur_var_gxg = np.var(gxg, ddof=0)
    scale_gxg = np.sqrt(s2gxg / cur_var_gxg)
    gxg = gxg * scale_gxg

    cur_var_a = np.var(a, ddof=0)
    scale_a = np.sqrt(s2a / cur_var_a)
    a = a * scale_a

    cur_var_e = np.var(e, ddof=0)
    scale_e = np.sqrt(s2e / cur_var_e)
    e = e * scale_e

    # Dominance component (skip rescaling when s2d == 0 to avoid 0/0)
    if s2d > 0:
        cur_var_d = np.var(d, ddof=0)
        scale_d = np.sqrt(s2d / cur_var_d)
        d = d * scale_d
    else:
        d = np.zeros(n)

    # Phenotype
    y = gxg + a + d + e
    y -= y.mean()

    return Z, y


####################################################################
# Exact AI-REML for the general k-component model V = sum_i s_i * K_i
# (here V = s2a*K + s2d*D + s2gxg*W + s2e*I).
# No fixed effects: the phenotype is mean-centred up front (E[y] = 0),
# so there is no intercept to project out and the REML projection P
# collapses to V^{-1}.  This is ML on the mean-zero Gaussian y ~ N(0, V)
# via the average-information algorithm.  No stochastic trace estimation
# -- every trace / quadratic form is computed exactly on the dense
# matrices.
####################################################################
def reml(y, Ks, iters=12, jitter=1e-8, tol=1e-8, verbose=False):
    """Average-Information REML.

    Parameters
    ----------
    y      : (n,) phenotype vector.
    Ks     : list of (n, n) covariance-component matrices [K, W, I].
    iters  : maximum number of AI-REML iterations.
    jitter : ridge added to V and to the AI matrix for stability.
    tol    : stop early once the largest absolute update drops below tol.

    Returns
    -------
    s  : (len(Ks),) estimated variance components, same order as Ks.
    AI : (len(Ks), len(Ks)) final average-information matrix.
    """
    y = np.asarray(y, dtype=float).flatten()
    n = y.shape[0]
    k = len(Ks)
    I = np.eye(n)

    # Initialise each component at an equal share of the phenotypic variance.
    s = np.full(k, y.var() / k)

    AI = np.eye(k)
    for it in range(iters):
        V = sum(si * Ki for si, Ki in zip(s, Ks))
        Vi = np.linalg.inv(V + jitter * I)      # V^{-1}; no fixed-effect projection

        Viy = Vi @ y                            # u       = V^{-1} y
        KiViy = [Ki @ Viy for Ki in Ks]         # Ki u    = Ki V^{-1} y
        ViKiViy = [Vi @ v for v in KiViy]       # V^{-1} Ki V^{-1} y

        # Score: 0.5 * ( y' V^-1 Ki V^-1 y - tr(V^-1 Ki) )
        score = np.array([
            0.5 * (Viy @ KiViy[i] - np.sum(Vi * Ks[i]))
            for i in range(k)
        ])

        # Average information: 0.5 * (Ki V^-1 y)' V^-1 (Kj V^-1 y)
        AI = np.array([
            [0.5 * (KiViy[i] @ ViKiViy[j]) for j in range(k)]
            for i in range(k)
        ])

        step = np.linalg.solve(AI + jitter * np.eye(k), score)
        s = np.clip(s + step, 1e-9, None)

        if verbose:
            print(f"iter {it:2d}  s={s}  max|step|={np.abs(step).max():.3e}")

        if np.abs(step).max() < tol:
            break

    return s, AI


def AI_REML(Z, y, SNP=None, iters=12):
    """Exact AI-REML wrapper mirroring MoM_std's interface.

    Z is the *standardized* genotype matrix (as stored by the pipeline).
    Builds the explicit additive (K = ZZ'/m) and epistatic (W) GRMs, then
    runs exact AI-REML.

    If SNP (raw dosage matrix) is provided, a dominance GRM D is also built and
    the four-component model V = s2a*K + s2d*D + s2gxg*W + s2e*I is fitted --
    fitting D jointly keeps the additive/interaction estimates identifiable when
    SNPs are in LD (Hivert et al. 2021).  In that case the return is
    (s2a_hat, s2d_hat, s2gxg_hat, s2e_hat, AI).

    Otherwise the original three-component model V = s2a*K + s2gxg*W + s2e*I is
    fitted and the return is (s2a_hat, s2gxg_hat, s2e_hat, AI).
    """
    n, m = Z.shape

    K = (Z @ Z.T) / m
    W = build_W_batched(Z)

    if SNP is not None:
        D = build_dominance_grm(SNP)
        Ks = [K, D, W, np.eye(n)]
        s, AI = reml(y, Ks, iters=iters)
        s2a_hat, s2d_hat, s2gxg_hat, s2e_hat = s
        return s2a_hat, s2d_hat, s2gxg_hat, s2e_hat, AI

    Ks = [K, W, np.eye(n)]
    s, AI = reml(y, Ks, iters=iters)
    s2a_hat, s2gxg_hat, s2e_hat = s
    return s2a_hat, s2gxg_hat, s2e_hat, AI
