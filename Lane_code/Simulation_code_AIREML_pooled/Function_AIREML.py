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


def simulate_Cholesky_from_std_withadd(real_data, s2a=0.2, s2gxg=0.7, s2e=0.1, stability=1e-10):

    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0)
    n, m = Z.shape

    # Epistatic GRM
    W = build_W_batched(Z)
    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)

    K = (Z @ Z.T) / m
    La = cholesky(s2a * K + stability * np.eye(n), lower=True)

    return Lgxg, La


def simulate_remove_sampling_err(real_data, Lgxg, La, s2a=0.2, s2gxg=0.7, s2e=0.1):
    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0)
    n, m = Z.shape

    u1 = np.random.randn(n)
    u2 = np.random.randn(n)
    u3 = np.random.randn(n)

    # Components
    gxg = Lgxg @ u1
    a = La @ u2
    e = np.sqrt(s2e) * u3

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

    # Phenotype
    y = gxg + a + e
    y -= y.mean()

    return Z, y


####################################################################
# Exact AI-REML for the 3-component model:  V = s2a*K + s2gxg*W + s2e*I
# Fixed effects: intercept only (a column of ones), handled through
# the projection matrix P.  No stochastic trace estimation -- every
# trace / quadratic form is computed exactly on the dense matrices.
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
    one = np.ones(n)
    I = np.eye(n)

    # Initialise each component at an equal share of the phenotypic variance.
    s = np.full(k, y.var() / k)

    AI = np.eye(k)
    for it in range(iters):
        V = sum(si * Ki for si, Ki in zip(s, Ks))
        Vi = np.linalg.inv(V + jitter * I)

        # P = Vinv - Vinv 1 (1' Vinv 1)^-1 1' Vinv   (REML projection on the mean)
        Vio = Vi @ one
        q = one @ Vio
        P = Vi - np.outer(Vio, Vio) / q

        Py = P @ y
        KiPy = [Ki @ Py for Ki in Ks]          # Ki P y
        PKiPy = [P @ v for v in KiPy]           # P Ki P y

        # Score: 0.5 * ( y' P Ki P y - tr(P Ki) )
        score = np.array([
            0.5 * (Py @ KiPy[i] - np.sum(P * Ks[i]))
            for i in range(k)
        ])

        # Average information: 0.5 * (Ki P y)' P (Kj P y)
        AI = np.array([
            [0.5 * (KiPy[i] @ PKiPy[j]) for j in range(k)]
            for i in range(k)
        ])

        step = np.linalg.solve(AI + jitter * np.eye(k), score)
        s = np.clip(s + step, 1e-9, None)

        if verbose:
            print(f"iter {it:2d}  s={s}  max|step|={np.abs(step).max():.3e}")

        if np.abs(step).max() < tol:
            break

    return s, AI


def AI_REML(Z, y, iters=12):
    """Exact AI-REML wrapper mirroring MoM_std's interface.

    Z is the *standardized* genotype matrix (as stored by the pipeline).
    Builds the explicit additive (K = ZZ'/m) and epistatic (W) GRMs, then
    runs exact AI-REML on V = s2a*K + s2gxg*W + s2e*I.

    Returns (s2a_hat, s2gxg_hat, s2e_hat, AI) -- same ordering as MoM_std.
    """
    n, m = Z.shape

    K = (Z @ Z.T) / m
    W = build_W_batched(Z)
    Ks = [K, W, np.eye(n)]

    s, AI = reml(y, Ks, iters=iters)

    s2a_hat, s2gxg_hat, s2e_hat = s
    return s2a_hat, s2gxg_hat, s2e_hat, AI
