# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.linalg import cholesky
import time

####simulate phenotypes using Cholesky decomposition of GRMs built from standardized genotypes, with additive and epistatic components
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



def simulate_remove_sampling_err(real_data, Lgxg, La,  s2a=0.2, s2gxg=0.7,s2e=0.1):
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


def compute_weight_matrices(Z):
    n, m = Z.shape
    
    ZTZ = Z.T @ Z
    E_ZaZb = ZTZ / n
    
    D = Z * Z
    DTD = D.T @ D
    term1 = DTD / n

    term2 = E_ZaZb * E_ZaZb
    
    V = term1 - term2
    
    S = 1.0 / V
    R = E_ZaZb / V
    T = term2 / V
    
    np.fill_diagonal(S, 0)
    np.fill_diagonal(R, 0)
    np.fill_diagonal(T, 0)
    
    return S, R, T


def compute_WU(Z, U, S, R, T):

    n, m = Z.shape
    nmc = U.shape[1]
    p = m * (m - 1) // 2
    
    ZR = Z @ R
    term2_base = 0.5 * np.sum(Z * ZR, axis=1, keepdims=True)
    sum_T = 0.5 * np.sum(T)
    sum_U = np.sum(U, axis=0, keepdims=True)
    
    term2 = term2_base @ sum_U
    term4 = sum_T * np.ones((n, 1)) @ sum_U
    
    term1 = np.zeros((n, nmc))
    term3 = np.zeros((n, nmc))
    
    for k in range(nmc):
        u = U[:, k]
        M = Z.T @ (u[:, np.newaxis] * Z)
        SM = S * M
        term1[:, k] = 0.5 * np.sum(Z * (Z @ SM), axis=1)
        term3[:, k] = 0.5 * np.sum(R * M)
    
    return (term1 - term2 - term3 + term4) / p


def MoM_std(Z, y, nmc=100):
    
    n = Z.shape[0]
    m = Z.shape[1]
    p = m * (m - 1) // 2 
    S, R, T = compute_weight_matrices(Z)
    
    u = np.random.choice([-1, 1], size=(n, nmc))
    
    # Compute Wu using standardized W
    Wu = compute_WU(Z, u, S, R, T) 
    
    # Compute Ku (additive): K = ZZ'/m
    Ku = (Z @ (Z.T @ u)) / m 
    
    # --- Trace estimations ---
    # tr(W)
    trW_app = np.mean(np.sum(u * Wu, axis=0))
    
    # tr(W^2)
    trW2_app = np.mean(np.sum(Wu * Wu, axis=0))
    
    # tr(K)
    trK_app = np.mean(np.sum(u * Ku, axis=0))
    
    # tr(K^2)
    trK2_app = np.mean(np.sum(Ku * Ku, axis=0))
    
    # tr(WK)
    trWK_app = np.mean(np.sum(Wu * Ku, axis=0))
    
    # --- Quadratic forms ---
    
    yTy = y @ y
    
    y_col = y.reshape(-1, 1)
    Wy = compute_WU(Z, y_col, S, R, T).flatten()
    yTWy = y @ Wy
    
    yTKy = (y @ (Z @ (Z.T @ y))) / m
    
    
    A = np.array([
        [trK2_app, trWK_app, trK_app],
        [trWK_app, trW2_app, trW_app],
        [trK_app,  trW_app,  n]
    ])
    
    b = np.array([
        yTKy,
        yTWy,
        yTy
    ])
    
    x = np.linalg.solve(A, b)
    
    s2a_hat = x[0]
    s2gxg_hat = x[1]
    s2e_hat = x[2]
    
    return s2a_hat, s2gxg_hat, s2e_hat, A
