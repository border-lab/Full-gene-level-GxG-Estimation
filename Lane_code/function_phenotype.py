# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.linalg import cholesky
import scipy.stats as stats
from scipy.sparse import random as sparse_random, eye as sparse_eye
from scipy.sparse.linalg import splu




def simulate_Cholesky_from_std(real_data, s2gxg=0.9, s2e=0.1, stability=1e-8):
    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0)
    n, m = Z.shape
    W = build_W_batched(Z)

    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)
    
    return Lgxg 



def build_W_batched(Z, batch_size=5000):
    n, m = Z.shape
    p = m * (m - 1) // 2
    
    idx_i, idx_j = np.triu_indices(m, k=1)
    
    W = np.zeros((n, n))
    
    for start in range(0, p, batch_size):
        end = min(start + batch_size, p)
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



def simulate_Cholesky_from_std_sparse(real_data, s2gxg=0.9, s2e=0.1, stability=1e-8, density=0.1):

    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0)
    n, m = Z.shape
    
    m_sub = int(m * density)

    Z_sub = Z[:, :m_sub]
    
    W = build_W_batched(Z_sub)
    
    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)
    
    return Lgxg

######pair sparse mode

def build_W_batched_sparse_front(Z, batch_size=5000, density=0.1):
    """
    Build W by selecting front (contiguous) pairs (i,j) rather than random pairs.
    
    density: fraction of SNPs (pairs scale as density^2)
    """
    n, m = Z.shape
    p_total = m * (m - 1) // 2
    
    idx_i, idx_j = np.triu_indices(m, k=1)
    
    # Calculate number of causal pairs
    m_causal = int(density * m)
    p_causal = int(m_causal * (m_causal - 1) // 2)
    p_causal = max(p_causal, 1)
    p_causal = min(p_causal, p_total)
    
    idx_i = idx_i[:p_causal]
    idx_j = idx_j[:p_causal]
    
    W = np.zeros((n, n))
    
    for start in range(0, p_causal, batch_size):
        end = min(start + batch_size, p_causal)
        H_batch = Z[:, idx_i[start:end]] * Z[:, idx_j[start:end]]
        
        # Standardize
        mu = H_batch.mean(axis=0)
        sig = H_batch.std(axis=0, ddof=0)
        mask = sig > 1e-10
        H_batch[:, mask] = (H_batch[:, mask] - mu[mask]) / sig[mask]
        H_batch[:, ~mask] = 0.0
        
        W += H_batch @ H_batch.T
    
    W /= p_causal
    return W



def simulate_Cholesky_from_std_sparse_pairVersion(real_data, s2gxg=0.9, s2e=0.1, stability=1e-8, density=0.1):

    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0)
    n, m = Z.shape
    
    W = build_W_batched_sparse(Z, batch_size=5000, density=density) 
    
    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)
        
    return Lgxg


def simulate_from_raw_only_W_remove_sampling_err(real_data, Lgxg, s2gxg=0.9, s2e=0.1):

    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0)
    n, m = Z.shape

    u1 = np.random.randn(n)
    u2 = np.random.randn(n)

    gxg = Lgxg @ u1

    Le = np.sqrt(s2e) * np.eye(n)
    e = Le @ u2 

    # eliminate the sampling variances
    cur_var_gxg = np.var(gxg, ddof=0)                 
    scale_g   = np.sqrt(s2gxg / cur_var_gxg)           
    gxg     = gxg * scale_g

    cur_var_e = np.var(e, ddof=0)                 
    scale_e   = np.sqrt(s2e / cur_var_e)           
    e     = e * scale_e
   
    y = gxg + e
    y -= y.mean()

    return  Z, y