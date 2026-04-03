# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.linalg import cholesky

def simulate_data(n, m, s2a=0.1,s2gxg=0.2,s2e=0.7):
    Z = np.random.randn(n, m)
    Z = np.apply_along_axis(lambda x: (x - np.mean(x)) / np.std(x), 0, Z) 
    
    #form K
    K = Z @ Z.T / m
    
    ###form W
    KZ = Z @ Z.T
    D = Z*Z
    p = m*(m-1)/2
    W = 0.5*(KZ*KZ - D@D.T)/p

    y = np.random.multivariate_normal(mean=np.zeros(n), cov=s2a * K + s2gxg * W + s2e * np.eye(n))
    y -= np.mean(y)
    y /= np.std(y)
    return  Z, y

def simulate_from_raw(real_data, s2a=0.1, s2gxg=0.2, s2e=0.7):

    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0, ddof=1)
    n, m = Z.shape

    K = Z @ Z.T / m
    KZ = Z @ Z.T
    D = Z * Z
    p = m * (m - 1) / 2
    W = 0.5 * (KZ * KZ - D @ D.T) / p

    cov_y = s2a * K + s2gxg * W + s2e * np.eye(n)
    y = np.random.multivariate_normal(mean=np.zeros(n), cov=cov_y)
    y -= np.mean(y)
    y /= np.std(y)

    return Z, y

def simulate_data_sparse(n, m, density = 0.02, s2a=0.1,s2gxg=0.2,s2e=0.7):
    Z = np.random.randn(n, m)
    Z = np.apply_along_axis(lambda x: (x - np.mean(x)) / np.std(x), 0, Z) 
    S = np.random.choice(m, size=int(m * density), replace=False)
    Z_S = Z[:,S]
    M_s = len(S)

    #form K
    K = Z_S @ Z_S.T / M_s
    
    ###form W
    K_w = Z_S @ Z_S.T
    D = Z_S * Z_S
    q = M_s*(M_s-1)/2
    W = 0.5*(K_w*K_w - D@D.T)/q

    y = np.random.multivariate_normal(mean=np.zeros(n), cov=s2a * K + s2gxg * W + s2e * np.eye(n))
    y -= np.mean(y)
    y /= np.std(y)
    return  Z, y


def tr(mat):
    return np.sum(np.diag(mat))

    
def MoM(Z, y, nmc=40):
    n = Z.shape[0]
    m = Z.shape[1]
    p = m*(m-1)/2
    u = np.random.randn(n,nmc)
    ZZ = Z*Z
    
    ### three important elSements
    Ku = Z @ (Z.T @ u) / m
    Du = ZZ @ (ZZ.T @ u) 
    KKu = np.zeros((n, nmc))
    for i in range(nmc):
        Z_tilde = u[:, i].reshape(-1, 1) * Z   
        M = Z.T @ Z_tilde
        ZM = Z @ M
        KKu[:, i] = np.sum(ZM * Z, axis=1)

    
    ############Trace part
    ###trK_app
    trK_app = np.mean(np.sum(u * Ku, axis=0))
    ###trW_app
    trW_app = 0.5*(np.mean(np.sum(u * KKu, axis=0)) - np.mean(np.sum(u * Du, axis=0)))/p
    ###trK2_app
    trK2_app = np.mean(np.sum(Ku * Ku, axis=0))
    ###trW2_app
    trKK2 =np.mean(np.sum(KKu *KKu, axis=0))
    trDK =np.mean(np.sum(Du *KKu, axis=0))
    trDD = np.mean(np.sum(Du * Du, axis=0))
    trW2_app = (trKK2 - 2*trDK +trDD)/(4*p**2)
    ###trWK_app
    trWK_app = 0.5*(np.mean(np.sum(Ku * KKu, axis=0)) - np.mean(np.sum(Ku * Du, axis=0)))/p

    
    ############# Quadratic forms part
    yTy = y @ y
    yTWy = 0.5*(y.T@(np.sum((Z @(Z.T@(y.reshape(-1, 1) * Z)))*Z,axis = 1)) - (y @ (ZZ @ (ZZ.T @ y))))/p
    yTKy =(y @ (Z @ (Z.T @ y))) / m

 

    A = np.array([
        [trK2_app, trWK_app, trK_app],
        [trWK_app, trW2_app, trW_app],
        [trK_app,  trW_app,  n]
    ])
    
    
    b = np.array([
        yTKy,
        yTWy,
        yTy
    ]).reshape(-1, 1)
     
    x1 = np.linalg.solve(A, b)
    return x1[0],x1[1],x1[2],A


def SE(Z,y, s2a_hat,s2gxg_hat,s2e_hat,T):
    n = Z.shape[0]
    m = Z.shape[1]
    ZZ = Z*Z
    p = m*(m-1)/2
    
    def k_middle(left,right):
        return float(np.sum((left@Z)* (right@Z))/m)
        
    def w_middle(left,right):
        return float((0.5*(left@(np.sum((Z @(Z.T@(right.reshape(-1, 1) * Z)))*Z,axis = 1)) - (left @ (ZZ @ (ZZ.T @ right))))/p))
        
    def no_middle(left,right):
        return float(left@right)

    Ky = Z @ (Z.T @ y) / m
    Wy = (0.5*((np.sum((Z @(Z.T@(y.reshape(-1, 1) * Z)))*Z,axis = 1)) - ((ZZ @ (ZZ.T @ y))))/p)
    side = [Ky,Wy,y]
    middle = [k_middle,w_middle,no_middle]
    var_component= [s2a_hat,s2gxg_hat,s2e_hat]
    cov_q = np.zeros((3,3))
    for i in range(3):
         for j in range(i, 3):
            val = 0
            for l in range(3):
                val += 2*(var_component[l].item() * middle[l](side[i], side[j]))
            cov_q[i, j] = val
            cov_q[j,i] = val
             
    X = np.linalg.solve(T, cov_q)
    Cov_theta = np.linalg.solve(T, X.T).T  
    return Cov_theta[0,0],Cov_theta[1,1],Cov_theta[2,2]




########################################################################################################################################################real data without gxg
def simulate_from_raw_add(real_data, s2a=0.3,s2e=0.7):

    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0, ddof=1)
    n, m = Z.shape
    K = Z @ Z.T / m

    cov_y = s2a * K + s2e * np.eye(n)
    y = np.random.multivariate_normal(mean=np.zeros(n), cov=cov_y)
    y -= np.mean(y)
    y /= np.std(y)

    return Z, y

def MoM_only_add(Z, y, nmc=80):
    n = Z.shape[0]
    m = Z.shape[1]
    u = np.random.randn(n,nmc)

    ### three important elSements
    Ku = Z @ (Z.T @ u) / m
    ############Trace part
    ###trK_app
    trK_app = np.linalg.norm(Z,ord='fro')**2/m
    ###trK2_app
    trK2_app = np.mean(np.sum(Ku * Ku, axis=0))

    ############# Quadratic forms part
    yTy = y @ y
    yTKy =(y @ (Z @ (Z.T @ y))) / m

    A = np.array([
        [trK2_app, trK_app],
        [trK_app, n],
    ])
    
    b = np.array([
        yTKy,
        yTy
    ]).reshape(-1, 1)
     
    x1 = np.linalg.solve(A, b)
    return x1[0],x1[1]


###############################################scaler

def var_XiXj(X):
    n, m = X.shape
    R = np.corrcoef(X, rowvar=False)
    R2 = R**2
    X2 = X**2
    E2 = (X2.T @ X2) / n      
    Var_mat = E2 - R2
    S = (np.sum(Var_mat) - np.sum(np.diag(Var_mat))) / 2
    return S

def var_XiXj_columns(X):
    n_cols = X.shape[1]
    vars_ij = []
    for i in range(n_cols):
        for j in range(i+1, n_cols):
            prod = X[:, i] * X[:, j]
            vars_ij.append(np.var(prod, ddof=0))  # population variance
    return  np.array(vars_ij)

def simulate_from_raw_only_W(real_data, s2gxg=0.2, s2e=0.8):

    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0, ddof=1)
    n, m = Z.shape

    KZ = Z @ Z.T
    D = Z * Z
    p = m * (m - 1) / 2
    W = 0.5 * (KZ * KZ - D @ D.T) / p

    cov_y = s2gxg * W + s2e * np.eye(n)
    y = np.random.multivariate_normal(mean=np.zeros(n), cov=cov_y)
    y -= np.mean(y)
    y /= np.std(y)

    return Z, y


def simulate_from_raw_only_W_scale(real_data, s2gxg=0.2, s2e=0.8):

    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0, ddof=1)
    n, m = Z.shape
    p = m*(m-1)/2
    origin_var =var_XiXj(Z)
    s = (origin_var/p) ** 0.25
    Z = Z/s

    KZ = Z @ Z.T
    D = Z * Z
    p = m * (m - 1) / 2
    W = 0.5 * (KZ * KZ - D @ D.T) / p

    gxg = np.random.multivariate_normal(mean=np.zeros(n), cov= s2gxg * W)
    e = np.random.multivariate_normal(mean=np.zeros(n), cov= s2e * np.eye(n))
    y = gxg + e
    y -= np.mean(y)
    y /= np.std(y)

    return Z, y



def MoM_only_M(Z, y, nmc=40):
    n = Z.shape[0]
    m = Z.shape[1]
    p = m*(m-1)/2
    u = u = np.random.choice([-1, 1], size=(n, nmc))
    ZZ = Z*Z
    
    Du = ZZ @ (ZZ.T @ u) 

    KKu = np.zeros((n, nmc))
    for i in range(nmc):
        Z_tilde = u[:, i].reshape(-1, 1) * Z   
        M = Z.T @ Z_tilde
        ZM = Z @ M
        KKu[:, i] = np.sum(ZM * Z, axis=1)

    ###trW_app
    trW_app = 0.5*(np.mean(np.sum(u * KKu, axis=0)) - np.mean(np.sum(u * Du, axis=0)))/p
    trKK2 =np.mean(np.sum(KKu *KKu, axis=0))
    trDK =np.mean(np.sum(Du *KKu, axis=0))
    trDD = np.mean(np.sum(Du * Du, axis=0))
    trW2_app = (trKK2 - 2*trDK +trDD)/(4*p**2)
    ############# Quadratic forms part
    yTy = y @ y
    yTWy = 0.5*(y.T@(np.sum((Z @(Z.T@(y.reshape(-1, 1) * Z)))*Z,axis = 1)) - (y @ (ZZ @ (ZZ.T @ y))))/p

    A = np.array([
        [trW2_app, trW_app],
        [trW_app, n]
    ])
    
    b = np.array([
        yTWy,
        yTy
    ]).reshape(-1, 1)
     
    x1 = np.linalg.solve(A, b)
    return x1[0],x1[1],A



def ld_prune(Z, window=50, step=5, threshold=0.1, verbose=True):
    """
    Perform LD pruning using sliding window (r^2 threshold).
    Returns:
        Z_pruned: matrix after pruning
        keep: list of retained SNP column indices
    """
    M = Z.shape[1]
    keep = list(range(M))
    i = 0

    while i < M:
        start = i
        end = min(i + window, M)

        R_window = np.corrcoef(Z[:, start:end], rowvar=False)**2

        for j in range(start, end):
            if j not in keep:
                continue
            for k in range(j + 1, end):
                if k in keep and R_window[j - start, k - start] >= threshold:
                    keep.remove(k)

        i += step

    Z_pruned = Z[:, keep]
    return Z_pruned


###################################################################METHOD From W directly
def simulate_Cholesky_from_raw(real_data,s2gxg=0.9,s2e=0.1,stability=1e-8):
    """
    Compute the GxG interaction kernel W from raw genotype data.
    Parameters:
    real_data: (n, m) array of raw genotype data (n samples, m SNPs)
    Returns:
    W: (n, n) GxG interaction kernel matrix
    """

    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0, ddof=1)
    n,m = Z.shape

    KZ = Z @ Z.T
    D = Z*Z
    p = m*(m-1)/2
    W = 0.5*(KZ*KZ - D@D.T)/p
    Lgxg = cholesky(s2gxg * W + stability* np.eye(n), lower=True)
    Le = cholesky(s2e * np.eye(n), lower=True)
    return  Lgxg, Le 

def simulate_Cholesky_from_raw_correction(real_data, s2gxg=0.9, s2e=0.1, stability=1e-8):
    """
    Compute the GxG interaction kernel W from raw genotype data with LD correction.
    
    Parameters:
    -----------
    real_data: (n, m) array of raw genotype data (n samples, m SNPs)
    s2gxg: GxG variance component
    s2e: residual variance component
    stability: numerical stability term
    
    Returns:
    --------
    Lgxg: Cholesky factor of the GxG covariance component
    Le: Cholesky factor of the residual covariance component
    Z_corr: corrected standardized genotype matrix
    """
    # Calculate correction factor
    factor = calculate_correction_factor(real_data)
    
    # Standardize
    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0, ddof=1)
    Z_corr = Z / (factor ** 0.25)
    
    n, m = Z_corr.shape
    
    # Compute kernel using CORRECTED Z
    KZ = Z_corr @ Z_corr.T
    D = Z_corr * Z_corr
    p = m * (m - 1) / 2
    W = 0.5 * (KZ * KZ - D @ D.T) / p
    
    # Cholesky decomposition
    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)
    Le = cholesky(s2e * np.eye(n), lower=True)
    
    return Lgxg, Le




def simulate_from_raw_only_W(real_data, s2gxg=0.9, s2e=0.1):
    """
    Simulate phenotype y from real genotype data using only the GxG interaction kernel W.
    
    Parameters:
    real_data: (n, m) array of raw genotype data (n samples, m SNPs)
    s2gxg: variance component for GxG epistatic interactions (default 0.9)
    s2e: residual/environmental variance component (default 0.1)
    
    Returns:
    Z: (n, m) standardized genotype matrix
    y: (n,) simulated phenotype vector (centered)
    """
    # Standardize genotype data (zero mean, unit variance per SNP)
    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0, ddof=1)
    n, m = Z.shape
    
    # Compute linear kernel (genetic relationship matrix)
    KZ = Z @ Z.T
    
    # Element-wise squared genotypes for diagonal correction
    D = Z * Z
    
    # Number of unique SNP pairs
    p = m * (m - 1) / 2
    
    # GxG interaction kernel: W = 0.5 * (KZ ¡Ñ KZ - D @ D.T) / p
    W = 0.5*(KZ *KZ -D@D.T)/p
    
    # Phenotype covariance: Var(y) = s2gxg * W + s2e * I
    cov_y = s2gxg * W + s2e * np.eye(n)
    
    # Sample phenotype from multivariate normal
    y = np.random.multivariate_normal(mean=np.zeros(n), cov=cov_y)
    y -= np.mean(y)  # Center the phenotype
    
    return Z, y



def MoM_only_M_standardized(Z, y, nmc=40):
    """
    Method of Moments estimator with standardized GxG kernel.
    W is scaled so that tr(W) = n for unbiased estimation.
    """
    n = Z.shape[0]
    m = Z.shape[1]
    p = m * (m - 1) / 2 
    
    # Random vectors for stochastic trace estimation
    u = np.random.choice([-1, 1], size=(n, nmc))
    
    # D = Z ¡Ñ Z (element-wise square)
    ZZ = Z * Z
    
    # Compute D @ D.T @ u
    Du = ZZ @ (ZZ.T @ u)

    # Compute (K ¡Ñ K) @ u where K = Z @ Z.T
    KKu = np.zeros((n, nmc))
    for i in range(nmc):
        Z_tilde = u[:, i].reshape(-1, 1) * Z
        M = Z.T @ Z_tilde
        ZM = Z @ M
        KKu[:, i] = np.sum(ZM * Z, axis=1)

    # W_raw = 0.5 * (K¡ÑK - D@D.T) / p  (before standardization)
    trW_raw = 0.5 * (np.mean(np.sum(u * KKu, axis=0)) - np.mean(np.sum(u * Du, axis=0))) / p
    
    # Scaling factor to make tr(W) = n
    # W_std = W_raw * (n / tr(W_raw))
    scale = n / trW_raw
    
    # Apply scaling to all W-related terms
    # tr(W_std) = n (by construction)
    trW_app = n
    
    # tr(W_std^2) = tr(W_raw^2) * scale^2
    trKK2 = np.mean(np.sum(KKu * KKu, axis=0))
    trDK = np.mean(np.sum(Du * KKu, axis=0))
    trDD = np.mean(np.sum(Du * Du, axis=0))
    trW2_raw = (trKK2 - 2 * trDK + trDD) / (4 * p**2)
    trW2_app = trW2_raw * (scale ** 2)
    
    # Quadratic forms
    yTy = y @ y
    
    # y.T @ W_raw @ y
    yTWy_raw = 0.5 * (y.T @ (np.sum((Z @ (Z.T @ (y.reshape(-1, 1) * Z))) * Z, axis=1)) - (y @ (ZZ @ (ZZ.T @ y)))) / p
    
    # y.T @ W_std @ y = y.T @ W_raw @ y * scale
    yTWy = yTWy_raw * scale

    # ============ Solve moment equations ============
    # E[y.T @ W @ y] = s2gxg * tr(W^2) + s2e * tr(W)
    # E[y.T @ y] = s2gxg * tr(W) + s2e * n
    A = np.array([
        [trW2_app, trW_app],
        [trW_app, n]
    ])
    
    b = np.array([
        yTWy,
        yTy
    ]).reshape(-1, 1)
    
    x1 = np.linalg.solve(A, b)
    
    return x1[0], x1[1], A


###== new sparse

def simulate_Cholesky_from_raw_sparse(real_data, s2gxg=0.9, s2e=0.1, causal_fraction=0.25,stability=1e-8):
    """
    Compute the GxG interaction kernel W from raw genotype data using SPARSE causal SNPs.
    Parameters:
    -----------
    real_data: (n, m) array of raw genotype data (n samples, m SNPs)
    s2gxg: variance component for GxG epistatic interactions
    s2e: residual/environmental variance component
    causal_fraction: fraction of SNPs that are causal (default 0.25)
    stability: small value for numerical stability
    seed: random seed for reproducibility
    
    Returns:
    --------
    Lgxg: Cholesky factor of the GxG covariance component (from causal SNPs)
    Le: Cholesky factor of the residual covariance component
    causal_indices: indices of causal SNPs
    """
    
    # Standardize ALL genotype data
    Z_full = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0)
    n, m = Z_full.shape
    
    # Select causal SNPs
    m_causal = int(m * causal_fraction)
    causal_indices = np.random.choice(m, size=m_causal, replace=False)
    causal_indices = np.sort(causal_indices)
    
    # Extract causal SNPs
    Z_causal = Z_full[:, causal_indices]
    
    # Compute kernel using ONLY causal SNPs
    KZ_causal = Z_causal @ Z_causal.T
    D_causal = Z_causal * Z_causal
    p_causal = m_causal * (m_causal - 1) / 2
    
    # GxG interaction kernel from causal SNPs only
    W_causal = 0.5 * (KZ_causal * KZ_causal - D_causal @ D_causal.T) / p_causal
    
    # Cholesky decomposition
    Lgxg = cholesky(s2gxg * W_causal + stability * np.eye(n), lower=True)
    Le = cholesky(s2e * np.eye(n), lower=True)
    
    return Lgxg, Le


### ======== LD correction
def calculate_correction_factor(X):
    """
    Calculate the correction factor E[Var(Zi*Zj)] from genotype matrix X.
    
    Returns:
    --------
    correction_factor: mean of Var(Zi*Zj) for all pairs
    """
    n, m = X.shape
    # Allele frequency and covariance
    p = np.sum(X, axis=0) / (2 * n)
    sigma_matrix = np.cov(X, rowvar=False)
    
    # Broadcasting
    pi = p[:, np.newaxis]
    pj = p[np.newaxis, :]
    
    # Var(Zi*Zj) matrix - CORRECT FORMULA (no sqrt)
    numerator = sigma_matrix * (1 - 2*pi) * (1 - 2*pj)
    denominator = 4 * pi * pj * (1 - pi) * (1 - pj)  
    
    var_mtx = np.where(denominator > 1e-10, 
                       1 + numerator / denominator, 
                       1)

    # Mean of upper triangle
    rows, cols = np.triu_indices(m, k=1)
    correction_factor = np.mean(var_mtx[rows, cols])
    
    return correction_factor

def calculate_per_snp_correction(real_data):
    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0, ddof=1)

    n, m = Z.shape
    p = real_data.mean(axis=0) / 2
    R = np.corrcoef(Z.T)
    v = (1 - 2 * p)
    V = np.outer(v, v)
    d = 4 * np.sqrt(np.outer(p * (1-p), p * (1-p)))
    delta = R * V / d
    VarZiZj = 1 + delta  # (m, m)
    # Per-SNP mean variance
    per_snp_factor = VarZiZj.mean(axis=1)  # (m,)
    return per_snp_factor


def MoM_only_M_per_snp_corrected(real_data, y, nmc=40):
    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0, ddof=1)
    per_snp_factor = calculate_per_snp_correction(real_data)
    # Correct each SNP individually
    Z_corr = Z / (per_snp_factor ** 0.25)
    return MoM_only_M(Z_corr, y, nmc=nmc)

def empirical_var_zizj(Z):
    """
    Compute empirical Var(Z_i Z_j) for all pairs (i, j).
    
    Parameters:
    Z: (n, m) standardized genotype matrix
    
    Returns:
    VarZiZj: (m, m) matrix of empirical pairwise product variances
    """
    n, m = Z.shape
    # Outer product of columns: (n, m, m)
    # VarZiZj[i,j] = Var(Z[:,i] * Z[:,j])
    VarZiZj = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            prod = Z[:, i] * Z[:, j]
            VarZiZj[i, j] = np.var(prod, ddof=1)
    return VarZiZj

def MoM_implicit_corrected(Z, y, nmc=40):
    """
    MoM estimator with pairwise product standardization,
    computed implicitly without building W explicitly.
    
    For each pair (i,j):
        tilde_zij = (Z[:,i] * Z[:,j] - mean) / std
    
    W u = (1/p) * sum_{i<j} tilde_zij * (tilde_zij^T u)
    """
    n, m = Z.shape
    p = m * (m - 1) // 2
    
    u = np.random.choice([-1, 1], size=(n, nmc))  # (n, nmc)
    
    # Accumulators
    Wu        = np.zeros((n, nmc))   # W @ u
    Wy        = np.zeros(n)          # W @ y
    trW_acc   = 0.0                  # tr(W)
    trW2_acc  = 0.0                  # tr(W^2) = ||W||_F^2 approx

    for i in range(m):
        for j in range(i + 1, m):
            # Pairwise product (n,)
            zij = Z[:, i] * Z[:, j]
            
            # Standardize to zero mean, unit variance
            mu  = zij.mean()
            sig = zij.std()
            if sig < 1e-10:
                continue           # skip constant pairs
            zij = (zij - mu) / sig
            
            # zij^T u : (nmc,)
            ztu = zij @ u          # (nmc,)
            
            # Accumulate W @ u
            Wu += np.outer(zij, ztu)          # (n, nmc)
            
            # Accumulate W @ y
            Wy += zij * (zij @ y)             # (n,)
            
            # tr(W) += zij^T zij / p^2 (after dividing by p)
            # but we divide by p at the end, so:
            trW_acc  += (zij @ zij)           # = n (since unit variance)
            
            # tr(W^2) approximation via Hutchinson:
            # tr(W^2) = E[u^T W^2 u] = E[||Wu||^2] / nmc
            # accumulated after loop

    # Divide by p
    Wu   /= p
    Wy   /= p
    trW_acc /= p ** 2

    # Trace estimates
    trW_app  = trW_acc                        # = n/p ¡Ö constant
    trW2_app = np.mean(np.sum(Wu * Wu, axis=0))  # Hutchinson tr(W^2)

    # Quadratic forms
    yTy  = float(y @ y)
    yTWy = float(y @ Wy)

    # MoM system
    A = np.array([
        [trW2_app, trW_app],
        [trW_app,  n]
    ])
    b = np.array([yTWy, yTy])

    x = np.linalg.solve(A, b)
    return x[0], x[1], A

##=== explicit way


def simulate_from_raw_only_standard_W_remove_sampling_err(real_data, s2gxg=0.9, s2e=0.1, stability=1e-8):
    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0, ddof=1)
    n, m = Z.shape
    p = m * (m - 1) // 2

    # Build standardized H and W
    idx_i, idx_j = np.triu_indices(m, k=1)
    H = Z[:, idx_i] * Z[:, idx_j]
    mu  = H.mean(axis=0)
    sig = H.std(axis=0, ddof=1)
    mask = sig > 1e-10
    H[:, mask]  = (H[:, mask] - mu[mask]) / sig[mask]
    H[:, ~mask] = 0.0
    W = (H @ H.T) / p

    # Option 2: Cholesky (efficient - O(n^3) once, O(n^2) per sample)
    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)
    gxg  = Lgxg @ np.random.randn(n)              # ~ N(0, s2gxg * W)
    e    = np.sqrt(s2e) * np.random.randn(n)      # ~ N(0, s2e * I)

    # Eliminate sampling variances
    gxg = gxg * np.sqrt(s2gxg / np.var(gxg, ddof=0))
    e   = e   * np.sqrt(s2e   / np.var(e,   ddof=0))

    y = gxg + e
    y -= y.mean()
    return Z, y,W



def MoM_useM(W, y, nmc=40):

    n = W.shape[0]
    
    # Hutchinson trace estimator
    u = np.random.choice([-1, 1], size=(n, nmc))   # (n, nmc)
    Wu = W @ u                                      # (n, nmc)
    
    # Trace estimates
    trW_app  = float(np.mean(np.sum(u * Wu, axis=0)))        # tr(W)
    trW2_app = float(np.mean(np.sum(Wu * Wu, axis=0)))       # tr(W^2)
    
    # Quadratic forms
    yTy  = float(y @ y)
    yTWy = float(y @ (W @ y))
    
    # MoM system
    A = np.array([
        [trW2_app, trW_app],
        [trW_app,  n      ]
    ])
    b = np.array([yTWy, yTy])
    
    x = np.linalg.solve(A, b)
    return x[0], x[1], A

###############################################################################################PSC mode(updated)
import numpy as np
import pandas as pd
from scipy.linalg import cholesky
import scipy.stats as stats
from scipy.sparse import random as sparse_random, eye as sparse_eye
from scipy.sparse.linalg import splu


def simulate_Cholesky_from_raw(real_data,s2gxg=0.9,s2e=0.1,stability=1e-8):
    """
    Compute the GxG interaction kernel W from raw genotype data.
    Parameters:
    real_data: (n, m) array of raw genotype data (n samples, m SNPs)
    Returns:
    Lgxg: Cholesky factor of the GxG covariance component   
    Le: Cholesky factor of the residual covariance component
    """

    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0)
    n,m = Z.shape

    KZ = Z @ Z.T
    D = Z*Z
    p = m*(m-1)/2
    W = 0.5*(KZ*KZ - D@D.T)/p
    Lgxg = cholesky(s2gxg * W + stability* np.eye(n), lower=True)
    Le = cholesky(s2e * np.eye(n), lower=True)
    return  Lgxg, Le 

def simulate_from_raw_only_W_remove_sampling_err(real_data, Lgxg, Le, s2gxg=0.9, s2e=0.1):
    """
    Simulate phenotype y from real genotype data using precomputed GxG interaction kernel W and removing sampling errors.
    Parameters:
    real_data: (n, m) array of raw genotype data (n samples, m SNPs)
    Lgxg: Cholesky factor of the GxG covariance component   
    Le: Cholesky factor of the residual covariance component
    s2gxg: variance component for GxG epistatic interactions (default 0.9)
    s2e: residual/environmental variance component (default 0.1)
    Returns:
    Z: (n, m) standardized genotype matrix
    y: (n,) simulated phenotype vector (centered)
    """

    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0)
    n, m = Z.shape

    u1 = np.random.randn(n)
    u2 = np.random.randn(n)

    gxg = Lgxg @ u1
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

def MoM_only_M(Z, y, nmc=40):
    """
    Method of Moments estimator for GxG variance component using only the epistatic kernel.

    Parameters:
    Z: (n, m) standardized genotype matrix
    y: (n,) phenotype vector
     nmc: number of Monte Carlo samples for trace estimation
    
    Returns:
    s2gxg: estimated GxG variance component
    s2e: estimated residual variance component
    A: 2x2 coefficient matrix from the moment equations
    """
    n = Z.shape[0]
    m = Z.shape[1]
    p = m * (m - 1) / 2 
    
    # Random vectors for stochastic trace estimation
    u = np.random.choice([-1, 1], size=(n, nmc))
    
    # D = Z ¡Ñ Z (element-wise square)
    ZZ = Z * Z
    
    # Compute D @ D.T @ u for trace estimation
    Du = ZZ @ (ZZ.T @ u)

    # Compute (K ¡Ñ K) @ u where K = Z @ Z.T
    # Done without forming K explicitly to save memory
    KKu = np.zeros((n, nmc))
    for i in range(nmc):
        # Z_tilde = diag(u_i) @ Z
        Z_tilde = u[:, i].reshape(-1, 1) * Z
        # M = Z.T @ Z_tilde = Z.T @ diag(u_i) @ Z
        M = Z.T @ Z_tilde
        # ZM @ Z.T gives (K ¡Ñ K) @ u_i via the identity
        ZM = Z @ M
        KKu[:, i] = np.sum(ZM * Z, axis=1)

    # tr(W) ¡Ö u.T @ W @ u, where W = 0.5 * (K¡ÑK - D@D.T) / p
    trW_app = 0.5 * (np.mean(np.sum(u * KKu, axis=0)) - np.mean(np.sum(u * Du, axis=0))) / p
    
    # tr(W^2) estimation via tr((K¡ÑK)^2), tr((K¡ÑK)(D@D.T)), tr((D@D.T)^2)
    trKK2 = np.mean(np.sum(KKu * KKu, axis=0))
    trDK = np.mean(np.sum(Du * KKu, axis=0))
    trDD = np.mean(np.sum(Du * Du, axis=0))
    trW2_app = (trKK2 - 2 * trDK + trDD) / (4 * p**2)
    
    #Quadratic forms
    yTy = y @ y  # y.T @ y
    
    # y.T @ W @ y = 0.5 * (y.T @ (K¡ÑK) @ y - y.T @ D @ D.T @ y) / p
    yTWy = 0.5 * (y.T @ (np.sum((Z @ (Z.T @ (y.reshape(-1, 1) * Z))) * Z, axis=1)) - (y @ (ZZ @ (ZZ.T @ y)))) / p

    # solve moment equations
    # E[y.T @ W @ y] = s2gxg * tr(W^2) + s2e * tr(W)
    # E[y.T @ y] = s2gxg * tr(W) + s2e * n
    A = np.array([
        [trW2_app, trW_app],
        [trW_app, n]
    ])
    
    b = np.array([
        yTWy,
        yTy
    ]).reshape(-1, 1)
    
    # Solve A @ [s2gxg, s2e].T = b
    x1 = np.linalg.solve(A, b)
    
    return x1[0], x1[1], A 

###############STD method
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

def build_W(Z):

    n, m = Z.shape
    p = m * (m - 1) // 2
    
    # Get indices for upper triangle (i < j)
    idx_i, idx_j = np.triu_indices(m, k=1)
    
    # Hadamard product of all pairs: H matrix (n x p)
    H = Z[:, idx_i] * Z[:, idx_j]
    
    # Standardize each column of H
    mu = H.mean(axis=0)
    sig = H.std(axis=0, ddof=0)
    
    mask = sig > 1e-10
    H[:, mask] = (H[:, mask] - mu[mask]) / sig[mask]
    H[:, ~mask] = 0.0
    
    # W = p^{-1} H H^T
    W = (H @ H.T) / p
    
    return W


def simulate_Cholesky_from_std(real_data, s2gxg=0.9, s2e=0.1, stability=1e-8):
    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0)
    n, m = Z.shape
    W = build_W_batched(Z)
    
    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)
    Le = np.sqrt(s2e) * np.eye(n)  
    
    return Lgxg, Le

def compute_weight_matrices(Z):
    n, m = Z.shape
    
    ZTZ = Z.T @ Z
    E_ZaZb = ZTZ / n
    
    D = Z * Z
    DTD = D.T @ D
    E_ZaZb_sq = DTD / n
    
    V = E_ZaZb_sq - E_ZaZb ** 2
    
    S = 1.0 / V
    R = E_ZaZb / V
    T = E_ZaZb ** 2 / V
    
    np.fill_diagonal(S, 0)
    np.fill_diagonal(R, 0)
    np.fill_diagonal(T, 0)
    
    return S, R, T


def compute_WU(Z, U, S, R, T):

    n, m = Z.shape
    nmc = U.shape[1]
    p = m * (m - 1) // 2
    
    # Precompute shared terms
    ZR = Z @ R
    term2_base = 0.5 * np.sum(Z * ZR, axis=1, keepdims=True)
    sum_T = 0.5 * np.sum(T)
    sum_U = np.sum(U, axis=0, keepdims=True)
    
    # Term 2 & 4 
    term2 = term2_base @ sum_U
    term4 = sum_T * np.ones((n, 1)) @ sum_U
    
    # Term 1 & 3 
    term1 = np.zeros((n, nmc))
    term3 = np.zeros((n, nmc))
    
    for k in range(nmc):
        u = U[:, k]
        M = Z.T @ (u[:, np.newaxis] * Z)
        SM = S * M
        term1[:, k] = 0.5 * np.sum(Z * (Z @ SM), axis=1)
        term3[:, k] = 0.5 * np.sum(R * M)
    
    return (term1 - term2 - term3 + term4) / p


def MoM_only_M_std(Z, y, nmc=40):
    n = Z.shape[0]
    m = Z.shape[1]
    p = m * (m - 1) / 2 
    
    S, R, T = compute_weight_matrices(Z)
    
    # Random vectors for stochastic trace estimation
    u = np.random.choice([-1, 1], size=(n, nmc))
    
    # Compute Wu using standardized W
    Wu = compute_WU(Z, u, S, R, T)  # This should return W @ u
    
    # tr(W) estimation: tr(W) = mean(u.T @ W @ u)
    trW_app = np.mean(np.sum(u * Wu, axis=0))
    
    # tr(W^2) estimation: tr(W^2) = mean(||W @ u||^2)
    trW2_app = np.mean(np.sum(Wu * Wu, axis=0))
    
    # Quadratic forms
    yTy = y @ y
    
    # y.T @ W @ y - need to compute W @ y first
    y_col = y.reshape(-1, 1)
    Wy = compute_WU(Z, y_col, S, R, T).flatten()  # W @ y
    yTWy = y @ Wy
    
    # Solve moment equations
    # E[y.T @ W @ y] = s2gxg * tr(W^2) + s2e * tr(W)
    # E[y.T @ y] = s2gxg * tr(W) + s2e * n
    A = np.array([
        [trW2_app, trW_app],
        [trW_app, n]
    ])
    
    b = np.array([yTWy, yTy])
    
    # Solve A @ [s2gxg, s2e].T = b
    x = np.linalg.solve(A, b)
    
    return x[0], x[1], A 

###sparse version

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
    idx = np.random.choice(m, m_sub, replace=False)
    Z_sub = Z[:, idx]
    
    W = build_W_batched(Z_sub)
    
    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)
    Le = np.sqrt(s2e) * np.eye(n)
    
    return Lgxg, Le



######pair sparse mode

def build_W_batched_sparse(Z, batch_size=5000, density=0.1):
    """
    Build W by sampling pairs (i,j) rather than sampling SNPs.
    
    density: fraction of SNPs (pairs scale as density^2)
    """
    n, m = Z.shape
    p_total = m * (m - 1) // 2
    
    idx_i, idx_j = np.triu_indices(m, k=1)
    
    # Sample pairs (density^2 to match SNP sampling)
    m_causal = density * m 
    p_causal = m_causal * (m_causal - 1) // 2
    causal_pairs = np.random.choice(p_total, size=p_causal, replace=False)
    idx_i = idx_i[causal_pairs]
    idx_j = idx_j[causal_pairs]
    
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
    Le = np.sqrt(s2e) * np.eye(n)
    
    return Lgxg, Le
