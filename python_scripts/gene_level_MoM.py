import numpy as np
import pandas as pd
from scipy.linalg import cholesky

#===============  Core method   ===============#
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
    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0)
    n, m = Z.shape
    
    # Compute linear kernel (genetic relationship matrix)
    KZ = Z @ Z.T
    
    # Element-wise squared genotypes for diagonal correction
    D = Z * Z
    
    # Number of unique SNP pairs
    p = m * (m - 1) / 2
    
    # GxG interaction kernel: W = 0.5 * (KZ ⊙ KZ - D @ D.T) / p
    # This captures pairwise epistatic interactions between SNPs
    W = 0.5 * (KZ * KZ - D @ D.T) / p
    
    # Phenotype covariance: Var(y) = s2gxg * W + s2e * I
    cov_y = s2gxg * W + s2e * np.eye(n)
    
    # Sample phenotype from multivariate normal
    y = np.random.multivariate_normal(mean=np.zeros(n), cov=cov_y)
    y -= np.mean(y)  # Center the phenotype
    
    return Z, y


###raw method
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
    
    # D = Z ⊙ Z (element-wise square)
    ZZ = Z * Z
    
    # Compute D @ D.T @ u for trace estimation
    Du = ZZ @ (ZZ.T @ u)

    # Compute (K ⊙ K) @ u where K = Z @ Z.T
    # Done without forming K explicitly to save memory
    KKu = np.zeros((n, nmc))
    for i in range(nmc):
        # Z_tilde = diag(u_i) @ Z
        Z_tilde = u[:, i].reshape(-1, 1) * Z
        # M = Z.T @ Z_tilde = Z.T @ diag(u_i) @ Z
        M = Z.T @ Z_tilde
        # ZM @ Z.T gives (K ⊙ K) @ u_i via the identity
        ZM = Z @ M
        KKu[:, i] = np.sum(ZM * Z, axis=1)

    # tr(W) ≈ u.T @ W @ u, where W = 0.5 * (K⊙K - D@D.T) / p
    trW_app = 0.5 * (np.mean(np.sum(u * KKu, axis=0)) - np.mean(np.sum(u * Du, axis=0))) / p
    
    # tr(W^2) estimation via tr((K⊙K)^2), tr((K⊙K)(D@D.T)), tr((D@D.T)^2)
    trKK2 = np.mean(np.sum(KKu * KKu, axis=0))
    trDK = np.mean(np.sum(Du * KKu, axis=0))
    trDD = np.mean(np.sum(Du * Du, axis=0))
    trW2_app = (trKK2 - 2 * trDK + trDD) / (4 * p**2)
    
    #Quadratic forms
    yTy = y @ y  # y.T @ y
    
    # y.T @ W @ y = 0.5 * (y.T @ (K⊙K) @ y - y.T @ D @ D.T @ y) / p
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
    """
    # Calculate correction factor
    factor = calculate_correction_factor(real_data)
    
    # Standardize
    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0)
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


def simulate_Cholesky_from_raw_sparse(real_data, s2gxg=0.9, s2e=0.1, causal_fraction=0.25,stability=1e-8, seed=None):
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
    
    # D = Z ⊙ Z (element-wise square)
    ZZ = Z * Z
    
    # Compute D @ D.T @ u
    Du = ZZ @ (ZZ.T @ u)

    # Compute (K ⊙ K) @ u where K = Z @ Z.T
    KKu = np.zeros((n, nmc))
    for i in range(nmc):
        Z_tilde = u[:, i].reshape(-1, 1) * Z
        M = Z.T @ Z_tilde
        ZM = Z @ M
        KKu[:, i] = np.sum(ZM * Z, axis=1)

    # W_raw = 0.5 * (K⊙K - D@D.T) / p  (before standardization)
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
    
    return x1[0], x1[1], A, scale


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
    p = np.clip(p, 0.01, 0.99)
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
    Z = (real_data - real_data.mean(axis=0)) / real_data.std(axis=0)
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


def MoM_only_M_LD_Correction(X, Z, y, nmc=40):
    """
    Method of Moments estimator for GxG variance component with LD correction.
    
    Parameters:
    -----------
    X: (n, m) original genotype matrix (0, 1, 2) for LD correction
    Z: (n, m) standardized genotype matrix
    y: (n,) phenotype vector
    nmc: number of Monte Carlo samples for trace estimation
    
    Returns:
    --------
    s2gxg: estimated GxG variance component
    s2e: estimated residual variance component
    A: 2x2 coefficient matrix from the moment equations
    """
    n = Z.shape[0]
    m = Z.shape[1]
    p = m * (m - 1) / 2 
    
    # ====== LD Correction Factor ======
    correction_factor = calculate_correction_factor(X)
 
    # ====== Random vectors for stochastic trace estimation ======
    u = np.random.choice([-1, 1], size=(n, nmc))
    
    # D = Z ⊙ Z (element-wise square)
    ZZ = Z * Z
    
    # Compute D @ D.T @ u for trace estimation
    Du = ZZ @ (ZZ.T @ u)
    
    # Compute (K ⊙ K) @ u where K = Z @ Z.T
    KKu = np.zeros((n, nmc))
    for i in range(nmc):
        Z_tilde = u[:, i].reshape(-1, 1) * Z
        M = Z.T @ Z_tilde
        ZM = Z @ M
        KKu[:, i] = np.sum(ZM * Z, axis=1)
    
    # ====== Raw trace estimates ======
    trW_raw = 0.5 * (np.mean(np.sum(u * KKu, axis=0)) - np.mean(np.sum(u * Du, axis=0))) / p
    
    trKK2 = np.mean(np.sum(KKu * KKu, axis=0))
    trDK = np.mean(np.sum(Du * KKu, axis=0))
    trDD = np.mean(np.sum(Du * Du, axis=0))
    trW2_raw = (trKK2 - 2 * trDK + trDD) / (4 * p**2)
    
    # Quadratic forms
    yTy = y @ y
    yTWy_raw = 0.5 * (y.T @ (np.sum((Z @ (Z.T @ (y.reshape(-1, 1) * Z))) * Z, axis=1)) - (y @ (ZZ @ (ZZ.T @ y)))) / p
    
    # ====== Apply LD correction ======
    trW_app = trW_raw / correction_factor
    trW2_app = trW2_raw / (correction_factor ** 2)
    yTWy = yTWy_raw / correction_factor
    
    # ====== Solve moment equations ======
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



###############
# New function
###############

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
    
    # D = Z ⊙ Z (element-wise square)
    ZZ = Z * Z
    
    # Compute D @ D.T @ u for trace estimation
    Du = ZZ @ (ZZ.T @ u)

    # Compute (K ⊙ K) @ u where K = Z @ Z.T
    # Done without forming K explicitly to save memory
    KKu = np.zeros((n, nmc))
    for i in range(nmc):
        # Z_tilde = diag(u_i) @ Z
        Z_tilde = u[:, i].reshape(-1, 1) * Z
        # M = Z.T @ Z_tilde = Z.T @ diag(u_i) @ Z
        M = Z.T @ Z_tilde
        # ZM @ Z.T gives (K ⊙ K) @ u_i via the identity
        ZM = Z @ M
        KKu[:, i] = np.sum(ZM * Z, axis=1)

    # tr(W) ≈ u.T @ W @ u, where W = 0.5 * (K⊙K - D@D.T) / p
    trW_app = 0.5 * (np.mean(np.sum(u * KKu, axis=0)) - np.mean(np.sum(u * Du, axis=0))) / p
    
    # tr(W^2) estimation via tr((K⊙K)^2), tr((K⊙K)(D@D.T)), tr((D@D.T)^2)
    trKK2 = np.mean(np.sum(KKu * KKu, axis=0))
    trDK = np.mean(np.sum(Du * KKu, axis=0))
    trDD = np.mean(np.sum(Du * Du, axis=0))
    trW2_app = (trKK2 - 2 * trDK + trDD) / (4 * p**2)
    
    #Quadratic forms
    yTy = y @ y  # y.T @ y
    
    # y.T @ W @ y = 0.5 * (y.T @ (K⊙K) @ y - y.T @ D @ D.T @ y) / p
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
    p_causal = int(p_total * density * density)
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
