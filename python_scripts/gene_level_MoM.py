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


def run_experiment_MC(X, n_mc=100, s2gxg=0.9, s2e=0.1):
    """
    Run Monte Carlo simulations to estimate GxG and residual variance components.
    
    parameters:
    X: (n, m) raw genotype data
    n_mc: number of Monte Carlo replicates
    s2gxg: true GxG variance component
    s2e: true residual variance component   

    returns:
    DataFrame with estimated GxG and residual variance components across replicates.
    """
    gxgs = []
    es = []
    
    for _ in range(n_mc):
        Z, y = simulate_from_raw_only_W(X, s2gxg=s2gxg, s2e=s2e)
        gxg, e, _ = MoM_only_M(Z, y, nmc=40)
        gxgs.append(gxg)
        es.append(e)
    
    gxgs = np.array(gxgs).flatten()
    es = np.array(es).flatten()
    
    return pd.DataFrame({"gxgs": gxgs, "es": es})




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
    
    # Var(Zi*Zj) matrix
    numerator = sigma_matrix * (1 - 2*pi) * (1 - 2*pj)
    denominator = 4 * np.sqrt(pi * pj * (1 - pi) * (1 - pj))
    
    var_mtx = np.where(denominator > 1e-10, 
                       1 + numerator / denominator, 
                       1)

    # Mean of upper triangle
    rows, cols = np.triu_indices(m, k=1)
    correction_factor = np.mean(var_mtx[rows, cols])
    
    return correction_factor

def MoM_only_M_LD_Correction(Z, y, nmc=40, X=None, ld_correction=False):
    """
    Method of Moments estimator for GxG variance component using only the epistatic kernel.

    Parameters:
    -----------
    Z: (n, m) standardized genotype matrix
    y: (n,) phenotype vector
    nmc: number of Monte Carlo samples for trace estimation
    X: (n, m) original genotype matrix (0, 1, 2) for LD correction (required if ld_correction=True)
    ld_correction: whether to apply LD bias correction
    
    Returns:
    --------
    s2gxg: estimated GxG variance component
    s2e: estimated residual variance component
    A: 2x2 coefficient matrix from the moment equations
    correction_factor: LD correction factor (if ld_correction=True)
    """
    n = Z.shape[0]
    m = Z.shape[1]
    p = m * (m - 1) / 2 
    
    #  Correction
    if ld_correction:
        if X is None:
            raise ValueError("X (original genotype matrix) is required for LD correction")
        correction_factor = calculate_correction_factor(X)
    else:
        correction_factor = 1.0
    
    # Random vectors for stochastic trace estimation
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

    #  Raw trace estimates 
    trW_raw = 0.5 * (np.mean(np.sum(u * KKu, axis=0)) - np.mean(np.sum(u * Du, axis=0))) / p
    
    trKK2 = np.mean(np.sum(KKu * KKu, axis=0))
    trDK = np.mean(np.sum(Du * KKu, axis=0))
    trDD = np.mean(np.sum(Du * Du, axis=0))
    trW2_raw = (trKK2 - 2 * trDK + trDD) / (4 * p**2)
    
    # Quadratic forms
    yTy = y @ y
    yTWy_raw = 0.5 * (y.T @ (np.sum((Z @ (Z.T @ (y.reshape(-1, 1) * Z))) * Z, axis=1)) - (y @ (ZZ @ (ZZ.T @ y)))) / p

    #  Apply correction
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
    
    if ld_correction:
        return x1[0], x1[1], A, correction_factor
    else:
        return x1[0], x1[1], A