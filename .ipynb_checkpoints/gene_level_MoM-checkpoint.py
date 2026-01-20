# %%
import pandas as pd
import numpy as np
from scipy import sparse
import random
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
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


def MoM_only_M(Z, y, nmc=40):
    """
    Method of Moments estimator for GxG variance component using only the W kernel.
    
    Estimates s2gxg and s2e by solving a linear system derived from moment conditions.
    Uses randomized trace estimation for computational efficiency.
    
    Parameters:
        Z: (n, m) standardized genotype matrix
        y: (n,) phenotype vector
        nmc: number of Monte Carlo samples for trace estimation (default 40)
    
    Returns:
        s2gxg: estimated GxG variance component
        s2e: estimated residual variance component
        A: 2x2 coefficient matrix from the moment equations
    """
    n = Z.shape[0]
    m = Z.shape[1]
    p = m * (m - 1) / 2  # Number of unique SNP pairs
    
    # Random vectors for stochastic trace estimation
    u = np.random.randn(n, nmc)
    
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


#===============  weight function   ===============#

def weight_vector(X):
    """
    Compute optimal weights for standardized genotypes to minimize variance
    of pairwise products.

    The goal is to find weights w_i such that Var(w_i * Z_i * w_j * Z_j) ≈ 1
    for all pairs (i, j), where Z_i is the standardized genotype at SNP i.

    Parameters:
        X: (n, m) genotype matrix with n samples and m SNPs (values in {0, 1, 2})

    Returns:
        w: (m,) weight vector where w_i = exp(u_hat_i)
    """
    n, m = X.shape

    # Estimate allele frequencies p_i from E[X_i] = 2*p_i (Hardy-Weinberg)
    p_i = np.clip(X.mean(axis=0) / 2, 0.01, 0.99)

    # Estimate covariance matrix of genotypes
    sigma_ij = np.cov(X, rowvar=False)

    def var_ZiZj(pi, pj, sij):
        """Compute Var(Z_i * Z_j) given allele frequencies and covariance."""
        numerator = sij * (1 - 2*pi) * (1 - 2*pj)
        denominator = 4 * pi * pj * (1 - pi) * (1 - pj)
        return 1 + numerator / denominator

    # Set up least squares problem to find log-weights u_i such that
    # u_i + u_j ≈ -0.5 * log(Var(Z_i Z_j)) for all pairs (i < j)
    # This ensures w_i * w_j * Var(Z_i Z_j) ≈ 1
    p = m * (m - 1) // 2  # number of unique pairs (i < j)
    A = np.zeros((p, m))
    b = np.zeros(p)

    row = 0
    for i in range(m):
        for j in range(i + 1, m):
            # Constraint: u_i + u_j = c_ij
            A[row, i] = 1
            A[row, j] = 1
            # Target: c_ij = -0.5 * log(Var(Z_i Z_j))
            v = var_ZiZj(p_i[i], p_i[j], sigma_ij[i, j])
            v = max(v, 1e-10)  # numerical stability: ensure positive before log
            b[row] = -0.5 * np.log(v)
            row += 1

    # Solve overdetermined system via least squares: min ||Au - b||^2
    u_hat, *_ = np.linalg.lstsq(A, b, rcond=None)

    # Transform log-weights back to weights: w_i = exp(u_i)
    w = np.exp(u_hat)

    return w


def check_weighted_products(Z, w, min_dist=1):
    """
    Compute variance matrix of weighted pairwise products and extract statistics
    for distant pairs.

    Parameters:
        Z: (n, p) genotype matrix with n samples and p SNPs
        w: (n,) weight vector for each sample
        min_dist: minimum distance between SNP indices to be considered "distant"

    Returns:
        var_matrix: (p, p) symmetric matrix where entry (i,j) is the variance of
                    the elementwise product (Z[:,i]*w) * (Z[:,j]*w) across samples
        distant_vars: 1D array of variances for SNP pairs with index distance >= min_dist
        stats: dict with mean, std, and count of distant pair variances
    """
    n, p = Z.shape
    Zw = Z * w  # Weight each column of Z by the sample weights
    var_matrix = np.zeros((p, p))

    # Compute pairwise product variances (upper triangular, then mirror)
    for i in range(p):
        for j in range(i, p):
            prod = Zw[:, i] * Zw[:, j]
            var_matrix[i, j] = prod.var()
            var_matrix[j, i] = var_matrix[i, j]

    # Extract pairs with distance >= min_dist (off-diagonal elements)
    distant_vars = []
    for i in range(p):
        for j in range(i + min_dist, p):
            distant_vars.append(var_matrix[i, j])

    distant_vars = np.array(distant_vars)

    stats = {
        'offdiag_mean': distant_vars.mean(),
        'offdiag_std': distant_vars.std(),
        'n_pairs': len(distant_vars),
    }

    return var_matrix, distant_vars, stats

# %%



