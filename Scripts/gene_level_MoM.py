import numpy as np

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






