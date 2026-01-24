import numpy as np

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