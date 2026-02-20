import numpy as np

def calculate_varXiXj(X):
    """
    Calculate the pairwise product variance for each pair of SNPs in genotype matrix X.
    
    parameters:
    X: 2D numpy array of shape (n, m) with genotype data coded as 0, 1, 2.

    returns:
    A 1D numpy array containing the pairwise product variances for each SNP pair (i < j).
    """
    n, m = X.shape

    # Calculate AF and Covariance Matrix
    p = np.sum(X, axis=0) / (2 * n)
    sigma_matrix = np.cov(X, rowvar=False)
    
    # Broadcasting for all pairs i, j
    pi = p[:, np.newaxis]
    pj = p[np.newaxis, :]
    
    # The Formula applied to the whole matrix
    var_mtx = (
        -12 * pi**2 * pj**2
        + 4 * pi**2 * pj
        + 4 * pi * pj**2
        + 4 * pi * pj
        + 2 * pi * sigma_matrix
        + 2 * pj * sigma_matrix
        - 4 * pi * pj * sigma_matrix
        + sigma_matrix
    )
    # Extract upper triangular part (i < j) for array output
    rows, cols = np.triu_indices(m, k=1)
    return var_mtx[rows, cols]


def calculate_varZiZj(X):
    """
    Calculate the pairwise product variance for each pair of SNPs in Standardlized matrix Z using X.

    parameters:
    X: 2D numpy array of shape (n, m) with genotype data coded as 0, 1, 2.

    returns:
    A 1D numpy array containing the pairwise product variances for each SNP pair (i < j).
    """
    n, m = X.shape

    # Calculate AF and Covariance Matrix
    p = np.sum(X, axis=0) / (2 * n)
    sigma_matrix = np.cov(X, rowvar=False)
    
    # Broadcasting for all pairs i, j
    pi = p[:, np.newaxis]
    pj = p[np.newaxis, :]
    
    # The Formula applied to the whole matrix
    numerator = sigma_matrix * (1 - 2*pi) * (1 - 2*pj)
    denominator = 4 * pi * pj * (1 - pi) * (1 - pj)
    var_mtx = 1 + numerator / denominator

    # Extract upper triangular part (i < j) for array output
    rows, cols = np.triu_indices(m, k=1)
    return var_mtx[rows, cols]


# Calculate mean r² for both datasets
def get_mean_r2(genotype_matrix):
    Z = (genotype_matrix - np.mean(genotype_matrix, axis=0)) / np.std(genotype_matrix, axis=0)
    corr_matrix = np.corrcoef(Z, rowvar=False)
    r2_matrix = corr_matrix ** 2
    upper_tri = r2_matrix[np.triu_indices(r2_matrix.shape[0], k=1)]
    return np.mean(upper_tri)
