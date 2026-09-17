"""Verify the O(nm) standardized W operator.

Core claim:  if the weight matrix factorizes as  B = sum_a b_a c_a^T  (rank q),
then the quartic term
      A_{ts} = sum_{ij} B_ij Z_ti Z_tj Z_si Z_sj
equals    A = sum_a K_{b_a} .* K_{c_a},      K_b := Z diag(b) Z^T,
and every Hadamard product is applied by  (K_b .* K_c) u = diag(K_b D_u K_c).
"""
import numpy as np

rng = np.random.default_rng(0)


# ---------------------------------------------------------------- genotypes
def sim_geno(n, m, rho=0.85, maf_lo=0.05, maf_hi=0.5, seed=1):
    """AR(1)-LD haplotype-ish genotypes with varying MAF, column-standardized."""
    r = np.random.default_rng(seed)
    L = r.normal(size=(n, m))
    for j in range(1, m):
        L[:, j] = rho * L[:, j - 1] + np.sqrt(1 - rho**2) * L[:, j]
    f = r.uniform(maf_lo, maf_hi, m)
    thr = np.array([np.quantile(L[:, j], [(1 - f[j])**2, 1 - f[j]**2]) for j in range(m)])
    Gc = (L > thr[:, 0]).astype(float) + (L > thr[:, 1]).astype(float)
    Z = (Gc - Gc.mean(0)) / np.where(Gc.std(0) < 1e-12, 1.0, Gc.std(0))
    return Z


def dense_W(Z):
    """Exact standardized pairwise-interaction GRM (one gene)."""
    n, m = Z.shape
    i, j = np.triu_indices(m, k=1)
    H = Z[:, i] * Z[:, j]
    mu, sd = H.mean(0), H.std(0)
    ok = sd > 1e-10
    H[:, ok] = (H[:, ok] - mu[ok]) / sd[ok]
    H[:, ~ok] = 0.0
    return H @ H.T / len(i), len(i)


def weights(Z):
    n = Z.shape[0]
    G = Z.T @ Z / n
    D = Z * Z
    S = D.T @ D / n
    sig2 = S - G * G
    Bf = 1.0 / sig2                      # full (diagonal included)
    return G, D, sig2, Bf


# ------------------------------------------------- the primitive: K_b .* K_c
def hadamard_apply_exact(Z, b, c, u):
    """(K_b .* K_c) u  computed exactly = diag(K_b D_u K_c).  O(n^2 m), test only."""
    Kb = (Z * b) @ Z.T
    Kc = (Z * c) @ Z.T
    return np.einsum('ts,ts->t', Kb, Kc * u[None, :])


def hadamard_apply_mc(Z, Bfac, Cfac, U, Nmc, rng):
    """sum_a (K_{b_a} .* K_{c_a}) U  by the Rademacher diag estimator.  O(n m q Nmc).

    Bfac, Cfac : (m, q).   U : (n, c).
    """
    n, m = Z.shape
    q = Bfac.shape[1]
    out = np.zeros_like(U)
    for _ in range(Nmc):
        v = rng.choice([-1.0, 1.0], size=n)
        zeta = Z.T @ v                                  # (m,)      O(nm)
        Y = Z @ (Cfac * zeta[:, None])                  # (n, q)    O(nmq)
        for k in range(U.shape[1]):
            Wk = U[:, k][:, None] * Y                   # (n, q)
            Mk = Z.T @ Wk                               # (m, q)    O(nmq)
            out[:, k] += (Z @ (Bfac * Mk)) @ np.ones(q) * 0 + ((Z @ (Bfac * Mk)) * v[:, None]).sum(1)
    return out / Nmc


def hadamard_apply_mc_fast(Z, Bfac, Cfac, U, Nmc, rng):
    """Same, batched over the q factors properly (one gemm per stage)."""
    n, m = Z.shape
    q = Bfac.shape[1]
    nc = U.shape[1]
    out = np.zeros_like(U)
    for _ in range(Nmc):
        v = rng.choice([-1.0, 1.0], size=n)
        zeta = Z.T @ v                                   # (m,)
        Y = Z @ (Cfac * zeta[:, None])                   # (n,q)   K_{c_a} v
        for k in range(nc):
            M = Z.T @ (U[:, k][:, None] * Y)             # (m,q)
            P = Z @ (Bfac * M)                           # (n,q)   K_{b_a}(u.*K_{c_a}v)
            out[:, k] += (P * v[:, None]).sum(1)
    return out / Nmc


# =========================================================== C1: exact algebra
print("=" * 74)
print("C1  rank-q identity   A = sum_a K_{b_a} .* K_{c_a}   (exact eigen-factor)")
print("=" * 74)
n, m = 120, 40
Z = sim_geno(n, m, seed=3)
G, D, sig2, Bf = weights(Z)

# exact symmetric factorization of the FULL B (diagonal included)
lam, Q = np.linalg.eigh(Bf)
Fac = Q * np.sqrt(np.abs(lam))
sgn = np.sign(lam)
u = rng.normal(size=n)

# direct quartic sum
A_direct = np.einsum('ij,ti,tj,si,sj,s->t', Bf, Z, Z, Z, Z, u, optimize=True)
A_had = sum(sgn[a] * hadamard_apply_exact(Z, Fac[:, a], Fac[:, a], u) for a in range(m))
print(f"  ||A_had - A_direct|| / ||A_direct||  = {np.linalg.norm(A_had-A_direct)/np.linalg.norm(A_direct):.3e}")


# ================================================ C2: full W u from primitive
print()
print("=" * 74)
print("C2  full standardized  W u  assembled from the SAME primitive")
print("=" * 74)


def W_apply_factored(Z, Bfac, Cfac, sgn, U, exact=True, Nmc=0, rng=None):
    """W U for one gene from a rank-q factorization  B_full = sum_a s_a b_a c_a^T.

    Everything -- the quartic term, v_R and s_T -- is one primitive.
    """
    n, m = Z.shape
    p = m * (m - 1) // 2
    Dz = Z * Z
    beta = np.einsum('ia,ia,a->i', Bfac, Cfac, sgn)      # diag(B_full)
    q = Bfac.shape[1]

    def Phi(V):                                          # sum_a (K_b .* K_c) V, full B
        if exact:
            return np.column_stack([
                sum(sgn[a] * hadamard_apply_exact(Z, Bfac[:, a], Cfac[:, a], V[:, k])
                    for a in range(q)) for k in range(V.shape[1])])
        return hadamard_apply_mc_fast(Z, Bfac * sgn, Cfac, V, Nmc, rng)

    ones = np.ones((n, 1))
    vRf = Phi(ones)[:, 0] / n                            # v_R with diagonal
    v_R = vRf - Dz @ beta                                # zero the diagonal
    s_Tf = ones[:, 0] @ vRf / n
    s_T = s_Tf - beta.sum()

    U = np.atleast_2d(U.T).T
    s_u = U.sum(0)
    A = Phi(U) - Dz @ (beta[:, None] * (Dz.T @ U))
    return (A - np.outer(v_R, s_u) + np.outer(np.ones(n), s_T * s_u - U.T @ v_R)) / (2 * p)


Wd, p = dense_W(Z)
U = rng.normal(size=(n, 3))
Wf = W_apply_factored(Z, Fac, Fac, sgn, U, exact=True)
print(f"  ||W_fac U - W_dense U|| / ||W_dense U|| = "
      f"{np.linalg.norm(Wf-Wd@U)/np.linalg.norm(Wd@U):.3e}")

# sanity: unstandardized special case  B = 11^T  reproduces  K.*K - D D^T
one = np.ones((m, 1))
K = Z @ Z.T
A_ref = (K * K) @ U - (Z * Z) @ ((Z * Z).T @ U)
A_un = np.column_stack([hadamard_apply_exact(Z, one[:, 0], one[:, 0], U[:, k])
                        for k in range(U.shape[1])]) - (Z*Z) @ ((Z*Z).T @ U)
print(f"  unstandardized special case B=11^T -> K.*K - DD^T : rel.err "
      f"{np.linalg.norm(A_un-A_ref)/np.linalg.norm(A_ref):.3e}")

print()
print("=" * 74)
print("C3  MC Hadamard estimator: accuracy vs Nmc  (rank-q exact factorization)")
print("=" * 74)
u1 = rng.normal(size=(n, 1))
truth = Wd @ u1
for Nmc in [10, 50, 200, 1000]:
    est = W_apply_factored(Z, Fac * sgn, Fac, np.ones(m), u1, exact=False,
                           Nmc=Nmc, rng=np.random.default_rng(7))
    print(f"   Nmc={Nmc:5d}   rel.err = {np.linalg.norm(est-truth)/np.linalg.norm(truth):.4f}")
