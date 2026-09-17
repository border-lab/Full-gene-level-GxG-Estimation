"""Rank-q factorization of B = 1/sigma^2 WITHOUT forming B, in O(n m s).

Landmark SNPs S (|S| = s).  The exact columns B[:,S] need only
    G[:,S] = Z' Z_S / n        (O(nms))
    S[:,S] = D' D_S / n        (O(nms)),  D = Z .* Z
so  B[:,S] = 1/(S[:,S] - G[:,S]^2)  costs O(ms) more.  Then
    B ~ C M^+ C',   C = B[:,S], M = B[S,S]
and an eigendecomposition of the s-by-s M gives  B ~ sum_a s_a b_a b_a'.
Storage O(ms), never O(m^2).
"""
import numpy as np


def sim_geno(n, m, rho=0.9, seed=1):
    r = np.random.default_rng(seed)
    L = r.normal(size=(n, m))
    for j in range(1, m):
        L[:, j] = rho * L[:, j - 1] + np.sqrt(1 - rho**2) * L[:, j]
    f = r.uniform(0.05, 0.5, m)
    thr = np.array([np.quantile(L[:, j], [(1 - f[j])**2, 1 - f[j]**2]) for j in range(m)])
    Gc = (L > thr[:, 0]).astype(float) + (L > thr[:, 1]).astype(float)
    sd = Gc.std(0)
    return (Gc - Gc.mean(0)) / np.where(sd < 1e-12, 1.0, sd)


def nystrom_B(Z, S, ridge=1e-8):
    """Rank-|S| factors (Fac, sgn) with  B ~ (Fac*sgn) Fac' ,  cost O(n m s)."""
    n, m = Z.shape
    D = Z * Z
    Gs = Z.T @ Z[:, S] / n                       # (m,s)   O(nms)
    Ss = D.T @ D[:, S] / n                       # (m,s)   O(nms)
    C = 1.0 / (Ss - Gs * Gs)                     # (m,s)   exact columns of B
    M = C[S, :]                                  # (s,s)   = B[S,S]
    M = (M + M.T) / 2
    lam, Q = np.linalg.eigh(M)
    keep = np.abs(lam) > ridge * np.abs(lam).max()
    lam, Q = lam[keep], Q[:, keep]
    # B ~ C M^+ C' = (C Q |lam|^{-1/2}) diag(sign lam) (...)'
    Fac = C @ (Q / np.sqrt(np.abs(lam)))
    return Fac, np.sign(lam)


n, m = 400, 200
Z = sim_geno(n, m, seed=5)
p = m * (m - 1) // 2
Dz = Z * Z
Gm = Z.T @ Z / n
Bf = 1.0 / (Dz.T @ Dz / n - Gm * Gm)

iu, ju = np.triu_indices(m, k=1)
H = Z[:, iu] * Z[:, ju]
mu, sd = H.mean(0), H.std(0)
ok = sd > 1e-10
H[:, ok] = (H[:, ok] - mu[ok]) / sd[ok]
H[:, ~ok] = 0.0
W = H @ H.T / p

rng = np.random.default_rng(0)
u = rng.normal(size=n)
truth = W @ u
one_n = np.ones(n)
s_u = u.sum()


def Wu_from_factors(Fac, sgn):
    """Exact O(n m^2) apply, but with the weights carried as rank-q factors."""
    beta = np.einsum('ia,ia,a->i', Fac, Fac, sgn)
    Bq = (Fac * sgn) @ Fac.T
    B0 = Bq.copy(); np.fill_diagonal(B0, 0.0)
    M = Z.T @ (u[:, None] * Z)
    A = np.sum(Z * (Z @ (B0 * M)), axis=1)
    R0 = Gm * B0; T0 = Gm * Gm * B0
    v_R = np.sum(Z * (Z @ R0), axis=1)
    return (A - s_u * v_R + (T0.sum() * s_u - u @ v_R) * one_n) / (2 * p)


print("=" * 74)
print(f"n={n}  m={m}   Nystrom rank-s factorization of B (never forms B)")
print("=" * 74)
print(f"{'s':>5} {'||B_s-B||_F/||B||_F':>22} {'||W_s u - W u||/||W u||':>26}")
for s in [2, 4, 8, 16, 32, 64]:
    # uniformly spaced landmarks (LD is local -> spread them along the gene)
    S = np.unique(np.linspace(0, m - 1, s).astype(int))
    Fac, sgn = nystrom_B(Z, S)
    Bq = (Fac * sgn) @ Fac.T
    est = Wu_from_factors(Fac, sgn)
    print(f"{len(S):>5} {np.linalg.norm(Bq-Bf)/np.linalg.norm(Bf):>22.4f}"
          f" {np.linalg.norm(est-truth)/np.linalg.norm(truth):>26.4f}")

print("\nreference: q=1 rank-one eigen-truncation of the FULL B")
lam, Q = np.linalg.eigh(Bf)
o = np.argsort(-np.abs(lam))[:1]
Fac = Q[:, o] * np.sqrt(np.abs(lam[o]))
print(f"   exact-eigen q=1 -> ||W_q u - W u||/||W u|| = "
      f"{np.linalg.norm(Wu_from_factors(Fac, np.sign(lam[o]))-truth)/np.linalg.norm(truth):.4f}")

print("\nand the plain centring-only kernel (B = 11', no scaling at all):")
Fac1 = np.ones((m, 1))
print(f"   B = 11'         -> ||W_q u - W u||/||W u|| = "
      f"{np.linalg.norm(Wu_from_factors(Fac1, np.ones(1))-truth)/np.linalg.norm(truth):.4f}")
