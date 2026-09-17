"""Is the weight matrix a dot-product kernel, i.e. sigma_ij^2 ~ psi(r_ij)?

If yes, then B = 1/psi(r), R = r/psi(r), T = r^2/psi(r) are all univariate
functions of the LD correlation, hence polynomials in the Hadamard powers of
G = Z'Z/n, hence rank-N factorizable by a Rademacher tensor sketch in O(nmkN).
Then v_R and s_T -- the whole u-INDEPENDENT setup -- cost O(nm), with NO
n-dimensional Hutchinson error.
"""
import numpy as np

rng = np.random.default_rng(0)


def sim_geno(n, m, rho=0.9, seed=1, maf=(0.05, 0.5)):
    r = np.random.default_rng(seed)
    L = r.normal(size=(n, m))
    for j in range(1, m):
        L[:, j] = rho * L[:, j - 1] + np.sqrt(1 - rho**2) * L[:, j]
    f = r.uniform(*maf, m)
    thr = np.array([np.quantile(L[:, j], [(1 - f[j])**2, 1 - f[j]**2]) for j in range(m)])
    Gc = (L > thr[:, 0]).astype(float) + (L > thr[:, 1]).astype(float)
    sd = Gc.std(0)
    return (Gc - Gc.mean(0)) / np.where(sd < 1e-12, 1.0, sd), f


n, m = 2000, 300
Z, f = sim_geno(n, m, seed=5)
G = Z.T @ Z / n
Dz = Z * Z
S = Dz.T @ Dz / n
sig2 = S - G * G

iu, ju = np.triu_indices(m, k=1)
r_ij = G[iu, ju]
s2_ij = sig2[iu, ju]

print("=" * 78)
print("Q1  is sigma_ij^2 a function of r_ij alone?")
print("=" * 78)
print(f"  r  range [{r_ij.min():+.3f}, {r_ij.max():+.3f}]   "
      f"sigma^2 range [{s2_ij.min():.3f}, {s2_ij.max():.3f}]")
# best univariate predictor: bin on r, compare within-bin scatter to total
nb = 40
edges = np.quantile(r_ij, np.linspace(0, 1, nb + 1))
idx = np.clip(np.digitize(r_ij, edges[1:-1]), 0, nb - 1)
pred = np.array([s2_ij[idx == b].mean() for b in range(nb)])[idx]
r2 = 1 - ((s2_ij - pred)**2).sum() / ((s2_ij - s2_ij.mean())**2).sum()
print(f"  R^2 of the best univariate psi(r)      : {r2:.4f}")

# polynomial fit  sigma^2 = 1 + a2 r^2 + a4 r^4 ...  (must ->1 at r=0)
for deg in [2, 4, 6]:
    ks = list(range(2, deg + 1))
    Xd = np.column_stack([r_ij**k for k in ks])
    coef, *_ = np.linalg.lstsq(Xd, s2_ij - 1.0, rcond=None)
    fit = 1.0 + Xd @ coef
    rr = 1 - ((s2_ij - fit)**2).sum() / ((s2_ij - s2_ij.mean())**2).sum()
    print(f"  R^2 of  1 + sum_(k=2..{deg}) a_k r^k       : {rr:.4f}   coef={np.round(coef,4)}")
    if deg == 4:
        coef4, ks4 = coef, ks

print("\n  Gaussian/HWE prediction sigma^2 = 1 + r^2 : "
      f"R^2 = {1-((s2_ij-(1+r_ij**2))**2).sum()/((s2_ij-s2_ij.mean())**2).sum():.4f}")


# ==================================================== tensor sketch of G^{ok}
print()
print("=" * 78)
print("Q2  Rademacher tensor sketch of the Hadamard powers  G^{ok}")
print("=" * 78)


def sketch_hadamard_power(Z, k, N, rng):
    """Phi (m,N) with  (1/N) Phi Phi^T ~ G^{ok},  G = Z'Z/n.   O(n m k N)."""
    n, m = Z.shape
    Phi = np.ones((m, N))
    for _ in range(k):
        Wr = rng.choice([-1.0, 1.0], size=(n, N))
        Phi *= (Z.T @ Wr) / np.sqrt(n)          # (m,N)
    return Phi


for k in [1, 2, 3]:
    Gk = G ** k
    print(f"  k={k}: ||G^ok||_F = {np.linalg.norm(Gk):8.2f}")
    for N in [64, 256, 1024]:
        Phi = sketch_hadamard_power(Z, k, N, np.random.default_rng(k * 100 + N))
        approx = Phi @ Phi.T / N
        print(f"       N={N:>5}  relFro = "
              f"{np.linalg.norm(approx-Gk)/np.linalg.norm(Gk):.4f}")


# ============================================ Q3  v_R and s_T in O(nm)
print()
print("=" * 78)
print("Q3  v_R and s_T from the sketch -- NO n-dimensional Hutchinson")
print("=" * 78)

# exact reference with the TRUE empirical weights
Bf = 1.0 / sig2
Rw, Tw = G * Bf, G * G * Bf
beta = np.diag(Bf).copy()
R0, T0 = Rw.copy(), Tw.copy()
np.fill_diagonal(R0, 0.0)
np.fill_diagonal(T0, 0.0)
vR_exact = np.sum(Z * (Z @ R0), axis=1)
sT_exact = T0.sum()
print(f"  exact  ||v_R|| = {np.linalg.norm(vR_exact):.4f}   s_T = {sT_exact:.4f}")

# model route: psi(r) = 1 + a2 r^2 + a4 r^4  ->  R(r) = r/psi,  T(r) = r^2/psi
rg = np.linspace(-1, 1, 4001)
psi = 1.0 + sum(c * rg**k for c, k in zip(coef4, ks4))
psi = np.maximum(psi, 1e-3)
for deg in [5, 9, 13]:
    cR = np.polyfit(rg, rg / psi, deg)[::-1]          # R(r) ~ sum_k cR_k r^k
    cT = np.polyfit(rg, rg**2 / psi, deg)[::-1]
    # v_R = sum_k cR_k (Z .* (Z G^{ok})) 1 ,  sketched G^{ok} = Phi_k Phi_k'/N
    N = 512
    vR = np.zeros(n)
    sT = 0.0
    for k, ck in enumerate(cR):
        if abs(ck) < 1e-14:
            continue
        if k == 0:
            vR += ck * (Z.sum(1) ** 2)                 # G^{o0} = 11'
            continue
        Phi = sketch_hadamard_power(Z, k, N, np.random.default_rng(1000 + k))
        vR += ck * ((Z @ Phi) ** 2).sum(1) / N
    for k, ck in enumerate(cT):
        if abs(ck) < 1e-14:
            continue
        if k == 0:
            sT += ck * m * m
            continue
        Phi = sketch_hadamard_power(Z, k, N, np.random.default_rng(1000 + k))
        sT += ck * (Phi.sum(0) ** 2).sum() / N
    # exact O(nm) diagonal removal:  R_ii = beta_i, T_ii = beta_i
    vR -= Dz @ beta
    sT -= beta.sum()
    print(f"  poly deg {deg:>2}, N={N}:  ||v_R|| rel.err = "
          f"{np.linalg.norm(vR-vR_exact)/np.linalg.norm(vR_exact):.4f}   "
          f"s_T rel.err = {abs(sT-sT_exact)/abs(sT_exact):.4f}")
