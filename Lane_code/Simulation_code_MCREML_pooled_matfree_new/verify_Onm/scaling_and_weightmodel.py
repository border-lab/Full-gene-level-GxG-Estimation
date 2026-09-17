"""Two things:
 (A) the error law of the O(nm) Hadamard sketch:  relFro ~ sqrt(n/N),
     independent of m  -> N must scale with n, so fixed-accuracy cost is O(n^2 m).
 (B) a weight model that DOES factorize: sigma_ij^2 ~ alpha_i alpha_j psi(r_ij).
"""
import numpy as np


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


# ============================================================ (A) error law
print("=" * 78)
print("(A)  relative Frobenius error of  K.*K  by the Hutchinson-diag sketch")
print("     prediction:  relFro ~ sqrt(n/N)   (NO m dependence)")
print("=" * 78)
print(f"{'n':>6} {'m':>6} {'N':>6} {'relFro':>10} {'relFro*sqrt(N/n)':>18}")
for n in [200, 400, 800]:
    for m in [40, 80, 160]:
        Z, _ = sim_geno(n, m, seed=n + m)
        K = Z @ Z.T
        target = K * K
        N = 256
        r = np.random.default_rng(3)
        Ah = np.zeros((n, n))
        for _ in range(N):
            v = r.choice([-1.0, 1.0], size=n)
            Ah += K * (v[:, None] * (K @ v)[None, :])
        Ah /= N
        Ah = (Ah + Ah.T) / 2
        rel = np.linalg.norm(Ah - target) / np.linalg.norm(target)
        print(f"{n:>6} {m:>6} {N:>6} {rel:>10.4f} {rel*np.sqrt(N/n):>18.3f}")

print("\n  -> the last column is flat in BOTH n and m: the constant is O(1),")
print("     so N ~ n is required for O(1) accuracy, i.e. O(n^2 m) at fixed error.")

# ====================================================== (B) weight model
print()
print("=" * 78)
print("(B)  which model for sigma_ij^2 = Var(Z_i Z_j) actually fits?")
print("=" * 78)
for (mlo, mhi), tag in [((0.05, 0.5), "MAF 0.05-0.50"), ((0.01, 0.5), "MAF 0.01-0.50"),
                        ((0.2, 0.5), "MAF 0.20-0.50")]:
    n, m = 3000, 250
    Z, f = sim_geno(n, m, seed=11, maf=(mlo, mhi))
    G = Z.T @ Z / n
    Dz = Z * Z
    sig2 = Dz.T @ Dz / n - G * G
    iu, ju = np.triu_indices(m, k=1)
    r_ij, y = G[iu, ju], sig2[iu, ju]
    ss = ((y - y.mean())**2).sum()

    def r2(fit):
        return 1 - ((y - fit)**2).sum() / ss

    # model 0: constant 1
    m0 = r2(np.ones_like(y))
    # model 1: 1 + r^2       (Gaussian / HWE)
    m1 = r2(1 + r_ij**2)
    # model 2: separable only,  sigma^2 = a_i a_j   (a_i = sqrt(E[Z_i^4]) style)
    #          fit by least squares in log space
    Lg = np.log(np.maximum(y, 1e-6))
    Amat = np.zeros((len(y), m))
    Amat[np.arange(len(y)), iu] = 1.0
    Amat[np.arange(len(y)), ju] += 1.0
    la, *_ = np.linalg.lstsq(Amat, Lg, rcond=None)
    m2 = r2(np.exp(Amat @ la))
    # model 3: separable x psi(r):  log sigma^2 = la_i + la_j + poly(r)
    Xd = np.column_stack([Amat] + [r_ij**k for k in (2, 3, 4)])
    cd, *_ = np.linalg.lstsq(Xd, Lg, rcond=None)
    m3 = r2(np.exp(Xd @ cd))
    # kurtosis-based closed form:  a_i = sqrt(E[Z_i^4])
    k4 = (Z**4).mean(0)
    m4 = r2(np.sqrt(k4[iu] * k4[ju]))
    print(f"\n  {tag}:  sigma^2 in [{y.min():.2f}, {y.max():.2f}]")
    print(f"    sigma^2 = 1                        R^2 = {m0:7.4f}")
    print(f"    sigma^2 = 1 + r^2   (Gaussian)     R^2 = {m1:7.4f}")
    print(f"    sigma^2 = a_i a_j   (separable)    R^2 = {m2:7.4f}   <- rank-1 in B")
    print(f"    sigma^2 = sqrt(E Z_i^4 E Z_j^4)    R^2 = {m4:7.4f}   <- closed form, O(nm)")
    print(f"    sigma^2 = a_i a_j psi(r)           R^2 = {m3:7.4f}")

    # what the separable model buys: EXACT reduction to a rescaled Z
    alpha = np.sqrt(k4)
    Bsep = 1.0 / np.outer(alpha, alpha)
    Btrue = 1.0 / sig2
    print(f"    ==> ||B_sep - B_true||_F/||B_true||_F = "
          f"{np.linalg.norm(Bsep-Btrue)/np.linalg.norm(Btrue):.4f}")
