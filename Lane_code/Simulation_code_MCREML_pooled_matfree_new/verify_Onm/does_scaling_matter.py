"""Does the 1/sigma_ij^2 pair-scaling matter, or is centring enough?

Three kernels, all on the same genotypes:
  W_std   : full pair standardization   B=1/sig2, R=mu/sig2, T=mu^2/sig2
  W_cen   : centring only               B=11', R=mu,       T=mu^2
  W_raw   : nothing                     B=11', R=0,        T=0
Simulate y from W_std, fit exact REML with each.
"""
import numpy as np
from scipy.optimize import minimize_scalar


def sim_geno(n, m, rho=0.9, seed=1, maf=(0.05, 0.5)):
    r = np.random.default_rng(seed)
    L = r.normal(size=(n, m))
    for j in range(1, m):
        L[:, j] = rho * L[:, j - 1] + np.sqrt(1 - rho**2) * L[:, j]
    f = r.uniform(*maf, m)
    thr = np.array([np.quantile(L[:, j], [(1 - f[j])**2, 1 - f[j]**2]) for j in range(m)])
    Gc = (L > thr[:, 0]).astype(float) + (L > thr[:, 1]).astype(float)
    sd = Gc.std(0)
    return (Gc - Gc.mean(0)) / np.where(sd < 1e-12, 1.0, sd)


def kernels(Z):
    n, m = Z.shape
    p = m * (m - 1) // 2
    K = Z @ Z.T
    Dz = Z * Z
    G = Z.T @ Z / n
    sig2 = Dz.T @ Dz / n - G * G
    one = np.ones(n)

    def build(B, R, T):
        B, R, T = B.copy(), R.copy(), T.copy()
        np.fill_diagonal(B, 0.0); np.fill_diagonal(R, 0.0); np.fill_diagonal(T, 0.0)
        A = np.einsum('ij,ti,tj,si,sj->ts', B, Z, Z, Z, Z, optimize=True)
        vR = np.sum(Z * (Z @ R), axis=1)
        sT = T.sum()
        return (A - np.outer(vR, one) - np.outer(one, vR) + sT * np.outer(one, one)) / (2 * p)

    Bi = 1.0 / sig2
    W_std = build(Bi, G * Bi, G * G * Bi)
    W_cen = build(np.ones((m, m)), G, G * G)
    W_raw = build(np.ones((m, m)), np.zeros((m, m)), np.zeros((m, m)))
    return W_std, W_cen, W_raw


def reml_h2(y, W):
    lam, Q = np.linalg.eigh((W + W.T) / 2)
    Qy = Q.T @ y
    n = len(y)

    def negll(h2):
        d = np.maximum(h2 * lam + (1 - h2), 1e-12)
        s2 = (Qy**2 / d).sum() / n
        return 0.5 * (np.log(d).sum() + n * np.log(s2))

    return minimize_scalar(negll, bounds=(1e-6, 1 - 1e-6), method='bounded',
                           options={'xatol': 1e-7}).x


n, m = 500, 80
print(f"n={n} m={m}")
W_std, W_cen, W_raw = kernels(sim_geno(n, m, seed=5))
Z = sim_geno(n, m, seed=5)

for nm, Wk in [("W_std", W_std), ("W_cen", W_cen), ("W_raw", W_raw)]:
    ev = np.linalg.eigvalsh((Wk + Wk.T) / 2)
    print(f"  {nm}: tr/n={np.trace(Wk)/n:7.4f}  minEV={ev[0]:+.4f}  maxEV={ev[-1]:8.3f}  "
          f"eff.rank={np.trace(Wk)**2/np.trace(Wk@Wk):6.1f}")

print(f"\n  ||W_cen - W_std||_F / ||W_std||_F = "
      f"{np.linalg.norm(W_cen-W_std)/np.linalg.norm(W_std):.4f}")
print(f"  ||W_raw - W_std||_F / ||W_std||_F = "
      f"{np.linalg.norm(W_raw-W_std)/np.linalg.norm(W_std):.4f}")
# correlation of off-diagonal entries -- what REML actually keys on
iu = np.triu_indices(n, 1)
print(f"  corr(offdiag W_cen, W_std) = {np.corrcoef(W_cen[iu], W_std[iu])[0,1]:.4f}")
print(f"  corr(offdiag W_raw, W_std) = {np.corrcoef(W_raw[iu], W_std[iu])[0,1]:.4f}")

print("\n  h2 recovered (y simulated from W_std), 20 replicates:")
lam, Q = np.linalg.eigh((W_std + W_std.T) / 2)
L = Q @ np.diag(np.sqrt(np.maximum(lam, 0)))
res = {k: [] for k in ("W_std", "W_cen", "W_raw")}
for rep in range(20):
    r = np.random.default_rng(1000 + rep)
    g = L @ r.normal(size=n); g *= np.sqrt(0.5) / g.std()
    e = r.normal(size=n);     e *= np.sqrt(0.5) / e.std()
    y = g + e; y -= y.mean()
    for k, Wk in [("W_std", W_std), ("W_cen", W_cen), ("W_raw", W_raw)]:
        res[k].append(reml_h2(y, Wk))
for k in res:
    a = np.array(res[k])
    print(f"    fit with {k}: mean h2 = {a.mean():.4f}  sd = {a.std():.4f}   (true 0.5)")
