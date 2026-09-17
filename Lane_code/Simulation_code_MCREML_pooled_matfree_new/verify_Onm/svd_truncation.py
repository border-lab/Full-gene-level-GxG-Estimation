"""Section 2.3 of Some_update: the DETERMINISTIC low-rank route.

K_w = Z Z^T = sum_r lam_r q_r q_r^T  (rank <= m)  gives EXACTLY

    K_w .* K_w = sum_{r,k} lam_r lam_k (q_r .* q_k)(q_r .* q_k)^T

so truncating K_w at rank s gives an s^2-term apply costing O(n s^2), against
O(n m^2) exact.  Unlike the Hutchinson sketch this is data-dependent (it rides
the LD spectrum, not a random subspace) and every weight lam_r lam_k >= 0, so
the truncated K.*K stays PSD.

Question: (a) how small can s/m be, (b) does W_s stay PSD once the exact
-D D^T and the rank-one centring corrections are subtracted, (c) what does
REML do on it.
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
    return (Gc - Gc.mean(0)) / np.where(sd < 1e-12, 1.0, sd), f


def W_centred(Z):
    """The B = 11^T operator of section 2.6, formed densely as ground truth."""
    n, m = Z.shape
    p = m * (m - 1) / 2
    K = Z @ Z.T
    D = Z * Z
    return (K * K - D @ D.T) / (2 * p)


def W_trunc(Z, s):
    """Same, with K.*K replaced by its rank-s-in-K truncation. -D D^T kept exact."""
    n, m = Z.shape
    p = m * (m - 1) / 2
    K = Z @ Z.T
    lam, Q = np.linalg.eigh(K)
    idx = np.argsort(lam)[::-1][:s]
    lam, Q = lam[idx], Q[:, idx]
    Ks = (Q * lam) @ Q.T
    D = Z * Z
    return (Ks * Ks - D @ D.T) / (2 * p)


def reml_h2(y, W):
    """Profile REML for y ~ N(0, s2g W + s2e I), no fixed effects but the mean."""
    n = len(y)
    lam, U = np.linalg.eigh(W)
    Uty = U.T @ y

    def nll(logit):
        h2 = 1 / (1 + np.exp(-logit))
        d = h2 * lam + (1 - h2)
        d = np.maximum(d, 1e-12)
        s2 = (Uty**2 / d).sum() / n
        return 0.5 * (np.log(d).sum() + n * np.log(s2))

    r = minimize_scalar(nll, bounds=(-12, 12), method="bounded")
    return 1 / (1 + np.exp(-r.x))


print("=" * 86)
print("(A)  LD spectrum of K_w = Z Z^T  --  how fast does it decay?")
print("=" * 86)
print(f"{'n':>6} {'m':>6} {'rho':>5} {'s for 90% tr':>13} {'95%':>6} {'99%':>6} "
      f"{'eff.rank K':>11}")
for rho in (0.99, 0.9, 0.5):
    for n, m in [(400, 100), (400, 200), (800, 200)]:
        Z, _ = sim_geno(n, m, rho=rho, seed=n + m)
        lam = np.sort(np.linalg.eigvalsh(Z @ Z.T))[::-1]
        lam = np.maximum(lam, 0)
        c = np.cumsum(lam) / lam.sum()
        s90, s95, s99 = (int(np.searchsorted(c, t) + 1) for t in (.9, .95, .99))
        eff = lam.sum()**2 / (lam**2).sum()
        print(f"{n:>6} {m:>6} {rho:>5} {s90:>13} {s95:>6} {s99:>6} {eff:>11.1f}")

print()
print("=" * 86)
print("(B)  accuracy and definiteness of the rank-s truncated operator")
print("     cost ratio = s^2/m^2  (apply O(n s^2) vs exact O(n m^2))")
print("=" * 86)
n, m, rho = 400, 200, 0.9
Z, _ = sim_geno(n, m, rho=rho, seed=7)
Wex = W_centred(Z)
nrm = np.linalg.norm(Wex)
lmin_ex = np.linalg.eigvalsh(Wex).min()
print(f"  n={n} m={m} rho={rho}   exact: lam_min(W) = {lmin_ex:.3e}")
print(f"\n{'s':>5} {'s/m':>6} {'cost s^2/m^2':>13} {'relFro(W_s)':>12} "
      f"{'lam_min(W_s)':>13} {'PSD?':>5}")
for s in [5, 10, 20, 40, 80, 120, 160, 200]:
    Ws = W_trunc(Z, s)
    rel = np.linalg.norm(Ws - Wex) / nrm
    lm = np.linalg.eigvalsh(Ws).min()
    print(f"{s:>5} {s/m:>6.2f} {s*s/(m*m):>13.3f} {rel:>12.4f} {lm:>13.3e} "
          f"{'yes' if lm > -1e-8 * np.abs(Wex).max() else 'NO':>5}")

print()
print("=" * 86)
print("(C)  what REML does on the truncated operator   (true h2 = 0.5)")
print("=" * 86)
n, m = 400, 200
Z, _ = sim_geno(n, m, rho=0.9, seed=7)
Wex = W_centred(Z)
h2_true, nrep = 0.5, 20
lam_e, U_e = np.linalg.eigh(Wex)
half = U_e @ np.diag(np.sqrt(np.maximum(lam_e, 0))) @ U_e.T
rng = np.random.default_rng(0)
Y = [np.sqrt(h2_true) * (half @ rng.normal(size=n))
     + np.sqrt(1 - h2_true) * rng.normal(size=n) for _ in range(nrep)]

print(f"{'s':>7} {'relFro':>9} {'lam_min':>11} {'mean h2':>9} {'sd':>7}")
for s in [10, 20, 40, 80, 160, 200]:
    Ws = W_trunc(Z, s)
    rel = np.linalg.norm(Ws - Wex) / nrm
    lm = np.linalg.eigvalsh(Ws).min()
    est = np.array([reml_h2(y, Ws) for y in Y])
    print(f"{s:>7} {rel:>9.4f} {lm:>11.2e} {est.mean():>9.3f} {est.std():>7.3f}")
est = np.array([reml_h2(y, Wex) for y in Y])
print(f"{'exact':>7} {0.0:>9.4f} {lmin_ex:>11.2e} {est.mean():>9.3f} {est.std():>7.3f}")

print()
print("=" * 86)
print("(D)  head-to-head at matched cost:  sketch (Nmc probes, O(n m Nmc))")
print("     vs truncation (rank s, O(n s^2)).  Matched when Nmc*m = s^2.")
print("=" * 86)


def W_sketch(Z, N, seed=3):
    n, m = Z.shape
    p = m * (m - 1) / 2
    K = Z @ Z.T
    r = np.random.default_rng(seed)
    Ah = np.zeros((n, n))
    for _ in range(N):
        v = r.choice([-1.0, 1.0], size=n)
        Ah += K * (v[:, None] * (K @ v)[None, :])
    Ah /= N
    Ah = (Ah + Ah.T) / 2
    D = Z * Z
    return (Ah - D @ D.T) / (2 * p)


print(f"{'budget (flops/nm)':>18} {'method':>12} {'param':>8} {'relFro':>9} "
      f"{'lam_min':>11} {'h2':>7}")
for budget in [40, 80, 160, 320]:
    Nmc = budget
    s = int(np.sqrt(budget * m))
    Wsk = W_sketch(Z, Nmc)
    Wtr = W_trunc(Z, min(s, m))
    for tag, Wc, prm in [("sketch", Wsk, f"N={Nmc}"), ("trunc", Wtr, f"s={min(s,m)}")]:
        rel = np.linalg.norm(Wc - Wex) / nrm
        lm = np.linalg.eigvalsh(Wc).min()
        h2 = np.mean([reml_h2(y, Wc) for y in Y])
        print(f"{budget:>18} {tag:>12} {prm:>8} {rel:>9.4f} {lm:>11.2e} {h2:>7.3f}")
