"""
Joint (n,m) scaling of the epistatic MoM floor.
Prediction (independence, Marchenko-Pastur): with gamma = n/p, p = C(m,2),
   r4 = (trW^2)^2/trW^4  ~  n * (1+g)^2 / (1 + 6g + 6g^2 + g^3),
so floor relSE = sqrt(2/r4) ~ sqrt(2/min(n,p)).
Consequence: fixed m -> plateau at sqrt(2/p)=~2/m; grow m with sqrt(n) -> 1/sqrt(n).
Under LD, p must be replaced by the *effective* number of independent interaction dims.
"""
import numpy as np
from numpy.linalg import inv, eigvalsh
from scipy.stats import norm

def Z_indep(n, m, rng):
    p = rng.uniform(0.05, 0.5, m)
    X = rng.binomial(2, p, (n, m)).astype(float)
    X = X[:, X.std(0) > 0]
    return (X - X.mean(0)) / X.std(0)

def Z_LD(n, m, rng, rho):
    idx = np.arange(m); C = rho ** np.abs(idx[:, None] - idx[None, :])
    L = np.linalg.cholesky(C + 1e-8 * np.eye(m))
    pr = rng.uniform(0.05, 0.5, m); thr = norm.ppf(pr)
    hap = lambda: (rng.standard_normal((n, m)) @ L.T < thr).astype(float)
    X = hap() + hap(); X = X[:, X.std(0) > 0]
    return (X - X.mean(0)) / X.std(0)

def W_epi(Z):
    n, m = Z.shape
    K = Z @ Z.T; D = Z * Z
    p = m * (m - 1) / 2
    return (K * K - D @ D.T) / (2 * p)

def r4_and_relSE(W, s2=0.9, s2e=0.1):
    n = W.shape[0]
    lam = np.clip(eigvalsh(W), 0, None)
    trW, trW2, trW4 = lam.sum(), (lam**2).sum(), (lam**4).sum()
    r4 = trW2**2 / trW4
    d = s2 * lam + s2e
    T = np.array([[trW2, trW], [trW, n]])
    Cq = 2 * np.array([[(d**2*lam**2).sum(), (d**2*lam).sum()],
                       [(d**2*lam).sum(),    (d**2).sum()]])
    relSE = np.sqrt((inv(T) @ Cq @ inv(T))[0, 0]) / s2
    return r4, relSE

def mp_r4(n, m):
    p = m * (m - 1) / 2
    g = n / p
    return n * (1 + g)**2 / (1 + 6*g + 6*g**2 + g**3)

rng = np.random.default_rng(7)

print("=== MP prediction check (independent SNPs) ===")
print(f"{'n':>5} {'m':>4} {'p':>7} {'gamma':>7} {'r4_emp':>8} {'r4_MP':>8} {'relSE':>7} {'sqrt(2/r4)':>10}")
for n, m in [(2000,30),(2000,100),(2000,300),(500,60),(4000,60),(1000,200)]:
    W = W_epi(Z_indep(n, m, rng))
    r4, se = r4_and_relSE(W)
    p = m*(m-1)/2
    print(f"{n:>5} {m:>4} {int(p):>7} {n/p:>7.3f} {r4:>8.1f} {mp_r4(n,m):>8.1f} {se:>7.4f} {np.sqrt(2/r4):>10.4f}")

print("\n=== Regime A: FIXED m=30 (plateau expected, floor ~ sqrt(2/435)=0.068) ===")
print(f"{'n':>6} {'m':>4} {'r4':>7} {'relSE':>7}")
for n in [500,1000,2000,4000]:
    r4, se = r4_and_relSE(W_epi(Z_indep(n, 30, rng)))
    print(f"{n:>6} {30:>4} {r4:>7.1f} {se:>7.4f}")

print("\n=== Regime B: m grows ~ sqrt(n)  (gamma fixed -> relSE ~ 1/sqrt(n)) ===")
print(f"{'n':>6} {'m':>4} {'gamma':>6} {'r4':>7} {'relSE':>7} {'relSE*sqrt(n)':>13}")
for n in [500,1000,2000,4000]:
    m = int(round(np.sqrt(n)))
    r4, se = r4_and_relSE(W_epi(Z_indep(n, m, rng)))
    p = m*(m-1)/2
    print(f"{n:>6} {m:>4} {n/p:>6.2f} {r4:>7.1f} {se:>7.4f} {se*np.sqrt(n):>13.3f}")

print("\n=== Same sqrt(n) growth but UNDER LD (rho=.95): does growing m still help? ===")
print(f"{'n':>6} {'m':>4} {'r4':>7} {'relSE':>7} {'relSE*sqrt(n)':>13}")
for n in [500,1000,2000,4000]:
    m = int(round(2*np.sqrt(n)))   # a bit larger m so LD has room to act
    r4, se = r4_and_relSE(W_epi(Z_LD(n, m, rng, 0.95)))
    print(f"{n:>6} {m:>4} {r4:>7.1f} {se:>7.4f} {se*np.sqrt(n):>13.3f}")
