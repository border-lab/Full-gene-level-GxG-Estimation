"""
Is the SE-plateau an 'epistasis' problem or a 'kernel effective-rank' problem?
Compare additive kernel K = ZZ'/m vs epistatic kernel W = (K0.K0 - DD')/(2p)
on the SAME genotypes. For each: effective rank of the SQUARED kernel
r4 = (tr S^2)^2 / tr S^4, and the exact relative SE of the variance-component
estimate (validated formula). Component h2 = 0.5 for whichever component is estimated.
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
    p = rng.uniform(0.05, 0.5, m); thr = norm.ppf(p)
    hap = lambda: (rng.standard_normal((n, m)) @ L.T < thr).astype(float)
    X = hap() + hap(); X = X[:, X.std(0) > 0]
    return (X - X.mean(0)) / X.std(0)

def K_add(Z):
    n, m = Z.shape
    return (Z @ Z.T) / m

def W_epi(Z):
    n, m = Z.shape
    K = Z @ Z.T; D = Z * Z
    p = m * (m - 1) / 2
    return (K * K - D @ D.T) / (2 * p)

def relSE(S, s2, s2e):
    """exact Gaussian exact-trace MoM relative SE of the s2 component, from kernel S."""
    n = S.shape[0]
    lam = np.clip(eigvalsh(S), 0, None)
    d = s2 * lam + s2e
    trS, trS2 = lam.sum(), (lam**2).sum()
    T = np.array([[trS2, trS], [trS, n]])
    Cq = 2 * np.array([[(d**2*lam**2).sum(), (d**2*lam).sum()],
                       [(d**2*lam).sum(),    (d**2).sum()]])
    Ti = inv(T)
    var = (Ti @ Cq @ Ti)[0, 0]
    r4 = trS2**2 / (lam**4).sum()
    return np.sqrt(var) / s2, r4

def run(tag, Zfun, m, ns, seed=2):
    rng = np.random.default_rng(seed)
    print(f"\n=== {tag} (m={m}), estimating component with h2=0.5 ===")
    print(f"{'n':>6} | {'r4(K_add)':>9} {'relSE_add':>9} | {'r4(W_epi)':>9} {'relSE_epi':>9}")
    for n in ns:
        Z = Zfun(n, m, rng)
        se_a, r4a = relSE(K_add(Z), 0.5, 0.5)
        se_e, r4e = relSE(W_epi(Z), 0.5, 0.5)
        print(f"{n:>6} | {r4a:>9.1f} {se_a:>9.4f} | {r4e:>9.2f} {se_e:>9.4f}")

if __name__ == "__main__":
    ns = [500, 1000, 2000, 4000, 8000]
    run("INDEPENDENT SNPs", Z_indep, m=200, ns=ns)
    run("MODERATE LD (rho=.9)", lambda n,m,r: Z_LD(n,m,r,0.9), m=200, ns=ns)
    run("STRONG LD (rho=.97)", lambda n,m,r: Z_LD(n,m,r,0.97), m=200, ns=ns)
    # additive with FEW effective SNPs (founder/low-rank-like) -> does additive ALSO plateau?
    run("INDEP but only m=15 SNPs", Z_indep, m=15, ns=ns)
