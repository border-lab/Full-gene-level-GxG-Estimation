"""
When does the REML/CRLB floor actually appear?
Claim: rank(W) = min(n, p), p = C(m,2). When n > p, there are n-p exact-zero
eigenvalues -> r_star saturates at ~p -> CRLB floors at sqrt(2/p).
So: tiny gene (small m) -> REML floors at sqrt(2/p); large region (p>n) -> REML consistent.
Compare REML floor sqrt(2/p) vs MoM floor sqrt(2/r4) on the SAME LD kernels.
"""
import numpy as np
from numpy.linalg import inv, eigvalsh
from scipy.stats import norm

def Z_LD(n, m, rng, rho):
    idx = np.arange(m); C = rho ** np.abs(idx[:, None] - idx[None, :])
    L = np.linalg.cholesky(C + 1e-8*np.eye(m)); pr = rng.uniform(0.05,0.5,m); thr = norm.ppf(pr)
    hap = lambda: (rng.standard_normal((n, m)) @ L.T < thr).astype(float)
    X = hap()+hap(); X = X[:, X.std(0) > 0]; return (X - X.mean(0)) / X.std(0)

def W_epi(Z):
    n, m = Z.shape; K = Z @ Z.T; D = Z * Z; p = m*(m-1)/2
    return (K*K - D @ D.T) / (2*p)

def crlb(lam, s2g, s2e):
    d = s2g*lam + s2e
    I11=0.5*np.sum(lam**2/d**2); I12=0.5*np.sum(lam/d**2); I22=0.5*np.sum(1/d**2)
    return I22/(I11*I22 - I12**2)

def mom_var(lam, s2g, s2e, n):
    d=s2g*lam+s2e; trW=lam.sum(); trW2=(lam**2).sum()
    T=np.array([[trW2,trW],[trW,n]]); Cq=2*np.array([[(d**2*lam**2).sum(),(d**2*lam).sum()],[(d**2*lam).sum(),(d**2).sum()]])
    return (inv(T)@Cq@inv(T))[0,0]

rng=np.random.default_rng(13); s2g,s2e=0.9,0.1
for m in (8, 20):
    p = m*(m-1)//2
    print(f"\n=== LD rho=.97, m={m} (p={p}); REML floor sqrt(2/p)={np.sqrt(2/p):.4f} ===")
    print(f"{'n':>6} {'rank':>6} {'r_star':>7} {'r4':>7} {'relSE_REML':>10} {'relSE_MoM':>10}")
    for n in (200, 800, 3200, 12800):
        W=W_epi(Z_LD(n,m,rng,0.97)); lam=np.clip(eigvalsh(W),0,None)
        d=s2g*lam+s2e; r_star=np.sum((s2g*lam)**2/d**2)
        r4=(lam**2).sum()**2/(lam**4).sum(); rank=int((lam>1e-9*lam.max()).sum())
        print(f"{n:>6} {rank:>6} {r_star:>7.1f} {r4:>7.2f} "
              f"{np.sqrt(crlb(lam,s2g,s2e))/s2g:>10.4f} {np.sqrt(mom_var(lam,s2g,s2e,n))/s2g:>10.4f}")
