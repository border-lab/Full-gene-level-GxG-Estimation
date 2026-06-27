"""Check the corrected floor: relSE_floor = sqrt(2 / r4), r4 = (trW^2)^2 / trW^4."""
import numpy as np
from numpy.linalg import inv, eigvalsh
from scipy.stats import norm

def Z_LD(n, m, rng, rho=0.97):
    idx = np.arange(m); C = rho ** np.abs(idx[:,None]-idx[None,:])
    L = np.linalg.cholesky(C + 1e-8*np.eye(m))
    p = rng.uniform(0.05,0.5,m); thr = norm.ppf(p)
    hap = lambda: (rng.standard_normal((n,m))@L.T < thr).astype(float)
    X = hap()+hap(); X = X[:, X.std(0)>0]
    return (X - X.mean(0))/X.std(0)

def W_std(Z):
    n,m = Z.shape; i,j = np.triu_indices(m,1)
    H = Z[:,i]*Z[:,j]; mu=H.mean(0); sg=H.std(0); mk=sg>1e-10
    H[:,mk]=(H[:,mk]-mu[mk])/sg[mk]; H[:,~mk]=0.0
    return (H@H.T)/H.shape[1]

def analytic_relSE(lam, s2g, s2e, n):
    d = s2g*lam+s2e; trW=lam.sum(); trW2=(lam**2).sum()
    T=np.array([[trW2,trW],[trW,n]])
    Cq=2*np.array([[(d**2*lam**2).sum(),(d**2*lam).sum()],[(d**2*lam).sum(),(d**2).sum()]])
    Ti=inv(T); return np.sqrt((Ti@Cq@Ti)[0,0])/s2g

rng=np.random.default_rng(3); s2g,s2e=0.9,0.1
print(f"{'cfg':>14} {'n':>6} {'r_eff':>7} {'r4':>6} {'sqrt(2/reff)':>12} {'sqrt(2/r4)':>11} {'relSE_exact':>12}")
for m in (20,60):
    for n in (1600,6400,12800):
        Z=Z_LD(n,m,rng); W=W_std(Z); lam=eigvalsh(W)
        lam=np.clip(lam,0,None)
        trW2=(lam**2).sum(); trW4=(lam**4).sum()
        r_eff=lam.sum()**2/trW2; r4=trW2**2/trW4
        print(f"{'LD m='+str(m):>14} {n:>6} {r_eff:>7.1f} {r4:>6.2f} "
              f"{np.sqrt(2/r_eff):>12.4f} {np.sqrt(2/r4):>11.4f} {analytic_relSE(lam,s2g,s2e,n):>12.4f}")
