"""
Verify the CRLB derivation for y ~ N(0, s2g*W + s2e*I):
  (1) eigenvalue-sum CRLB  ==  matrix-trace CRLB  (algebra check)
  (2) CRLB  <=  MoM variance               (REML/ML efficient vs MoM)
  (3) CRLB floor  ~  2*s2g^2 / r_star       (r_star = sum (s2g*lam)^2/d^2)
  (4) CRLB plateaus in n under LD           (REML can't restore consistency)
"""
import numpy as np
from numpy.linalg import inv, eigvalsh
from scipy.stats import norm

def Z_indep(n, m, rng):
    p = rng.uniform(0.05, 0.5, m); X = rng.binomial(2, p, (n, m)).astype(float)
    X = X[:, X.std(0) > 0]; return (X - X.mean(0)) / X.std(0)

def Z_LD(n, m, rng, rho):
    idx = np.arange(m); C = rho ** np.abs(idx[:, None] - idx[None, :])
    L = np.linalg.cholesky(C + 1e-8*np.eye(m)); pr = rng.uniform(0.05,0.5,m); thr = norm.ppf(pr)
    hap = lambda: (rng.standard_normal((n, m)) @ L.T < thr).astype(float)
    X = hap()+hap(); X = X[:, X.std(0) > 0]; return (X - X.mean(0)) / X.std(0)

def W_epi(Z):
    n, m = Z.shape; K = Z @ Z.T; D = Z * Z; p = m*(m-1)/2
    return (K*K - D @ D.T) / (2*p)

def crlb_eig(lam, s2g, s2e):
    d = s2g*lam + s2e
    I11 = 0.5*np.sum(lam**2/d**2); I12 = 0.5*np.sum(lam/d**2); I22 = 0.5*np.sum(1/d**2)
    return I22/(I11*I22 - I12**2)

def crlb_mat(W, s2g, s2e):
    n = W.shape[0]; V = s2g*W + s2e*np.eye(n); Vi = inv(V)
    ViW = Vi @ W
    I11 = 0.5*np.einsum('ij,ji->', ViW, ViW)
    I12 = 0.5*np.einsum('ij,ji->', ViW, Vi)
    I22 = 0.5*np.einsum('ij,ji->', Vi, Vi)
    return I22/(I11*I22 - I12**2)

def mom_var(lam, s2g, s2e, n):
    d = s2g*lam + s2e; trW=lam.sum(); trW2=(lam**2).sum()
    T = np.array([[trW2, trW],[trW, n]])
    Cq = 2*np.array([[(d**2*lam**2).sum(),(d**2*lam).sum()],[(d**2*lam).sum(),(d**2).sum()]])
    Ti = inv(T); return (Ti@Cq@Ti)[0,0]

rng = np.random.default_rng(11); s2g, s2e = 0.9, 0.1

print("(1)+(2)+(3) checks (LD rho=.95, m=120):")
print(f"{'n':>5} {'CRLB_eig':>9} {'CRLB_mat':>9} {'MoM_var':>9} {'2s2g^2/r*':>10} {'r_star':>7} {'rank':>6}")
for n in [400, 800, 1600]:
    W = W_epi(Z_LD(n, 120, rng, 0.95)); lam = np.clip(eigvalsh(W), 0, None)
    d = s2g*lam + s2e; r_star = np.sum((s2g*lam)**2/d**2)
    rank = int(np.sum(lam > 1e-9*lam.max()))
    ce, cm, mv = crlb_eig(lam,s2g,s2e), crlb_mat(W,s2g,s2e), mom_var(lam,s2g,s2e,n)
    print(f"{n:>5} {ce:>9.5f} {cm:>9.5f} {mv:>9.5f} {2*s2g**2/r_star:>10.5f} {r_star:>7.2f} {rank:>6}")

print("\n(4) does the CRLB plateau in n? relSE = sqrt(CRLB)/s2g")
for tag, Zf in [("INDEP m=120", lambda n,r: Z_indep(n,120,r)),
                ("LD rho=.97 m=120", lambda n,r: Z_LD(n,120,r,0.97))]:
    print(f"  {tag}:")
    print(f"    {'n':>6} {'r_star':>7} {'relSE_CRLB':>10} {'relSE_MoM':>10}")
    for n in [400, 800, 1600, 3200]:
        W = W_epi(Zf(n, rng)); lam = np.clip(eigvalsh(W), 0, None)
        d = s2g*lam+s2e; r_star = np.sum((s2g*lam)**2/d**2)
        ce = crlb_eig(lam,s2g,s2e); mv = mom_var(lam,s2g,s2e,n)
        print(f"    {n:>6} {r_star:>7.2f} {np.sqrt(ce)/s2g:>10.4f} {np.sqrt(mv)/s2g:>10.4f}")
