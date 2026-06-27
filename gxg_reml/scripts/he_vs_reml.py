"""
Reconcile: literature says HE (MoM) has ~50% higher SE than REML; I showed 35x VARIANCE.
Claim: the efficiency gap is NOT a constant -- it depends on the spectrum (h2 and LD/skew).
Benign genome-wide-like regime -> ~1.2-1.5x SE (literature); pathological single-locus
high-LD high-h2 -> huge.  Same exact formulas (eq 6 = MoM, eq 7 = CRB).  Additive kernel
(matches the HE-vs-REML literature).
"""
import numpy as np
from numpy.linalg import inv, eigvalsh
from scipy.stats import norm

def Z_indep(n,m,rng):
    p=rng.uniform(0.05,0.5,m); X=rng.binomial(2,p,(n,m)).astype(float)
    X=X[:,X.std(0)>0]; return (X-X.mean(0))/X.std(0)
def Z_LD(n,m,rng,rho):
    idx=np.arange(m); C=rho**np.abs(idx[:,None]-idx[None,:]); L=np.linalg.cholesky(C+1e-8*np.eye(m))
    pr=rng.uniform(0.05,0.5,m); thr=norm.ppf(pr)
    hap=lambda:(rng.standard_normal((n,m))@L.T<thr).astype(float)
    X=hap()+hap(); X=X[:,X.std(0)>0]; return (X-X.mean(0))/X.std(0)
def Kadd(Z): return (Z@Z.T)/Z.shape[1]

def crb(lam,s2g,s2e):
    d=s2g*lam+s2e; I11=.5*np.sum(lam**2/d**2); I12=.5*np.sum(lam/d**2); I22=.5*np.sum(1/d**2)
    return I22/(I11*I22-I12**2)
def momvar(lam,s2g,s2e,n):
    d=s2g*lam+s2e; trK=lam.sum(); trK2=(lam**2).sum()
    T=np.array([[trK2,trK],[trK,n]]); Cq=2*np.array([[(d**2*lam**2).sum(),(d**2*lam).sum()],[(d**2*lam).sum(),(d**2).sum()]])
    return (inv(T)@Cq@inv(T))[0,0]

rng=np.random.default_rng(202)
print(f"{'regime':<34}{'r4':>8}{'M_e':>8}{'h2':>5}{'var ratio':>10}{'SE ratio':>9}")
def row(tag,Z,h2):
    lam=np.clip(eigvalsh(Kadd(Z)),0,None); n=len(lam); s2g,s2e=h2,1-h2
    r4=(lam**2).sum()**2/(lam**4).sum(); Me=lam.sum()**2/(lam**2).sum()
    vr=momvar(lam,s2g,s2e,n)/crb(lam,s2g,s2e)
    print(f"{tag:<34}{r4:>8.1f}{Me:>8.1f}{h2:>5.2f}{vr:>10.2f}{np.sqrt(vr):>9.2f}")

# benign, genome-wide-like: many independent SNPs (m>n), flat-ish spectrum, moderate h2
Zb=Z_indep(2000,8000,rng)
for h2 in (0.25,0.5,0.8): row("indep m=8000 (genome-wide-like)",Zb,h2)
# moderate: fewer independent SNPs
Zm=Z_indep(2000,400,rng)
for h2 in (0.5,0.8): row("indep m=400",Zm,h2)
# regional LD
Zr=Z_LD(2000,200,rng,0.8)
for h2 in (0.5,0.8): row("LD rho=.8 m=200 (regional)",Zr,h2)
# pathological single-locus high-LD high-h2 (Ziyan-like)
Zp=Z_LD(2000,20,rng,0.97)
for h2 in (0.5,0.9): row("LD rho=.97 m=20 (single locus)",Zp,h2)
