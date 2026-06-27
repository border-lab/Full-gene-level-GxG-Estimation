"""
Adversarial check of the closed-form SE results added to the doc:
  (a) Hadamard factor:  Var(offdiag of G.G) ~ 2 Var(G_ij)^2   [G additive GRM]
  (b) additive detection SE  ~ sqrt(2 M_e)/N      (M_e = 1/Var(G_ij))
  (c) AxA detection SE       ~ M_e/N
Both SEs computed EXACTLY from the kernel eigenvalues (h2->0), intercept projected,
and compared to the closed forms.  Also check whether LD/non-Gaussianity breaks (a).
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

def momvar(lam,s2g,s2e,n):
    d=s2g*lam+s2e; trK=lam.sum(); trK2=(lam**2).sum()
    T=np.array([[trK2,trK],[trK,n]]); Cq=2*np.array([[(d**2*lam**2).sum(),(d**2*lam).sum()],[(d**2*lam).sum(),(d**2).sum()]])
    return (inv(T)@Cq@inv(T))[0,0]

def center_norm(K,n):
    P=np.eye(n)-np.ones((n,n))/n
    Kc=P@K@P; Kc*=n/np.trace(Kc); return Kc

rng=np.random.default_rng(7)
print(f"{'geno':<12}{'n':>5}{'m':>6}{'M_e':>8}{'VarGAA/2VarG^2':>15}"
      f"{'SEadd pred':>11}{'SEadd obs':>10}{'SEaxa pred':>11}{'SEaxa obs':>10}")
for tag,Zf in [('indep',lambda n,m: Z_indep(n,m,rng)),
               ('LD .6',lambda n,m: Z_LD(n,m,rng,0.6)),
               ('LD .9',lambda n,m: Z_LD(n,m,rng,0.9))]:
    n,m=2500,1500
    Z=Zf(n,m); n,m=Z.shape
    G=Z@Z.T/m
    od=~np.eye(n,dtype=bool)
    v=G[od].var(); Me=1/v
    Q=G*G
    vAA=Q[od].var()
    lamG=eigvalsh(center_norm(G,n)); lamQ=eigvalsh(center_norm(Q,n))
    SEadd=np.sqrt(momvar(lamG,1e-9,1.0,n)); SEaxa=np.sqrt(momvar(lamQ,1e-9,1.0,n))
    print(f"{tag:<12}{n:>5}{m:>6}{Me:>8.0f}{vAA/(2*v**2):>15.2f}"
          f"{np.sqrt(2*Me)/n:>11.4f}{SEadd:>10.4f}{Me/n:>11.4f}{SEaxa:>10.4f}")
