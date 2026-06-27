"""
Verify the unified spectral framework for additive & epistatic VC estimation.

Kernel Sigma = (1/M) Phi Phi^T, Phi = n x M standardized features
  additive : features = m SNPs        -> M = m,      R = SNP corr matrix
  epistatic: features = SNP pairs ZaZb -> M = C(m,2), R = interaction-feature corr

Phenotype y ~ N(0, s2g*Sigma + s2e*I).  Eigenvalues lam_i of Sigma, d_i = s2g*lam_i+s2e.
Claims to check:
  (A) r4(Sigma) := (tr S^2)^2 / tr S^4  ==  tr(R^2)^2 / tr(R^4)   [feature-spectrum reduction]
  (B) mu2 := tr(Sigma^2)/n  ~  1 + n/M + n*rho2bar + (kappa-1)/M
  (C) M_e := n/mu2  ~  (1/n + 1/M + rho2bar)^(-1)
  (D) detection var  = 2 s2e^2 / (n (mu2-1))  == exact CRLB at h2->0
  (E) floors (h2=.9, large n):  MoM -> 2 s2g^2/r4 ,  REML -> 2 s2g^2/r_star
  (F) epistatic feature-LD rho2bar_int ~ order (rho2bar_SNP)   (<< 1, decorrelates)
"""
import numpy as np
from numpy.linalg import inv, eigvalsh
from scipy.stats import norm

def Z_LD(n, m, rng, rho):
    idx=np.arange(m); C=rho**np.abs(idx[:,None]-idx[None,:])
    L=np.linalg.cholesky(C+1e-8*np.eye(m)); pr=rng.uniform(0.05,0.5,m); thr=norm.ppf(pr)
    hap=lambda:(rng.standard_normal((n,m))@L.T<thr).astype(float)
    X=hap()+hap(); X=X[:,X.std(0)>0]; return (X-X.mean(0))/X.std(0)

def add_features(Z):           # additive: columns are the SNPs
    return Z / np.sqrt(Z.shape[1])           # Sigma = (1/m) Z Z^T  -> features Z, /sqrt(m) folds the 1/M
def epi_features(Z):           # epistatic: standardized pairwise products
    n,m=Z.shape; i,j=np.triu_indices(m,1)
    H=Z[:,i]*Z[:,j]; mu=H.mean(0); sg=H.std(0); mk=sg>1e-10
    H[:,mk]=(H[:,mk]-mu[mk])/sg[mk]; H[:,~mk]=0.0
    return H/np.sqrt(H.shape[1])

def kernel(Phi): return Phi@Phi.T            # Sigma = Phi Phi^T  (Phi already has 1/sqrt(M))

def spec(Sig):
    lam=np.clip(eigvalsh(Sig),0,None); return lam

def crlb(lam,s2g,s2e):
    d=s2g*lam+s2e; I11=.5*np.sum(lam**2/d**2); I12=.5*np.sum(lam/d**2); I22=.5*np.sum(1/d**2)
    return I22/(I11*I22-I12**2)
def momvar(lam,s2g,s2e,n):
    d=s2g*lam+s2e; trS=lam.sum(); trS2=(lam**2).sum()
    T=np.array([[trS2,trS],[trS,n]]); Cq=2*np.array([[(d**2*lam**2).sum(),(d**2*lam).sum()],[(d**2*lam).sum(),(d**2).sum()]])
    return (inv(T)@Cq@inv(T))[0,0]

rng=np.random.default_rng(17)

# ---- population feature-correlation summaries from a large reference sample ----
def pop_summaries(featfun, m, rho, n_big=40000):
    Zb=Z_LD(n_big,m,rng,rho); Phi=featfun(Zb)*np.sqrt(Phicols(featfun,m))  # unscale 1/sqrt(M) -> unit-var cols
    # standardize columns to unit variance (they already ~are), corr matrix:
    Phi=(Phi-Phi.mean(0))/Phi.std(0)
    R=np.corrcoef(Phi.T)
    M=R.shape[0]
    off=R[np.triu_indices(M,1)]
    rho2bar=np.mean(off**2)
    kappa=np.mean((Phi**4).mean(0))     # E[phi^4] for unit-var standardized features
    return rho2bar, kappa, M
def Phicols(featfun,m):
    return m if featfun is add_features else m*(m-1)//2

print("===== (A) r4 identity:  (trS^2)^2/trS^4  vs  tr(R^2)^2/tr(R^4) =====")
for tag,ff,m in [("ADD",add_features,60),("EPI",epi_features,30)]:
    Z=Z_LD(800,m,rng,0.9); Phi=ff(Z); Sig=kernel(Phi); lam=spec(Sig)
    r4_sigma=(lam**2).sum()**2/(lam**4).sum()
    Pn=(Phi-Phi.mean(0)); Pn=Pn/np.linalg.norm(Pn,axis=0); R=Pn.T@Pn      # feature gram (~corr)
    nu=np.clip(eigvalsh(R),0,None); r4_R=(nu**2).sum()**2/(nu**4).sum()
    print(f"  {tag}: r4(Sigma)={r4_sigma:8.2f}   r4(R)={r4_R:8.2f}")

print("\n===== (B,C) mu2 and effective markers (ADDITIVE) =====")
print(f"{'rho':>5} {'n':>6} {'m':>4} {'mu2_emp':>8} {'mu2_pred':>9} {'Me_emp':>8} {'Me_pred':>8}")
for rho in (0.0, 0.6, 0.9):
    rho2,kap,M = pop_summaries(add_features, 60, rho)
    for n in (500,1500):
        Z=Z_LD(n,60,rng,rho); lam=spec(kernel(add_features(Z)))
        mu2_emp=(lam**2).sum()/n
        mu2_pred=1 + n/M + n*rho2 + (kap-1)/M
        Me_emp=n/mu2_emp; Me_pred=1.0/(1/n+1/M+rho2)
        print(f"{rho:>5} {n:>6} {M:>4} {mu2_emp:>8.3f} {mu2_pred:>9.3f} {Me_emp:>8.1f} {Me_pred:>8.1f}")

print("\n===== (D) detection var 2 s2e^2/(n(mu2-1)) vs exact CRLB at h2=0.01 (ADD, rho=.9) =====")
s2e=0.99; s2g=0.01
for n in (1000,2000,4000):
    Z=Z_LD(n,60,rng,0.9); lam=spec(kernel(add_features(Z))); mu2=(lam**2).sum()/n
    print(f"  n={n:>5}  detect_pred={2*s2e**2/(n*(mu2-1)):.4e}   CRLB={crlb(lam,s2g,s2e):.4e}   MoM={momvar(lam,s2g,s2e,n):.4e}")

print("\n===== (E) floors at h2=0.9 (ADD, rho=.97, large n): MoM=2s2g^2/r4, REML=2s2g^2/r* =====")
s2g,s2e=0.9,0.1
for n in (800,3200):
    Z=Z_LD(n,30,rng,0.97); lam=spec(kernel(add_features(Z))); d=s2g*lam+s2e
    r4=(lam**2).sum()**2/(lam**4).sum(); rstar=np.sum((s2g*lam)**2/d**2)
    print(f"  n={n:>5}  relSE_MoM={np.sqrt(momvar(lam,s2g,s2e,n))/s2g:.4f} (floor {np.sqrt(2/r4):.4f})"
          f"   relSE_REML={np.sqrt(crlb(lam,s2g,s2e))/s2g:.4f} (floor {np.sqrt(2/rstar):.4f})")

print("\n===== (F) feature LD: epistatic rho2bar vs additive rho2bar (same SNPs) =====")
for rho in (0.6, 0.9):
    rho2_add,_,Ma = pop_summaries(add_features, 30, rho)
    rho2_epi,_,Me = pop_summaries(epi_features, 30, rho)
    print(f"  SNP rho={rho}:  rho2bar_SNP={rho2_add:.4f} (M={Ma})   rho2bar_INT={rho2_epi:.4f} (M={Me})"
          f"   ratio INT/SNP^2 ~ {rho2_epi/rho2_add**2:.2f}")
