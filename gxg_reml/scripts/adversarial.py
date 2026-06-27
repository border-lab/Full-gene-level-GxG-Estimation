"""
Adversarial checks of gxg_variance_theory.md.
  1. detection formula: CRLB at h2->0 actually converges to 2 s2e^2/(n(mu2-1)) ?
  2. ordering r4 <= M_e <= rank ?  (I claimed it; test it)
  3. 'LD helps detection, hurts estimation': detection var DOWN, est floor UP as rho rises?
  4. THE BIG ONE: does the MLE empirical variance actually achieve CRLB (<< MoM)?
     If CRLB is loose, the 'switch to REML' recommendation is wrong.
"""
import numpy as np
from numpy.linalg import inv, eigh, eigvalsh
from scipy.stats import norm

def Z_LD(n, m, rng, rho):
    idx=np.arange(m); C=rho**np.abs(idx[:,None]-idx[None,:])
    L=np.linalg.cholesky(C+1e-8*np.eye(m)); pr=rng.uniform(0.05,0.5,m); thr=norm.ppf(pr)
    hap=lambda:(rng.standard_normal((n,m))@L.T<thr).astype(float)
    X=hap()+hap(); X=X[:,X.std(0)>0]; return (X-X.mean(0))/X.std(0)
def W_epi(Z):
    n,m=Z.shape; K=Z@Z.T; D=Z*Z; p=m*(m-1)/2; return (K*K-D@D.T)/(2*p)

def crlb(lam,s2g,s2e):
    d=s2g*lam+s2e; I11=.5*np.sum(lam**2/d**2); I12=.5*np.sum(lam/d**2); I22=.5*np.sum(1/d**2)
    return I22/(I11*I22-I12**2)
def momvar(lam,s2g,s2e,n):
    d=s2g*lam+s2e; trW=lam.sum(); trW2=(lam**2).sum()
    T=np.array([[trW2,trW],[trW,n]]); Cq=2*np.array([[(d**2*lam**2).sum(),(d**2*lam).sum()],[(d**2*lam).sum(),(d**2).sum()]])
    return (inv(T)@inv(T)) if False else (inv(T)@Cq@inv(T))[0,0]

rng=np.random.default_rng(101)

print("=== 1. detection: CRLB as h2->0  vs  2 s2e^2/(n(mu2-1)) ===")
Z=Z_LD(1500,40,rng,0.9); lam=np.clip(eigvalsh(W_epi(Z)),0,None); n=len(lam); mu2=(lam**2).sum()/n
pred=lambda s2e: 2*s2e**2/(n*(mu2-1))
for h2 in (1e-2,1e-3,1e-4,1e-5):
    s2g,s2e=h2,1-h2
    print(f"  h2={h2:.0e}: CRLB={crlb(lam,s2g,s2e):.4e}  pred(2se^4/n(mu2-1))={pred(s2e):.4e}  ratio={crlb(lam,s2g,s2e)/pred(s2e):.3f}")

print("\n=== 2. ordering r4 <= M_e <= rank ? ===")
for rho in (0.0,0.6,0.9,0.97):
    Z=Z_LD(1000,40,rng,rho); lam=np.clip(eigvalsh(W_epi(Z)),0,None); n=len(lam)
    Me=lam.sum()**2/(lam**2).sum(); r4=(lam**2).sum()**2/(lam**4).sum(); rank=int((lam>1e-9*lam.max()).sum())
    ok="OK" if (r4<=Me+1e-6<=rank+1e-6) else "VIOLATED"
    print(f"  rho={rho}:  r4={r4:7.2f}  M_e={Me:7.2f}  rank={rank:5d}   {ok}")

print("\n=== 3. LD helps detection (var down), hurts estimation (floor up)? n=2000,m=40 ===")
for rho in (0.0,0.6,0.9,0.97):
    Z=Z_LD(2000,40,rng,rho); lam=np.clip(eigvalsh(W_epi(Z)),0,None); n=len(lam); mu2=(lam**2).sum()/n
    det=2*(0.99**2)/(n*(mu2-1))                       # detection var at h2~0
    r4=(lam**2).sum()**2/(lam**4).sum(); floor=np.sqrt(2/r4)  # MoM est floor (relSE)
    print(f"  rho={rho}:  detection_var={det:.3e}   est_floor(relSE)={floor:.3f}   M_e={lam.sum()**2/(lam**2).sum():.1f}")

print("\n=== 4. DOES THE MLE ACHIEVE CRLB?  (n=400, m=10, rho=.97, 2500 reps) ===")
Z=Z_LD(400,10,rng,0.97); W=W_epi(Z); n=W.shape[0]
w,Vc=eigh(W); w=np.clip(w,0,None)
s2g,s2e=0.9,0.1
hgrid=np.linspace(1e-4,1-1e-4,400)               # profile over h = s2g/(s2g+s2e)
Dgrid=np.outer(hgrid,w)+(1-hgrid)[:,None]        # (G, n): D_i(h)= h*lam_i+(1-h)
Vfull=s2g*W+s2e*np.eye(n); Lc=np.linalg.cholesky(Vfull)
reps=2500; s2g_mle=np.empty(reps); s2g_mom=np.empty(reps)
trW=w.sum(); trW2=(w**2).sum(); Tmom=np.array([[trW2,trW],[trW,n]]); Tmi=inv(Tmom)
for r in range(reps):
    y=Lc@rng.standard_normal(n)
    a2=(Vc.T@y)**2                               # squared projections onto eigvecs of W
    # profile ML over h:
    chat=(a2[None,:]/Dgrid).sum(1)/n             # (G,)
    nll=n*np.log(chat)+np.log(Dgrid).sum(1)      # -2*profile loglik (up to const)
    k=np.argmin(nll); h=hgrid[k]; c=chat[k]
    s2g_mle[r]=c*h
    # MoM on same y:
    q=np.array([y@(W@y), y@y]); s2g_mom[r]=(Tmi@q)[0]
print(f"  CRLB              = {crlb(w,s2g,s2e):.5f}   (relSE {np.sqrt(crlb(w,s2g,s2e))/s2g:.3f})")
print(f"  MLE empirical var = {s2g_mle.var():.5f}   (relSE {np.sqrt(s2g_mle.var())/s2g:.3f}, mean {s2g_mle.mean():.3f})")
print(f"  MoM analytic var  = {momvar(w,s2g,s2e,n):.5f}   (relSE {np.sqrt(momvar(w,s2g,s2e,n))/s2g:.3f})")
print(f"  MoM empirical var = {s2g_mom.var():.5f}   (relSE {np.sqrt(s2g_mom.var())/s2g:.3f}, mean {s2g_mom.mean():.3f})")
