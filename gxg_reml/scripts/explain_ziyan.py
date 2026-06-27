"""
Do our variance predictions explain Ziyan's existing simulation SDs?
Load his actual chr1 panels (random & contiguous, m=1000), build the raw epistatic
W = (K.K - DD^T)/(2p) he uses (Fig 3/4), and at each n compute:
  - exact MoM relative SE (eq 6, h2=0.9)  -> compare to observed relative-error SD
  - r4, M_e, rank, and the MoM floor sqrt(2/r4)
Observed SDs from thesis Table 2 / Figs 3-4:
  random  (Fig 3): n=1k..8k -> 0.232, 0.164, 0.106, 0.072
  contig. (Fig 4): n=1k..8k -> 0.324, 0.400, 0.332, 0.416   (biased; SD plateaus)
"""
import numpy as np, pandas as pd
from numpy.linalg import inv, eigvalsh

def load_raw(path, nmax):
    # PLINK .raw: 6 metadata cols then genotype dosages; first row header
    df = pd.read_csv(path, sep=r'\s+', nrows=nmax, engine='c')
    X = df.iloc[:, 6:].to_numpy(dtype=np.float64)
    # impute NA by column mean
    col_mean = np.nanmean(X, axis=0)
    inds = np.where(np.isnan(X)); X[inds] = np.take(col_mean, inds[1])
    return X

def W_raw(Z):
    K = Z @ Z.T; D = Z * Z; m = Z.shape[1]; p = m*(m-1)/2
    return (K*K - D @ D.T) / (2*p)

def momrelse(lam, s2g, s2e, n):
    d = s2g*lam + s2e; trW=lam.sum(); trW2=(lam**2).sum()
    T=np.array([[trW2,trW],[trW,n]]); Cq=2*np.array([[(d**2*lam**2).sum(),(d**2*lam).sum()],[(d**2*lam).sum(),(d**2).sum()]])
    return np.sqrt((inv(T)@Cq@inv(T))[0,0])/s2g

NMAX=8000
obs = {'random':{1000:0.232,2000:0.164,4000:0.106,8000:0.072},
       'contig':{1000:0.324,2000:0.400,4000:0.332,8000:0.416}}
paths={'random':'/tmp/geno_random/chr1_Random1KSNP.raw','contig':'/tmp/geno_contig/chr1_Contiguous1KSNP.raw'}
s2g,s2e=0.9,0.1

for tag in ('random','contig'):
    Xall=load_raw(paths[tag], NMAX)
    print(f"\n=== {tag}  (loaded {Xall.shape}) ===")
    print(f"{'n':>6} {'r4':>9} {'M_e':>9} {'rank':>6} {'relSE_eq6':>10} {'floor√(2/r4)':>13} {'observed_SD':>11}")
    for n in (1000,2000,4000,8000):
        X=Xall[:n]; X=X[:, X.std(0)>0]; Z=(X-X.mean(0))/X.std(0)
        W=W_raw(Z); lam=np.clip(eigvalsh(W),0,None)
        r4=(lam**2).sum()**2/(lam**4).sum(); Me=lam.sum()**2/(lam**2).sum()
        rank=int((lam>1e-9*lam.max()).sum())
        print(f"{n:>6} {r4:>9.1f} {Me:>9.1f} {rank:>6} {momrelse(lam,s2g,s2e,n):>10.4f} "
              f"{np.sqrt(2/r4):>13.4f} {obs[tag][n]:>11.3f}")
