import numpy as np, pandas as pd
from itertools import combinations
from scipy.optimize import brentq
np.seterr(all='ignore')

print("="*70); print("ANALYTIC: Hivert Eq6/7, biobank meta, WGS  (brief sections 8-11)"); print("="*70)
b=1/74000.0; a=2*b
sej=lambda N: np.sqrt(2/(4*N*a+3*N*(N-1)*b**2))
sem=lambda N,K:(lambda n0: np.sqrt((2/(4*n0*a+3*n0*(n0-1)*b**2))/K))(N/K)
N=254679
print(f"Hivert 255k: HE={1/(N*b):.3f}[0.29] REMLjoint={sej(N):.3f}[0.18] meta={sem(N,8):.3f}[0.25] pen={sem(N,8)/sej(N):.2f}[1.40]")
print(f"N for 76% power(SE=.021): joint={brentq(lambda N:sej(N)-.021,1e4,1e10)/1e6:.1f}M[2.8] meta32k={brentq(lambda N:sem(N,N/32000)-.021,1e5,1e10)/1e6:.0f}M[35]")
print("biobank meta:")
for lab,Nt,K in [("5x200k",1_000_000,5),("10x100k",1_000_000,10),("8x250k",2_000_000,8)]:
    print(f"  {lab:<9} joint={sej(Nt):.3f} meta={sem(Nt,K):.3f} pen={sem(Nt,K)/sej(Nt):.1f}x")
print(f"  N for SE.021: joint={brentq(lambda N:sej(N)-.021,1e4,1e10)/1e6:.1f}M meta250k={brentq(lambda N:sem(N,N/250000)-.021,1e5,1e10)/1e6:.0f}M")
print("WGS Me table @500k (joint/meta):")
def seJ(N,Me): bb=1/Me;aa=2*bb; return np.sqrt(2/(4*N*aa+3*N*(N-1)*bb**2))
def seM(N,Me,n0=32000): bb=1/Me;aa=2*bb;K=N/n0; return np.sqrt((2/(4*n0*aa+3*n0*(n0-1)*bb**2))/K)
for Me in (74000,200000,500000,1000000): print(f"  Me={Me:>8}: joint={seJ(5e5,Me):.2f} meta={seM(5e5,Me):.2f}")
print(f"  255k->500k joint: {sej(254679):.3f}->{sej(5e5):.3f} (factor {sej(254679)/sej(5e5):.2f}[1.75])")
print(f"additive-vs-AxA SE ratio sqrt(Me/2)={np.sqrt(74000/2):.0f}[190]")

print("\n"+"="*70); print("SIM: sparse vs infinitesimal  (brief section 4)"); print("="*70)
rng=np.random.default_rng(1); n,M=1200,200
Z=rng.standard_normal((n,M)); Z=(Z-Z.mean(0))/Z.std(0); Sig=Z@Z.T/M
lam=np.linalg.eigvalsh(Sig); trS,trS2=lam.sum(),(lam**2).sum()
Ti=np.linalg.inv(np.array([[trS2,trS],[trS,n]]))
def mvar(s2g,s2e):
    d=s2g*lam+s2e; Cq=2*np.array([[(d**2*lam**2).sum(),(d**2*lam).sum()],[(d**2*lam).sum(),(d**2).sum()]]); return (Ti@Cq@Ti)[0,0]
def mom(y): q=np.array([y@(Sig@y),y@y]); return (Ti@q)[0]
for h2 in (0.01,0.9):
    sd0=np.sqrt(mvar(h2,1-h2)); row=[]
    for pi in (1.0,0.02,0.005):
        tau=np.sqrt(h2/(M*pi)); est=np.empty(1500)
        for r in range(1500):
            mask=rng.random(M)<pi; gc=np.zeros(M); gc[mask]=rng.standard_normal(mask.sum())*tau
            est[r]=mom(Z@gc+rng.standard_normal(n)*np.sqrt(1-h2))
        row.append(est.std()/sd0)
    print(f"  h2={h2}: inflation pi=1/0.02/0.005 = {row[0]:.1f}/{row[1]:.1f}/{row[2]:.1f}x")

print("\n"+"="*70); print("SIM: centering bias  (brief section 6)"); print("="*70)
def load(path,nmax,mc):
    df=pd.read_csv(path,sep=r'\s+',nrows=nmax,engine='c'); X=df.iloc[:,6:6+mc].to_numpy(float)
    cm=np.nanmean(X,0); ii=np.where(np.isnan(X)); X[ii]=np.take(cm,ii[1]); X=X[:,X.std(0)>0]; return (X-X.mean(0))/X.std(0)
def Wraw(Z): K=Z@Z.T;D=Z*Z;m=Z.shape[1];p=m*(m-1)/2; return (K*K-D@D.T)/(2*p)
rng=np.random.default_rng(3)
for tag,path,mc in [('contig m=80','/tmp/geno_contig/chr1_Contiguous1KSNP.raw',80),
                    ('contig m=1000','/tmp/geno_contig/chr1_Contiguous1KSNP.raw',1000),
                    ('random m=1000','/tmp/geno_random/chr1_Random1KSNP.raw',1000)]:
    Z=load(path,2000,mc); n=Z.shape[0]; W=Wraw(Z); P=np.eye(n)-np.ones((n,n))/n; Wc=P@W@P
    trW,trW2,trWc,trWc2=np.trace(W),np.sum(W*W),np.trace(Wc),np.sum(Wc*Wc)
    L=np.linalg.cholesky(0.9*W+0.1*np.eye(n)+1e-8*np.eye(n)); eu=[];ec=[]
    for _ in range(200):
        y=L@rng.standard_normal(n); yc=y-y.mean()
        eu.append(np.linalg.solve(np.array([[trW2,trW],[trW,n]]),[yc@(W@yc),yc@yc])[0])
        ec.append(np.linalg.solve(np.array([[trWc2,trWc],[trWc,n]]),[yc@(Wc@yc),yc@yc])[0])
    print(f"  {tag:<14} offdiag={W[~np.eye(n,dtype=bool)].mean():+.2f} uncentered={np.mean(eu):.2f}(bias{np.mean(eu)-0.9:+.2f}) centered={np.mean(ec):.2f}")

print("\n"+"="*70); print("SIM: gene-burden epistasis  (brief section 12)"); print("="*70)
rng=np.random.default_rng(2); N,G,mg,maf=2500,200,20,0.01
X=rng.binomial(2,maf,size=(N,G*mg)).astype(float); gene=np.repeat(np.arange(G),mg)
keep=X.std(0)>1e-9; X=X[:,keep]; gene=gene[keep]; Z=(X-X.mean(0))/X.std(0)
B=np.column_stack([X[:,gene==g].sum(1) for g in range(G)]); B=B[:,B.std(0)>1e-9]; Bz=(B-B.mean(0))/B.std(0)
P=np.eye(N)-np.ones((N,N))/N; cen=lambda W:P@W@P; nrm=lambda K:K*(N/np.trace(K))
U=np.column_stack([(Z[:,gene==g].sum(1)**2-(Z[:,gene==g]**2).sum(1)) for g in range(G)]); Uc=(U-U.mean(0))/U.std(0); KA=nrm(cen(Uc@Uc.T/U.shape[1]))
pairs=[Z[:,a]*Z[:,b] for g in range(G) for a,b in combinations(np.where(gene==g)[0],2)]
H=np.column_stack(pairs); Hc=(H-H.mean(0))/H.std(0); KB=nrm(cen(Hc@Hc.T/H.shape[1]))
Gb=Bz@Bz.T/Bz.shape[1]; Wb=nrm(cen(Gb*Gb)); Gv=Z@Z.T/Z.shape[1]; Wv=nrm(cen(Gv*Gv))
def mf(W):
    Ai=np.linalg.inv(np.array([[np.sum(W*W),np.trace(W)],[np.trace(W),N]]))
    return lambda y:(Ai@np.array([(y-y.mean())@(W@(y-y.mean())),(y-y.mean())@(y-y.mean())]))[0]
fb,fv=mf(Wb),mf(Wv)
for fc,lab in [(1.0,'coordinated'),(0.5,'50/50'),(0.0,'random pairwise')]:
    L=np.linalg.cholesky(0.3*nrm(fc*KA+(1-fc)*KB)+0.7*np.eye(N)+1e-8*np.eye(N)); eb=[];ev=[]
    for _ in range(200):
        y=L@rng.standard_normal(N); eb.append(fb(y)); ev.append(fv(y))
    print(f"  {lab:<16}: burden={np.mean(eb):.2f}+/-{np.std(eb):.2f}  variant={np.mean(ev):.2f}+/-{np.std(ev):.2f}")
print("\nALL_VERIFY_DONE")
