"""
Adversarial stress-test of the feasibility claim. The headline rides on two extrapolations and a
power claim; attack all three:
  (1) N-scaling: is SE ∝ 1/N? (extrapolating n=1500 -> ~25k). Exact CRB across N, + REML tracking.
  (2) G-scaling: is Me_within = sqrt(G)*Me_gene as G -> large? Does between-gene LD break it?
  (3) Null calibration: LRT type-I error at h2_AA=0 (is the 0.58 power real or inflated?).
If any fails, the N~25k headline is wrong.
"""
import numpy as np
from numpy.linalg import inv, slogdet
from scipy.stats import norm as Nz
rng = np.random.default_rng(3)
norm = lambda K: K * (K.shape[0] / np.trace(K))

# ---------- (1) N-scaling: exact CRB ----------
Gn, m, maf, h2_A = 50, 6, 0.3, 0.30
nmax = 3000
X = rng.binomial(2, maf, (nmax, Gn * m)).astype(float); gene = np.repeat(np.arange(Gn), m)
keep = X.std(0) > 0; X, gene = X[:, keep], gene[keep]
Z = (X - X.mean(0)) / X.std(0)
Afull = Z @ Z.T / Z.shape[1]
Pn = np.eye(nmax) - np.ones((nmax, nmax)) / nmax
Wfull = sum(norm(Pn @ ((Zg := Z[:, gene == g]) @ Zg.T / Zg.shape[1])**2 @ Pn) for g in range(Gn)) / Gn

def kernels(n):
    return norm(Afull[:n, :n].copy()), norm(Wfull[:n, :n].copy())

def crb_se(n, h2AA):
    A, W = kernels(n); I = np.eye(n)
    V = h2_A * A + h2AA * W + (1 - h2_A - h2AA) * I
    Mi = [inv(V) @ Gk for Gk in (A, W, I)]
    F = np.array([[0.5 * np.sum(Mi[a] * Mi[b].T) for b in range(3)] for a in range(3)])
    return np.sqrt(inv(F)[1, 1])

print("=== (1) N-scaling: exact CRB at h2_AA=0.05 (SE*N should be constant if SE ∝ 1/N) ===", flush=True)
print(f"{'n':>6} {'CRB SE':>9} {'CRB*n':>8}", flush=True)
for n in (500, 800, 1200, 1800, 2500, 3000):
    se = crb_se(n, 0.05); print(f"{n:6d} {se:9.4f} {se*n:8.1f}", flush=True)

# ---------- (2) G-scaling: Me_within = sqrt(G)*Me_gene ----------
print("\n=== (2) G-scaling: Me_within vs sqrt(G)*Me_gene (independent vs between-gene LD) ===", flush=True)
ng, mg = 1500, 6
od = ~np.eye(ng, dtype=bool)
Me = lambda W: np.sqrt(2 / W[od].var())

def gene_kernels(Gtot, between_rho):
    # Gtot genes, mg SNPs each. between_rho>0: AR(1) correlation ACROSS all SNPs (genes leak into neighbours).
    if between_rho == 0:
        Zall = rng.binomial(2, maf, (ng, Gtot * mg)).astype(float)
    else:
        Msnp = Gtot * mg
        idx = np.arange(Msnp);
        # cheap AR(1) latent via cumulative blend, then threshold to genotypes
        L = rng.standard_normal((ng, Msnp))
        for j in range(1, Msnp):
            L[:, j] = between_rho * L[:, j - 1] + np.sqrt(1 - between_rho**2) * L[:, j]
        thr = Nz.ppf(maf)
        Zall = (L < thr).astype(float) + (rng.standard_normal((ng, Msnp)) < thr).astype(float)
    Zall = Zall[:, Zall.std(0) > 0]
    g = np.repeat(np.arange(Gtot), mg)[:Zall.shape[1]]
    Zs = (Zall - Zall.mean(0)) / Zall.std(0)
    return [norm((Zs[:, g == k] @ Zs[:, g == k].T / max(1, (g == k).sum()))**2) for k in range(Gtot)]

for rho_b, tag in [(0.0, "independent"), (0.3, "between-LD rho=.3")]:
    Kg = gene_kernels(1600, rho_b)
    Me_gene = np.mean([Me(K) for K in Kg[:50]])
    print(f"  {tag}: Me_gene={Me_gene:.2f}", flush=True)
    print(f"    {'G':>5} {'Me_within':>10} {'sqrt(G)*Me_gene':>15} {'ratio':>7}", flush=True)
    for G in (25, 100, 400, 1600):
        W = norm(sum(Kg[:G]) / G); mw = Me(W)
        print(f"    {G:5d} {mw:10.2f} {np.sqrt(G)*Me_gene:15.2f} {mw/(np.sqrt(G)*Me_gene):7.3f}", flush=True)

# ---------- (3) REML tracking at a 2nd N, and NULL type-I calibration ----------
def reml(y, Ks, n, iters=12):
    one = np.ones(n); I = np.eye(n); s = np.full(len(Ks), y.var() / len(Ks))
    for _ in range(iters):
        V = sum(si * Ki for si, Ki in zip(s, Ks)); Vi = inv(V + 1e-8 * I)
        Vio = Vi @ one; Pm = Vi - np.outer(Vio, Vio) / (one @ Vio); Py = Pm @ y
        sc = np.array([0.5 * (Py @ Ki @ Py - np.sum(Pm * Ki)) for Ki in Ks])
        AI = np.array([[0.5 * (Py @ Ki @ Pm @ Kj @ Py) for Kj in Ks] for Ki in Ks])
        s = np.clip(s + np.linalg.solve(AI + 1e-8 * np.eye(len(Ks)), sc), 1e-9, None)
    return s

def loglik(y, s, Ks, n):
    one = np.ones(n); I = np.eye(n); V = sum(si * Ki for si, Ki in zip(s, Ks))
    Vi = inv(V + 1e-8 * I); _, ld = slogdet(V + 1e-8 * I)
    Vio = Vi @ one; q = one @ Vio; Py = Vi @ y - Vio * ((Vio @ y) / q)
    return -0.5 * (ld + np.log(q) + y @ Py)

print("\n=== (3a) REML tracks exact CRB at a 2nd sample size (n=2500, h2_AA=0.05) ===", flush=True)
n = 2500; A, W = kernels(n); LA = np.linalg.cholesky(A + 1e-8 * np.eye(n)); LW = np.linalg.cholesky(W + 1e-8 * np.eye(n))
est = []
for _ in range(80):
    y = np.sqrt(0.30) * (LA @ rng.standard_normal(n)) + np.sqrt(0.05) * (LW @ rng.standard_normal(n)) \
        + np.sqrt(0.65) * rng.standard_normal(n); y -= y.mean(); y /= y.std()
    s = reml(y, [A, W, I := np.eye(n)], n); est.append(s[1] / s.sum())
est = np.array(est); cr = crb_se(2500, 0.05)
print(f"  empirical SD={est.std():.4f}  exact CRB={cr:.4f}  ratio={est.std()/cr:.2f}  mean={est.mean():.4f}", flush=True)

print("\n=== (3b) NULL type-I: LRT at h2_AA=0 should reject ~5% (crit 2.706) ===", flush=True)
n = 1500; A, W = kernels(n); LA = np.linalg.cholesky(A + 1e-8 * np.eye(n)); I = np.eye(n)
rej = 0; R = 300
for r in range(R):
    y = np.sqrt(0.30) * (LA @ rng.standard_normal(n)) + np.sqrt(0.70) * rng.standard_normal(n)
    y -= y.mean(); y /= y.std()
    s = reml(y, [A, W, I], n); s0 = reml(y, [A, I], n)
    lr = 2 * (loglik(y, s, [A, W, I], n) - loglik(y, np.array([s0[0], 0.0, s0[1]]), [A, W, I], n))
    rej += lr > 2.706
    if (r + 1) % 100 == 0:
        print(f"  [{r+1}/{R}] type-I = {rej/(r+1):.3f}", flush=True)
print(f"  NULL type-I error = {rej/R:.3f}  (target 0.05; >0.08 => power was inflated)", flush=True)

