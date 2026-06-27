"""
GO/NO-GO end-to-end feasibility: simulate y = additive + cis-AxA + e and fit the ACTUAL
3-component AI-REML (V = sA A + sAA W + se I). Checks:
  (i)  E[h2_AA_hat] = h2_AA           (unbiased)
  (ii) SD(h2_AA_hat) ~ detection CRB  SE_pred = sqrt(2 (Gram^{-1})_{AA}), Gram_kl = tr(G_k G_l)
  (iii) empirical power (LRT vs 2-comp null) ~ analytic power.
Scale n=1500, G=50 chosen so Me_within/n ~ 0.03 MATCHES the full-scale operating point
(Me_within=523, N=18k -> 0.029): a faithful replica of the regime position. If this passes, the
N~18k claim is an extrapolation along SE ∝ 1/N (the verified detection scaling).
"""
import numpy as np, sys
rng = np.random.default_rng(1)
n, G, m, maf = 1500, 50, 6, 0.3
h2_A, h2_AA = 0.30, 0.05
reps = int(sys.argv[1]) if len(sys.argv) > 1 else 150
P0 = np.eye(n) - np.ones((n, n)) / n
norm = lambda K: K * (n / np.trace(K))
od = ~np.eye(n, dtype=bool)

# fixed design: additive GRM A, pooled cis AxA kernel W
X = rng.binomial(2, maf, size=(n, G * m)).astype(float); gene = np.repeat(np.arange(G), m)
keep = X.std(0) > 0; X, gene = X[:, keep], gene[keep]
Z = (X - X.mean(0)) / X.std(0)
A = norm(Z @ Z.T / Z.shape[1])
W = norm(sum(norm(P0 @ ((Zg := Z[:, gene == g]) @ Zg.T / Zg.shape[1])**2 @ P0)
             for g in range(G)) / G)
LA = np.linalg.cholesky(A + 1e-8 * np.eye(n)); LW = np.linalg.cholesky(W + 1e-8 * np.eye(n))
I = np.eye(n); one = np.ones(n)

Me_within = np.sqrt(2 / W[od].var())
G3 = np.array([[np.sum(Gi * Gj) for Gj in (A, W, I)] for Gi in (A, W, I)])
se_det = np.sqrt(2 * np.linalg.inv(G3)[1, 1])      # sigma->0 detection CRB
# exact CRB at the true operating point h2_AA=0.05 (the proper benchmark)
Vt = h2_A * A + h2_AA * W + (1 - h2_A - h2_AA) * I
Mi = [np.linalg.inv(Vt) @ Gk for Gk in (A, W, I)]
Ifish = np.array([[0.5 * np.sum(Mi[a] * Mi[b].T) for b in range(3)] for a in range(3)])
se_exact = np.sqrt(np.linalg.inv(Ifish)[1, 1])
print(f"design: n={n}, G={G}, Me_within={Me_within:.1f}, Me_within/n={Me_within/n:.3f}")
print(f"detection CRB (sigma->0) = {se_det:.4f} (implied c={Me_within/(n*se_det):.2f}); "
      f"EXACT CRB at h2=0.05 = {se_exact:.4f}")
se_pred = se_exact

def reml(y, Ks, iters=12):
    s = np.full(len(Ks), y.var() / len(Ks))
    for _ in range(iters):
        V = sum(si * Ki for si, Ki in zip(s, Ks))
        Vi = np.linalg.inv(V + 1e-8 * I)
        Vio = Vi @ one; q = one @ Vio
        Pm = Vi - np.outer(Vio, Vio) / q             # P = Vinv - Vinv 1 (1'Vinv1)^-1 1'Vinv
        Py = Pm @ y
        sc = np.array([0.5 * (Py @ Ki @ Py - np.sum(Pm * Ki)) for Ki in Ks])
        AI = np.array([[0.5 * (Py @ Ki @ Pm @ Kj @ Py) for Kj in Ks] for Ki in Ks])
        step = np.linalg.solve(AI + 1e-8 * np.eye(len(Ks)), sc)
        s = np.clip(s + step, 1e-9, None)
    return s

def loglik(y, s, Ks):
    V = sum(si * Ki for si, Ki in zip(s, Ks))
    Vi = np.linalg.inv(V + 1e-8 * I); _, ld = np.linalg.slogdet(V + 1e-8 * I)
    Vio = Vi @ one; q = one @ Vio; Py = Vi @ y - Vio * ((Vio @ y) / q)
    return -0.5 * (ld + np.log(q) + y @ Py)

ests, rej, floored = [], [], 0
for r in range(reps):
    y = np.sqrt(h2_A) * (LA @ rng.standard_normal(n)) + np.sqrt(h2_AA) * (LW @ rng.standard_normal(n)) \
        + np.sqrt(1 - h2_A - h2_AA) * rng.standard_normal(n)
    y -= y.mean(); y /= y.std()
    s = reml(y, [A, W, I]); ests.append(s[1] / s.sum()); floored += s[1] < 1e-6
    s0 = reml(y, [A, I])
    lr = 2 * (loglik(y, s, [A, W, I]) - loglik(y, np.array([s0[0], 0.0, s0[1]]), [A, W, I]))
    rej.append(lr > 2.706)                           # 0.5 chi2_0 + 0.5 chi2_1, alpha=0.05 one-sided
    if (r + 1) % 50 == 0:
        e = np.array(ests); print(f"  [{r+1}/{reps}] mean={e.mean():.4f} sd={e.std():.4f} "
                                  f"power={np.mean(rej):.2f} floored={floored}", flush=True)

e = np.array(ests)
from scipy.stats import norm as N
analytic_power = 1 - N.cdf(N.ppf(0.95) - h2_AA / se_pred)
print(f"\nRESULT  true h2_AA={h2_AA}  ({reps} reps)")
print(f"  mean h2_AA_hat = {e.mean():.4f}  median = {np.median(e):.4f}  bias = {e.mean()-h2_AA:+.4f}")
print(f"  empirical SD   = {e.std():.4f}   vs EXACT CRB = {se_pred:.4f}  (ratio {e.std()/se_pred:.2f})")
print(f"  reps at boundary (s_AA~0) = {floored}/{reps}  ({100*floored/reps:.0f}%)")
print(f"  empirical power= {np.mean(rej):.2f}   vs analytic = {analytic_power:.2f}")
print("GO" if abs(e.mean()-h2_AA) < 0.012 and 0.8 < e.std()/se_pred < 1.3 else "CHECK")
