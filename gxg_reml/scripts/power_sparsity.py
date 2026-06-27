"""
Figure 1 variation: detection power vs sample size for the within-gene pooled AxA estimator under
DIFFERENT AxA architectures -- varying the causal fraction pi (hence the number of causal interaction
pairs k = pi*P), variance-matched so total h2_AA = 0.05 is held fixed. The analyst's test always uses
the infinitesimal kernel W = HH^T/P; the true signal g = H gamma uses a sparse gamma (Bernoulli(pi)
support, Gaussian magnitudes), so excess kurtosis = 3(1-pi)/pi. Tests the architecture-robustness of
detection: do the power curves move as the variance is concentrated into fewer pairs?

H0 (type-I) is architecture-free (no AxA signal), so one null curve. Reference: infinitesimal analytic
prediction (score-test non-centrality). Writes within_gene_power_sparsity.{pdf,png}.
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numpy.linalg import inv, cholesky
from scipy.stats import norm

import sys
rng = np.random.default_rng(13)
G, maf, h2_A, h2_AA = 50, 0.3, 0.30, 0.05
reps = int(sys.argv[1]) if len(sys.argv) > 1 else 1500   # argv: reps [m] [tag]
m = int(sys.argv[2]) if len(sys.argv) > 2 else 6
TAG = sys.argv[3] if len(sys.argv) > 3 else ""
nmax = 4000
z95 = norm.ppf(0.95)
nrm = lambda K: K * (K.shape[0] / np.trace(K))

# genotypes; additive GRM; within-gene product-feature matrix H (n x P)
X = rng.binomial(2, maf, (nmax, G * m)).astype(float); gene = np.repeat(np.arange(G), m)
keep = X.std(0) > 0; X, gene = X[:, keep], gene[keep]; Z = (X - X.mean(0)) / X.std(0)
Afull = Z @ Z.T / Z.shape[1]
cols = []
for g in range(G):
    Zg = Z[:, gene == g]; mg = Zg.shape[1]
    for a in range(mg):
        for b in range(a + 1, mg):
            p = Zg[:, a] * Zg[:, b]; s = p.std()
            if s > 0:
                cols.append((p - p.mean()) / s)
Hfull = np.column_stack(cols); P = Hfull.shape[1]
print(f"P (within-gene causal pairs available) = {P}", flush=True)

def setup(n):
    A = nrm(Afull[:n, :n].copy()); Hn = Hfull[:n]; W = nrm(Hn @ Hn.T / P)
    I = np.eye(n); one = np.ones(n)
    V0 = h2_A * A + (1 - h2_A) * I; V0i = inv(V0); v = V0i @ one; P0 = V0i - np.outer(v, v) / (one @ v)
    B = P0 @ W @ P0; cterm = np.trace(P0 @ W); Mv = B @ V0; sd0 = np.sqrt(0.5 * np.sum(Mv * Mv.T))
    V1 = h2_A * A + h2_AA * W + (1 - h2_A - h2_AA) * I; lam = 0.5 * (np.sum(B * V1) - cterm) / sd0
    return A, Hn, I, B, cterm, sd0, lam, cholesky(A + 1e-8 * I)

Ns = [500, 800, 1200, 1800, 2500, 3200, 4000]
pis = [1.0, 0.1, 0.02, 0.005]
powers = {pi: [] for pi in pis}; typeI = []; pred_inf = []
for n in Ns:
    A, Hn, I, B, cterm, sd0, lam, LA = setup(n); se1 = np.sqrt(1 - h2_A - h2_AA)
    r0 = 0
    for _ in range(reps):
        y0 = np.sqrt(h2_A) * (LA @ rng.standard_normal(n)) + np.sqrt(1 - h2_A) * rng.standard_normal(n)
        r0 += 0.5 * (y0 @ B @ y0 - cterm) / sd0 > z95
    typeI.append(r0 / reps); pred_inf.append(1 - norm.cdf(z95 - lam))
    for pi in pis:
        scale = np.sqrt(h2_AA / (pi * P)); r1 = 0
        for _ in range(reps):
            g = scale * (Hn @ ((rng.random(P) < pi) * rng.standard_normal(P)))
            y1 = np.sqrt(h2_A) * (LA @ rng.standard_normal(n)) + g + se1 * rng.standard_normal(n)
            r1 += 0.5 * (y1 @ B @ y1 - cterm) / sd0 > z95
        powers[pi].append(r1 / reps)
    print(f"n={n:5d} typeI={typeI[-1]:.3f} " + " ".join(f"pi{pi}={powers[pi][-1]:.3f}" for pi in pis), flush=True)

np.savez(f"scripts/_power_sparsity{TAG}.npz", Ns=Ns, typeI=typeI, pred_inf=pred_inf,
         pis=pis, P=P, reps=reps, m=m, **{f"pow_{pi}": powers[pi] for pi in pis})

eb = lambda p: np.sqrt(np.array(p) * (1 - np.array(p)) / reps)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(Ns, pred_inf, "-", color="gray", lw=1.5, label="infinitesimal (analytic)")
mk = ["o", "s", "^", "v"]
for pi, mm in zip(pis, mk):
    k = int(round(pi * P))
    lab = f"infinitesimal, $k={k}$" if pi == 1.0 else f"$\\pi={pi}$, $k={k}$"
    ax[0].errorbar(Ns, powers[pi], yerr=eb(powers[pi]), fmt=mm + "-", capsize=2, ms=5, label=lab)
ax[0].axhline(0.8, ls=":", color="gray"); ax[0].set_ylim(0, 1.02)
ax[0].set_xlabel("sample size $N$"); ax[0].set_ylabel("power")
ax[0].set_title("Power vs causal fraction $\\pi$ ($h^2_{AA}=0.05$, variance-matched)")
ax[0].legend(fontsize=8, loc="lower right")
ax[1].errorbar(Ns, typeI, yerr=eb(typeI), fmt="o", color="C3", capsize=3, label="observed (any $\\pi$)")
ax[1].axhline(0.05, ls="--", color="C0", label="nominal $\\alpha=0.05$")
ax[1].set_ylim(0, 0.12); ax[1].set_xlabel("sample size $N$"); ax[1].set_ylabel("type-I error")
ax[1].set_title("Under $H_0$ (architecture-free)"); ax[1].legend(fontsize=8)
fig.suptitle(f"Within-gene pooled AxA detection vs architecture (test scale $G{{=}}50$, $m{{=}}{m}$, $P{{=}}{P}$ pairs)")
fig.tight_layout()
fig.savefig(f"within_gene_power_sparsity{TAG}.pdf"); fig.savefig(f"within_gene_power_sparsity{TAG}.png", dpi=150)
print(f"wrote within_gene_power_sparsity{TAG}.{{pdf,png}}")
