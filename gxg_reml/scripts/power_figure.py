"""
Predicted vs observed detection power for the within-gene pooled AxA estimator, under H1
(h2_AA=0.05) and H0 (h2_AA=0), across sample size N. Test = variance-component score test for
sigma^2_AA=0 in V = sigma^2_A A + sigma^2_AA W + sigma^2_e I, additive/residual nuisance at null
plug-in values (locally most powerful). The null covariance V0 is inverted ONCE per N, so each
replicate is a quadratic form.

  U(y) = 1/2 ( y' B y - tr(P0 W) ),  B = P0 W P0,  P0 the null GLS projector.
  H0: E[U]=0, Var(U)=1/2 tr(B V0 B V0); reject (one-sided 5%) if U/sd0 > z_0.95.
  PREDICTED power = analytic non-centrality of this test: lam = E[U|H1]/sd0,
    E[U|H1] = 1/2 ( tr(B V1) - tr(P0 W) ),  power = 1 - Phi(z_0.95 - lam).
  The exact 3-component CRB power is plotted as the efficiency BOUND (an efficient estimator);
  the full REML LRT point from feasibility_sim.py anchors where REML sits between the two.
Writes within_gene_power.{pdf,png}.
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numpy.linalg import inv, cholesky
from scipy.stats import norm

rng = np.random.default_rng(11)
G, m, maf, h2_A, h2_AA = 50, 6, 0.3, 0.30, 0.05
nmax, reps = 4000, 2500
z95 = norm.ppf(0.95)
nrm = lambda K: K * (K.shape[0] / np.trace(K))

X = rng.binomial(2, maf, (nmax, G * m)).astype(float); gene = np.repeat(np.arange(G), m)
keep = X.std(0) > 0; X, gene = X[:, keep], gene[keep]; Z = (X - X.mean(0)) / X.std(0)
Afull = Z @ Z.T / Z.shape[1]
Pc = np.eye(nmax) - np.ones((nmax, nmax)) / nmax
Wfull = sum(nrm(Pc @ ((Zg := Z[:, gene == g]) @ Zg.T / Zg.shape[1])**2 @ Pc) for g in range(G)) / G

def setup(n):
    A = nrm(Afull[:n, :n].copy()); W = nrm(Wfull[:n, :n].copy()); I = np.eye(n); one = np.ones(n)
    V0 = h2_A * A + (1 - h2_A) * I
    V0i = inv(V0); v = V0i @ one; P0 = V0i - np.outer(v, v) / (one @ v)
    B = P0 @ W @ P0; cterm = np.trace(P0 @ W)
    Mv = B @ V0; sd0 = np.sqrt(0.5 * np.sum(Mv * Mv.T))
    V1 = h2_A * A + h2_AA * W + (1 - h2_A - h2_AA) * I
    lam = 0.5 * (np.sum(B * V1) - cterm) / sd0            # E[U|H1]/sd0
    return A, W, I, B, cterm, sd0, lam

def crb_power(A, W, n):
    I = np.eye(n); V = h2_A * A + h2_AA * W + (1 - h2_A - h2_AA) * I
    Mi = [inv(V) @ Gk for Gk in (A, W, I)]
    F = np.array([[0.5 * np.sum(Mi[a] * Mi[b].T) for b in range(3)] for a in range(3)])
    return 1 - norm.cdf(z95 - h2_AA / np.sqrt(inv(F)[1, 1]))

Ns = [500, 800, 1200, 1800, 2500, 3200, 4000]
typeI, powObs, powPred, powCRB = [], [], [], []
for n in Ns:
    A, W, I, B, cterm, sd0, lam = setup(n)
    LA = cholesky(A + 1e-8 * I); LW = cholesky(W + 1e-8 * I); se1 = np.sqrt(1 - h2_A - h2_AA)
    r0 = r1 = 0
    for _ in range(reps):
        y0 = np.sqrt(h2_A) * (LA @ rng.standard_normal(n)) + np.sqrt(1 - h2_A) * rng.standard_normal(n)
        r0 += 0.5 * (y0 @ B @ y0 - cterm) / sd0 > z95
        y1 = np.sqrt(h2_A) * (LA @ rng.standard_normal(n)) + np.sqrt(h2_AA) * (LW @ rng.standard_normal(n)) + se1 * rng.standard_normal(n)
        r1 += 0.5 * (y1 @ B @ y1 - cterm) / sd0 > z95
    typeI.append(r0 / reps); powObs.append(r1 / reps)
    powPred.append(1 - norm.cdf(z95 - lam)); powCRB.append(crb_power(A, W, n))
    print(f"n={n:5d}  typeI={r0/reps:.3f}  power_obs={r1/reps:.3f}  power_pred={powPred[-1]:.3f}  "
          f"CRB_bound={powCRB[-1]:.3f}", flush=True)

ng = np.linspace(450, 4000, 16).astype(int)
pg_pred, pg_crb = [], []
for n in ng:
    A, W, I, B, cterm, sd0, lam = setup(int(n))
    pg_pred.append(1 - norm.cdf(z95 - lam)); pg_crb.append(crb_power(A, W, int(n)))

np.savez("scripts/_power_data.npz", Ns=Ns, typeI=typeI, powObs=powObs, powPred=powPred,
         powCRB=powCRB, ng=ng, pg_pred=pg_pred, pg_crb=pg_crb, reps=reps)

eb = lambda p: np.sqrt(np.array(p) * (1 - np.array(p)) / reps)
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].plot(ng, pg_pred, "-", color="C0", lw=2, label="predicted (score test)")
ax[0].plot(ng, pg_crb, "--", color="C2", lw=1.5, label="efficiency bound (exact CRB)")
ax[0].errorbar(Ns, powObs, yerr=eb(powObs), fmt="o", color="C3", capsize=3, label="observed (sim)")
ax[0].plot(1500, 0.58, "D", color="k", ms=8, label="REML LRT (feasibility$\\_$sim)")
ax[0].axhline(0.8, ls=":", color="gray"); ax[0].set_ylim(0, 1.02)
ax[0].set_xlabel("sample size $N$"); ax[0].set_ylabel("power"); ax[0].set_title("Under $H_1$:  $h^2_{AA}=0.05$")
ax[0].legend(fontsize=8, loc="lower right")
ax[1].errorbar(Ns, typeI, yerr=eb(typeI), fmt="o", color="C3", capsize=3, label="observed (sim)")
ax[1].axhline(0.05, ls="--", color="C0", label="nominal $\\alpha=0.05$")
ax[1].plot(1500, 0.063, "D", color="k", ms=8, label="REML LRT (feasibility$\\_$sim)")
ax[1].set_ylim(0, 0.16); ax[1].set_xlabel("sample size $N$"); ax[1].set_ylabel("type-I error")
ax[1].set_title("Under $H_0$:  $h^2_{AA}=0$"); ax[1].legend(fontsize=8)
fig.suptitle("Within-gene pooled AxA: predicted vs observed (test scale $G{=}50$, $M_e^{\\rm within}{\\approx}45$)")
fig.tight_layout()
fig.savefig("within_gene_power.pdf"); fig.savefig("within_gene_power.png", dpi=150)
print("wrote within_gene_power.{pdf,png}")
