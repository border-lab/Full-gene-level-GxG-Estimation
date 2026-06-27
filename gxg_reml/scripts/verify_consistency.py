"""
Verify the claims about the gene-level GxG MoM estimator's sampling variance.

Claims under test (all for y ~ N(0, s2g*W + s2e*I), exact-trace MoM):
  (E) Exact Gaussian MoM variance: Cov(theta_hat) = Tinv @ Cov(q) @ Tinv,
      with Cov(q)_{ij} = 2 tr(V G_i V G_j), G1=W, G2=I.
  (A) Regime-A approx:  Var(s2g_hat) ~ 2*sigma_y^4 / trW2_tilde,
      trW2_tilde = trW2 - (trW)^2/n.
  (B) Floor (Regime B):  Var(s2g_hat) -> 2*s2g^2 / r_eff,
      r_eff = (trW)^2 / trW2.
We check (E) against direct Monte-Carlo of the estimator, then check
(A),(B) against (E) across an n-sweep, for independent vs LD genotypes.
"""
import numpy as np
from numpy.linalg import inv, cholesky, eigvalsh
from scipy.stats import norm

# ---------- genotype generators ----------
def Z_indep(n, m, rng):
    p = rng.uniform(0.05, 0.5, m)
    X = rng.binomial(2, p, size=(n, m)).astype(float)
    s = X.std(0); keep = s > 0
    X = X[:, keep]
    return (X - X.mean(0)) / X.std(0)

def Z_LD(n, m, rng, rho=0.97):
    """AR(1) LD haplotypes -> genotypes in {0,1,2} with LD decaying along position."""
    idx = np.arange(m)
    C = rho ** np.abs(idx[:, None] - idx[None, :])
    L = cholesky(C + 1e-8 * np.eye(m))
    p = rng.uniform(0.05, 0.5, m)
    thr = norm.ppf(p)
    def hap():
        g = rng.standard_normal((n, m)) @ L.T
        return (g < thr).astype(float)
    X = hap() + hap()
    s = X.std(0); keep = s > 0
    X = X[:, keep]
    return (X - X.mean(0)) / X.std(0)

def W_standardized(Z):
    """Adjusted/standardized epistatic GRM: columns of H standardized, W = H H^T / p."""
    n, m = Z.shape
    i, j = np.triu_indices(m, 1)
    H = Z[:, i] * Z[:, j]
    mu = H.mean(0); sg = H.std(0)
    mask = sg > 1e-10
    H[:, mask] = (H[:, mask] - mu[mask]) / sg[mask]
    H[:, ~mask] = 0.0
    p = H.shape[1]
    return (H @ H.T) / p

# ---------- variance predictions ----------
def analytic_var(lam, s2g, s2e, n):
    """Exact Gaussian exact-trace MoM variance of s2g_hat from eigenvalues lam of W (length n)."""
    d = s2g * lam + s2e
    trW, trW2 = lam.sum(), (lam**2).sum()
    T = np.array([[trW2, trW], [trW, n]])
    Cq = 2 * np.array([[(d**2 * lam**2).sum(), (d**2 * lam).sum()],
                       [(d**2 * lam).sum(),    (d**2).sum()]])
    Ti = inv(T)
    return (Ti @ Cq @ Ti)[0, 0]

def mc_var(W, s2g, s2e, reps, rng):
    """Brute-force: sample y, run exact MoM, return empirical Var(s2g_hat)."""
    n = W.shape[0]
    V = s2g * W + s2e * np.eye(n)
    L = cholesky(V)
    trW = np.trace(W); trW2 = np.einsum('ij,ji->', W, W)
    Ti = inv(np.array([[trW2, trW], [trW, n]]))
    G = rng.standard_normal((n, reps))
    Y = L @ G                      # columns ~ N(0,V)
    WY = W @ Y
    q1 = np.einsum('ir,ir->r', Y, WY)   # y'Wy
    q2 = np.einsum('ir,ir->r', Y, Y)    # y'y
    s2g_hat = Ti[0, 0] * q1 + Ti[0, 1] * q2
    return s2g_hat.var(), s2g_hat.mean()

# ---------- run ----------
def summarize(tag, Zfun, m, ns, s2g=0.9, s2e=0.1, seed=1, mc_upto=1600, reps=8000):
    sy2 = s2g + s2e
    rng = np.random.default_rng(seed)
    print(f"\n=== {tag} (m={m}, s2g={s2g}) ===")
    print(f"{'n':>6} {'r_eff':>8} {'floor√(2/reff)':>14} {'relSE_exact':>12} "
          f"{'relSE_MC':>10} {'relSE_A':>9} {'mean_hat':>9}")
    for n in ns:
        Z = Zfun(n, m, rng)
        W = W_standardized(Z)
        lam = eigvalsh(W)
        trW, trW2 = lam.sum(), (lam**2).sum()
        r_eff = trW**2 / trW2
        trW2_tilde = trW2 - trW**2 / n
        v_exact = analytic_var(lam, s2g, s2e, n)
        v_A = 2 * sy2**2 / trW2_tilde
        floor = np.sqrt(2.0 / r_eff)                      # relative SE floor (Regime B)
        relSE_exact = np.sqrt(v_exact) / s2g
        relSE_A = np.sqrt(v_A) / s2g
        if n <= mc_upto:
            vmc, mhat = mc_var(W, s2g, s2e, reps, rng)
            mc_str, mean_str = f"{np.sqrt(vmc)/s2g:10.4f}", f"{mhat:9.4f}"
        else:
            mc_str, mean_str = f"{'--':>10}", f"{'--':>9}"
        print(f"{n:>6} {r_eff:>8.1f} {floor:>14.4f} {relSE_exact:>12.4f} "
              f"{mc_str} {relSE_A:>9.4f} {mean_str}")

if __name__ == "__main__":
    # Config 1: large m -> r_eff huge -> should see Regime-A DECLINE (consistent) for indep,
    #           and LD with same m has smaller r_eff -> higher floor / slower decline.
    summarize("INDEP, large m", Z_indep, m=60, ns=[100,200,400,800,1600,3200])
    summarize("LD rho=.97, large m", Z_LD,  m=60, ns=[100,200,400,800,1600,3200])
    # Config 2: small m -> r_eff small -> PLATEAU visible even for independent SNPs
    #           (proves mechanism is effective rank, not LD per se).
    summarize("INDEP, small m", Z_indep, m=20, ns=[200,400,800,1600,3200,6400])
    summarize("LD rho=.97, small m", Z_LD, m=20, ns=[200,400,800,1600,3200,6400])
