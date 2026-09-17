"""Which O(nm) stochastic estimator of the Hadamard term is usable?

Two candidates for   A = K_b .* K_c ,   K_b = Z diag(b) Z^T :

 (H) Hutchinson-diag :  A u = diag(K_b D_u K_c) ~ mean_l [ (K_b (u .* K_c v_l)) .* v_l ]
     -- sketches the n-dimension.  Asymmetric at finite N, indefinite.

 (P) Pair-probe      :  psi(w) = (Z(b.*w)) .* (Z(c.*w)),  w Rademacher in R^m
     E[psi psi'] = (K_{b^2}.*K_{c^2}) + (K_{bc}.*K_{bc}) + (D(bc))(D(bc))' - 2 D D_{b^2c^2} D'
     -- sketches the m^2 pair-space.  Symmetric and PSD by construction.
"""
import numpy as np

rng = np.random.default_rng(0)


def sim_geno(n, m, rho=0.85, seed=1):
    r = np.random.default_rng(seed)
    L = r.normal(size=(n, m))
    for j in range(1, m):
        L[:, j] = rho * L[:, j - 1] + np.sqrt(1 - rho**2) * L[:, j]
    f = r.uniform(0.05, 0.5, m)
    thr = np.array([np.quantile(L[:, j], [(1 - f[j])**2, 1 - f[j]**2]) for j in range(m)])
    Gc = (L > thr[:, 0]).astype(float) + (L > thr[:, 1]).astype(float)
    sd = Gc.std(0)
    return (Gc - Gc.mean(0)) / np.where(sd < 1e-12, 1.0, sd)


def dense_W(Z):
    n, m = Z.shape
    i, j = np.triu_indices(m, k=1)
    H = Z[:, i] * Z[:, j]
    mu, sd = H.mean(0), H.std(0)
    ok = sd > 1e-10
    H[:, ok] = (H[:, ok] - mu[ok]) / sd[ok]
    H[:, ~ok] = 0.0
    return H @ H.T / len(i), len(i)


n, m = 300, 60
Z = sim_geno(n, m, seed=5)
Dz = Z * Z
Wt, p = dense_W(Z)
K = Z @ Z.T

# --------------------------------------------------------- unstandardized W0
# W0 = 1/(2p) (K.*K - D D')   -- centring/scaling dropped
W0 = (K * K - Dz @ Dz.T) / (2 * p)

print("=" * 78)
print(f"n={n}  m={m}  p={p}     ||W_std - W_unstd||_F / ||W_std||_F = "
      f"{np.linalg.norm(Wt-W0)/np.linalg.norm(Wt):.3f}")
ev = np.linalg.eigvalsh(Wt)
print(f"  W spectrum: min={ev[0]:.4f} max={ev[-1]:.3f} "
      f"eff.rank tr(W)^2/tr(W^2)={ev.sum()**2/ (ev**2).sum():.1f}  (n={n})")
print("=" * 78)


# ------------------------------------------------------------- estimator (H)
def est_H(Z, b, c, N, rng):
    """Frozen-probe Hutchinson-diag operator, materialised as a matrix."""
    n, m = Z.shape
    Kb = (Z * b) @ Z.T
    Kc = (Z * c) @ Z.T
    Ahat = np.zeros((n, n))
    for _ in range(N):
        v = rng.choice([-1.0, 1.0], size=n)
        Ahat += Kb * (v[:, None] * (Kc @ v)[None, :])
    Ahat /= N
    return (Ahat + Ahat.T) / 2                      # symmetrised


# ------------------------------------------------------------- estimator (P)
def est_P(Z, b, N, rng):
    """Pair-probe operator for the SYMMETRIC case b = c, materialised.

    E[psi psi'] = 2 (K_{b^2}.*K_{b^2}) + (D b^2)(D b^2)' - 2 D D_{b^4} D'
    so   K_{b^2}.*K_{b^2} = 1/2 [ E psi psi' - (D b^2)(D b^2)' + 2 D D_{b^4} D' ].
    """
    n, m = Z.shape
    Dz = Z * Z
    Psi = np.zeros((n, N))
    for l in range(N):
        w = rng.choice([-1.0, 1.0], size=m)
        Psi[:, l] = (Z @ (b * w)) ** 2
    S = Psi @ Psi.T / N
    r1 = Dz @ (b * b)
    corr = Dz @ ((b**4)[:, None] * Dz.T)
    return 0.5 * (S - np.outer(r1, r1) + 2 * corr)


bb = np.ones(m)
target = K * K                                        # b = c = 1
print("\nUNSTANDARDIZED quartic  K.*K   (b = c = 1)")
print(f"{'N':>7} {'(H) relFro':>12} {'(P) relFro':>12} {'(H) minEV of W0hat':>20} {'(P) minEV':>12}")
for N in [16, 64, 256, 1024]:
    AH = est_H(Z, bb, bb, N, np.random.default_rng(11))
    AP = est_P(Z, bb, N, np.random.default_rng(11))
    eH = np.linalg.norm(AH - target) / np.linalg.norm(target)
    eP = np.linalg.norm(AP - target) / np.linalg.norm(target)
    WH = (AH - Dz @ Dz.T) / (2 * p)
    WP = (AP - Dz @ Dz.T) / (2 * p)
    print(f"{N:>7} {eH:>12.4f} {eP:>12.4f} {np.linalg.eigvalsh(WH)[0]:>20.4f}"
          f" {np.linalg.eigvalsh(WP)[0]:>12.4f}")

print(f"\n  reference: true W0 min eigenvalue = {np.linalg.eigvalsh(W0)[0]:.4f}")

# --------------------------------------- how well does a sketch CAN do at all
print("\nBest possible rank-N approximation of the pair-feature matrix H (SVD):")
i, j = np.triu_indices(m, k=1)
H = Z[:, i] * Z[:, j]
mu, sd = H.mean(0), H.std(0)
ok = sd > 1e-10
H[:, ok] = (H[:, ok] - mu[ok]) / sd[ok]
sv = np.linalg.svd(H, compute_uv=False)
tot = (sv**2).sum()
for N in [16, 64, 256, 1024]:
    keep = (sv[:min(N, len(sv))]**2).sum()
    print(f"   N={N:>5}   optimal rank-N captures {keep/tot:6.2%} of ||W||_*  "
          f"-> irreducible relFro >= {np.sqrt(1-keep/tot)*0 + np.linalg.norm(np.diag(sv[min(N,len(sv)):])) / np.linalg.norm(sv):.4f}")
