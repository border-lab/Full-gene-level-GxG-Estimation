"""Does the O(nm) sketch error in W actually move the REML estimate?

Simulate y from the TRUE standardized W, then fit exact (dense) REML with
 (a) the true W,
 (b) the frozen-probe sketch W_hat at several probe counts N.
Also report what MoM (scalars only) does with the same probe budget.
"""
import numpy as np
from scipy.optimize import minimize_scalar

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


def std_H(Z):
    m = Z.shape[1]
    i, j = np.triu_indices(m, k=1)
    H = Z[:, i] * Z[:, j]
    mu, sd = H.mean(0), H.std(0)
    ok = sd > 1e-10
    H[:, ok] = (H[:, ok] - mu[ok]) / sd[ok]
    H[:, ~ok] = 0.0
    return H, len(i)


def reml_1comp(y, W):
    """Profile REML for y ~ N(0, s2g W + s2e I) with no fixed effects, h2 grid."""
    lam, Q = np.linalg.eigh(W)
    Qy = Q.T @ y
    n = len(y)

    def negll(h2):
        h2 = min(max(h2, 1e-8), 1 - 1e-8)
        d = h2 * lam + (1 - h2)
        d = np.maximum(d, 1e-12)
        s2 = (Qy**2 / d).sum() / n
        return 0.5 * (np.log(d).sum() + n * np.log(s2) + n)

    r = minimize_scalar(negll, bounds=(1e-6, 1 - 1e-6), method='bounded',
                        options={'xatol': 1e-7})
    return r.x


# ------------------------------------------------------------------- setup
n, m = 400, 60
Z = sim_geno(n, m, seed=5)
Dz = Z * Z
H, p = std_H(Z)
W = H @ H.T / p
h2_true = 0.5

print(f"n={n} m={m} p={p}  tr(W)/n={np.trace(W)/n:.4f}  "
      f"eff.rank={np.trace(W)**2/np.trace(W@W):.1f}")

# phenotype from the TRUE W
lamW, QW = np.linalg.eigh(np.maximum(W, W.T))
L = QW @ np.diag(np.sqrt(np.maximum(lamW, 0)))
g = L @ rng.normal(size=n)
g *= np.sqrt(h2_true) / g.std()
e = rng.normal(size=n)
e *= np.sqrt(1 - h2_true) / e.std()
y = g + e
y -= y.mean()

print(f"\nh2 with the TRUE W : {reml_1comp(y, W):.4f}   (target {h2_true})")


# ------------------------------------------ frozen-probe sketches of W
def W_hat_hutch(Z, N, seed):
    """Hutchinson-diag sketch of the standardized W, materialised (test only)."""
    r = np.random.default_rng(seed)
    n, m = Z.shape
    Gm = Z.T @ Z / n
    Sm = (Z * Z).T @ (Z * Z) / n
    Bf = 1.0 / (Sm - Gm * Gm)
    lam, Q = np.linalg.eigh(Bf)
    Fac, sgn = Q * np.sqrt(np.abs(lam)), np.sign(lam)
    beta = np.einsum('ia,ia,a->i', Fac, Fac, sgn)
    Dz = Z * Z
    q = Fac.shape[1]

    Ah = np.zeros((n, n))
    for _ in range(N):
        v = r.choice([-1.0, 1.0], size=n)
        for a in range(q):
            Kb = (Z * (Fac[:, a] * sgn[a])) @ Z.T
            Kc = (Z * Fac[:, a]) @ Z.T
            Ah += Kb * (v[:, None] * (Kc @ v)[None, :])
    Ah /= N
    Ah = (Ah + Ah.T) / 2
    A = Ah - Dz @ (beta[:, None] * Dz.T)
    vRf = Ah @ np.ones(n) / n                      # NOTE: uses same frozen probes
    v_R = vRf - Dz @ beta
    s_T = np.ones(n) @ vRf / n - beta.sum()
    one = np.ones(n)
    return (A - np.outer(v_R, one) - np.outer(one, v_R)
            + s_T * np.outer(one, one)) / (2 * p)


def W_hat_pair(Z, N, seed):
    """Pair-probe sketch: psi = H w_tilde with w_tilde_(ij) = w_i w_j.

    Exact for the CENTRED-SCALED H by using the true B,R,T corrections; here we
    build it the honest cheap way -- weighted pair probe on the raw products,
    then the exact O(nm) centring corrections.
    """
    r = np.random.default_rng(seed)
    n, m = Z.shape
    Gm = Z.T @ Z / n
    Sm = (Z * Z).T @ (Z * Z) / n
    sig = np.sqrt(Sm - Gm * Gm)
    i, j = np.triu_indices(m, k=1)
    # cheap surrogate: sample pair probes in the m^2 space with rank-1 structure
    Psi = np.zeros((n, N))
    for l in range(N):
        w = r.choice([-1.0, 1.0], size=m)
        # weight each pair by 1/sigma_ij is NOT rank-1; use the separable part
        Psi[:, l] = 0.5 * ((Z @ w) ** 2 - (Z * Z) @ np.ones(m))
    S = Psi @ Psi.T / N
    return S / p


print("\n{:>7} {:>14} {:>10} {:>12}".format("N", "relFro(W)", "h2_hat", "minEV"))
for N in [32, 128, 512]:
    Wh = W_hat_hutch(Z, N, seed=100 + N)
    Wh = (Wh + Wh.T) / 2
    rel = np.linalg.norm(Wh - W) / np.linalg.norm(W)
    print(f"{N:>7} {rel:>14.4f} {reml_1comp(y, Wh):>10.4f} "
          f"{np.linalg.eigvalsh(Wh)[0]:>12.4f}")

# --------------------------------------------- MoM with the same probe budget
print("\nMoM (scalars only) with a pair-probe budget -- same O(nm) primitive")
yWy = y @ W @ y
trW, trW2 = np.trace(W), np.trace(W @ W)
h2_mom_exact = None


def mom(y, trW, trW2, yWy, n):
    """2x2 MoM for (s2g, s2e) with intercept-free model."""
    A = np.array([[trW2, trW], [trW, n]])
    b = np.array([yWy, y @ y])
    s = np.linalg.solve(A, b)
    return s[0] / (s[0] + s[1])


print(f"  exact scalars                    -> h2 = {mom(y, trW, trW2, yWy, n):.4f}")
for N in [32, 128, 512]:
    r = np.random.default_rng(7)
    # pair probes z with E[z z'] = p W  (unstandardized surrogate for the shape)
    zs = []
    for _ in range(N):
        w = r.choice([-1.0, 1.0], size=m)
        zs.append(0.5 * ((Z @ w) ** 2 - (Z * Z).sum(1)))
    Zp = np.array(zs).T / np.sqrt(p)
    trW_h = (Zp**2).sum() / N
    trW2_h = np.linalg.norm(Zp.T @ Zp / N, 'fro')**2 * 0 + ((Zp.T @ Zp / N)**2).sum()
    yWy_h = ((y @ Zp)**2).sum() / N
    print(f"  N={N:>4}: trW {trW_h/trW:6.3f}x  trW2 {trW2_h/trW2:6.3f}x  "
          f"yWy {yWy_h/yWy:6.3f}x  -> h2 = {mom(y, trW_h, trW2_h, yWy_h, n):.4f}")
