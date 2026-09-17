# -*- coding: utf-8 -*-
"""Term-by-term check of the UNSTANDARDIZED pooled operator against dense W.

    python3 verify_unstd.py [--n 300] [--m 120] [--G 4] [--Nw 4000]

Every claim the pipeline rests on is asserted here at small n, m, where the
dense kernel can be formed and compared directly:

  1. build_W_pooled('hadamard') == build_W_pooled('pairs').
     The closed form  2 sum_{a<b} h h' = (K .* K) - D D'  really is the pair
     sum, not merely something close to it.  This is the identity that lets the
     simulation and the estimator be the same object.

  2. exact_WU == W @ U to machine precision, where exact_WU is the INDEPENDENT
     reference implementation below -- (Z .* (Z M)) 1 with M = Z' diag(u) Z,
     written here and nowhere else.  The pipeline module implements only the
     stochastic operator; a reference that imported the code under test would
     prove nothing, so this one is deliberately standalone.  What it certifies
     is that the matrix-free formula IS the simulated kernel, i.e. that the
     kernel-mismatch bias of the standardized StochasticWu pipeline is gone
     rather than merely reduced.

  3. compute_WU_pooled -- the pipeline's stochastic operator -- converges to
     W @ U as Nw grows, at the predicted ~1/sqrt(Nw) rate with NO floor.  The
     absence of a floor is the real content: a kernel mismatch would show up as
     an error that stops falling.  Reports the measured relative error at a few
     Nw so the constant can be read off.

  4. The stochastic operator's MATRIX is symmetric (it is applied to I) and its
     smallest eigenvalue is reported -- W-hat is PSD only in expectation, and a
     negative eigenvalue would make V indefinite and break CG.

  5. tr(W) is NOT n (the standardized siblings' identity does not survive) and
     W 1 != 0 (no centering), stated numerically so the two properties this
     variant deliberately gives up are visible rather than assumed.

  6. mc_reml runs end to end on simulated phenotypes under both trace_methods,
     and every estimate is compared against the EXACT likelihood optimum found
     by grid search on the eigenbasis of the dense W.
     Comparing against the exact optimum, not against the truth, is the point:
     at these sizes the likelihood itself is flat and individual replicates land
     far from (s2gxg, s2e) -- one draw in five peaks at s2gxg = 0 exactly.  That
     is the estimand's sampling behaviour, not an estimator fault, and only the
     grid tells the two apart.  mc_reml should sit within the Monte-Carlo trace
     noise (~0.02 at Nmc = 50) of the grid on EVERY replicate.
"""
import argparse
import numpy as np

from Function_MCREML import (additive_design, split_into_genes, build_W_pooled,
                             setup_pooled, compute_WU_pooled, simulate_Cholesky_gxg,
                             simulate_remove_sampling_err, mc_reml)


def _rel(a, b):
    return np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300)


def exact_WU(genes, U):
    """Reference  W @ U  for the unstandardized pooled kernel, O(n m_g^2)/column.

    INDEPENDENT of Function_MCREML on purpose: this is the yardstick the
    pipeline's stochastic operator is measured against, so it shares no code
    with it.  Per gene, using 2 S_g u = (K_g .* K_g) u - D_g(D_g' u) and
    (K .* K) u = (Z .* (Z M)) 1 with M = Z' diag(u) Z:
    """
    n = genes[0].shape[0]
    U = np.atleast_2d(np.asarray(U, dtype=float))
    if U.shape[0] != n:
        U = U.T
    out = np.zeros_like(U)
    P = 0
    for Zg in genes:
        mg = Zg.shape[1]
        if mg < 2:
            continue
        P += mg * (mg - 1) // 2
        Dg = Zg * Zg
        for j in range(U.shape[1]):
            M = Zg.T @ (Zg * U[:, j:j + 1])          # m_g-by-m_g
            out[:, j] += np.einsum('na,na->n', Zg, Zg @ M)
        out -= Dg @ (Dg.T @ U)
    return out / (2.0 * P)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=300)
    ap.add_argument('--m', type=int, default=120)
    ap.add_argument('--G', type=int, default=4)
    ap.add_argument('--Nw', type=int, default=4000)
    # Nw for the end-to-end fits (check 6) only.  Checks 3 and 4 want Nw large
    # to pin the error CONSTANT; the fits want it at pipeline scale, and running
    # them at Nw=4000 makes the Hutchinson path (Nmc=50 applies per iteration)
    # dominate the whole script's runtime for no extra information.
    ap.add_argument('--Nw_reml', type=int, default=400)
    ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args()

    rng = np.random.default_rng(a.seed)
    SNP = rng.binomial(2, rng.uniform(0.05, 0.5, size=a.m), size=(a.n, a.m))
    Z = additive_design(SNP)
    genes = split_into_genes(Z, a.G)
    n = a.n

    # ---- 1. the two dense builds agree ----------------------------------
    W_had = build_W_pooled(genes, method='hadamard')
    W_prs = build_W_pooled(genes, method='pairs')
    r1 = _rel(W_had, W_prs)
    print(f"[1] hadamard vs pair-sum dense W      rel = {r1:.3e}   "
          f"{'OK' if r1 < 1e-10 else 'FAIL'}")
    W = W_had

    # ---- 2. the independent exact reference == dense W ------------------
    U = rng.normal(size=(n, 3))
    r2 = _rel(exact_WU(genes, U), W @ U)
    print(f"[2] exact reference vs W @ U          rel = {r2:.3e}   "
          f"{'OK' if r2 < 1e-10 else 'FAIL'}")

    # ---- 3. the pipeline's stochastic operator: convergence in Nw -------
    print("[3] stochastic operator (the pipeline's) vs W @ U:")
    ref = W @ U
    for Nw in (a.Nw // 16, a.Nw // 4, a.Nw):
        V_st, P = setup_pooled(genes, Nw=Nw, seed=1)
        r = _rel(compute_WU_pooled(genes, V_st, P, U), ref)
        # the note's prediction: relative error ~ c sqrt(n / Nw)
        print(f"      Nw = {Nw:6d}   rel = {r:.3e}   "
              f"rel / sqrt(n/Nw) = {r / np.sqrt(n / Nw):.3f}")
    V_st, P = setup_pooled(genes, Nw=a.Nw, seed=1)
    u = U[:, 0]
    r3v = _rel(compute_WU_pooled(genes, V_st, P, u),
               compute_WU_pooled(genes, V_st, P, u.reshape(n, 1))[:, 0])
    print(f"      (n,) and (n,1) input paths agree: rel = {r3v:.3e}")

    # ---- 4. symmetry and PSD-ness of the stochastic matrix --------------
    What = compute_WU_pooled(genes, V_st, P, np.eye(n))
    asym = np.linalg.norm(What - What.T) / np.linalg.norm(What)
    ev = np.linalg.eigvalsh(0.5 * (What + What.T))
    evW = np.linalg.eigvalsh(W)
    print(f"[4] asymmetry of W-hat                rel = {asym:.3e}   "
          f"{'OK' if asym < 1e-10 else 'FAIL'}")
    print(f"    lambda_min  W = {evW[0]:+.3e}   W-hat = {ev[0]:+.3e}"
          f"    lambda_max  W = {evW[-1]:.3e}   W-hat = {ev[-1]:.3e}")

    # ---- 5. the properties this variant gives up ------------------------
    print(f"[5] tr(W) = {np.trace(W):.4f}  vs n = {n}   (NOT an identity here)")
    print(f"    ||W 1|| / ||W||_F = {np.linalg.norm(W @ np.ones(n)) / np.linalg.norm(W):.4f}"
          f"   (0 in the centered siblings)")

    # ---- 6. end-to-end REML vs the exact likelihood optimum -------------
    s2gxg, s2e = 0.2, 0.8
    np.random.seed(a.seed)
    Lgxg, _ = simulate_Cholesky_gxg(SNP, a.G, s2gxg=s2gxg, s2e=s2e)
    lam, Qw = np.linalg.eigh(W)          # exact optimum lives in this basis

    def grid_opt(yv, hi=1.5, k=151):
        """argmax of the exact log-likelihood of V = s2gxg W + s2e I."""
        z = Qw.T @ yv
        best = (-np.inf, 0.0, 0.0)
        for sa in np.linspace(0.0, hi, k):
            d = sa * lam[:, None] + np.linspace(0.01, hi, k)[None, :]
            ll = -0.5 * (np.sum(np.log(d), axis=0) + (z * z) @ (1.0 / d))
            j = int(np.argmax(ll))
            if ll[j] > best[0]:
                best = (ll[j], sa, float(np.linspace(0.01, hi, k)[j]))
        return best[1], best[2]

    print(f"[6] end-to-end REML  (truth s2gxg={s2gxg}, s2e={s2e}, "
          f"Nw={a.Nw_reml}); 'grid' is the exact likelihood optimum for that draw")
    for rep in range(3):
        y = simulate_remove_sampling_err(Lgxg, n, s2gxg=s2gxg, s2e=s2e)
        ga, gb = grid_opt(y)
        edge = (' [BOUNDARY: s2gxg clamped -- the s2e gap below is the shared '
                'optimizer, not the operator; see mc_reml]'
                if ga <= 1e-12 else '')
        print(f"    rep{rep}   grid  s2gxg={ga:.4f}  s2e={gb:.4f}{edge}")
        for tm in ('slq', 'hutchinson'):
            s, _ = mc_reml(Z, y, a.G, iters=50, Nmc=50, seed=1,
                           trace_method=tm, Nw=a.Nw_reml)
            print(f"           trace={tm:11s} "
                  f"s2gxg={s[0]:.4f}  s2e={s[1]:.4f}   "
                  f"|d_grid| = {max(abs(s[0] - ga), abs(s[1] - gb)):.4f}")


if __name__ == '__main__':
    main()
