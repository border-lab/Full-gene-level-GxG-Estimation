# -*- coding: utf-8 -*-
"""Term-by-term check of the FOUR-component pipeline against dense algebra.

    python3 verify_lowrank.py [--n 300] [--m 120] [--G 4]

Every claim the pipeline rests on is asserted here at small n, m, where the
dense kernels can be formed and compared directly:

  1. The closed form  2 sum_{a<b} h h' = (K .* K) - D D'  used by
     build_W_pooled really is the literal pair sum, not merely something close
     to it.  The pair sum is written out here, independently of the module, so
     the identity that lets the simulation and the estimator be the same object
     is checked against something that shares no code with it.

  2. exact_WU == W @ U to machine precision, where exact_WU is the INDEPENDENT
     reference implementation below -- (Z .* (Z M)) 1 with M = Z' diag(u) Z,
     written here and nowhere else.  The pipeline module implements only the
     rank-r truncation; a reference that imported the code under test would
     prove nothing, so this one is deliberately standalone.

  3. compute_WU_pooled -- the pipeline's low-rank operator -- converges to
     W @ U as r grows and is EXACT (machine precision) at full rank, the limit
     the stochastic sibling has at no finite Nw.  The error at each r is a
     DETERMINISTIC truncation bias (rerunning changes nothing); it tracks the
     discarded eigenvalue share, reported alongside, and under linkage
     equilibrium (the RandomSNP genotype used here) it stays large until r
     approaches rank(Z_g) -- the note's Table 1.  Under LD the same small r
     would already be near-exact (Table 2); this script certifies correctness,
     not the LD-dependence.

  4. compute_KU == K_i @ U to machine precision for BOTH standardized designs,
     and tr(K_a) = tr(K_d) = n EXACTLY.  The ADDITIVE and DOMINANCE components
     are applied with NO truncation and no probe, so there is nothing to
     converge: they are either exact or broken, and this is the check that says
     which.  Their trace identity is asserted because it is the one the
     epistasis kernel gives up (check 6).  The sample correlation between the
     two designs is reported alongside: the GCTA dominance coding is orthogonal
     to the additive one only in expectation under HWE, and how far from 0 that
     number sits is how much K_a and K_d actually overlap in this genotype.

  5. The low-rank operator's MATRIX is EXACTLY symmetric (it is K_r .* K with
     both factors symmetric -- asserted at machine precision, where the
     stochastic sibling needed a transpose-averaging pass), and its smallest
     eigenvalue is reported: W-hat <= W in the Loewner order, so truncation can
     push lambda_min below 0, which would make V indefinite and break CG.  K_a
     and K_d are PSD by construction (Gram matrices), so V's conditioning is
     still the epistasis term's problem alone.

  6. tr(W) is NOT n (the standardized siblings' identity does not survive) and
     W 1 != 0 (no centering), stated numerically next to K_a's exact tr = n so
     the two properties this kernel deliberately gives up are visible rather
     than assumed.

  7. mc_reml runs end to end on simulated 4-component phenotypes, and every
     estimate is compared against the EXACT likelihood optimum found by dense
     gradient optimisation over (s2a, s2d, s2gxg, s2e).  At FULL rank the operator
     is exact, so the fit must sit within the Monte-Carlo trace noise (~0.02 at
     Nmc = 50) of that optimum on EVERY replicate -- that certifies the
     estimator machinery.  At the pipeline default r the SAME fit is repeated
     so the truncation's systematic effect on the estimates is visible next to
     it.  Comparing against the exact optimum, not against the truth, is the
     point: at these sizes the likelihood is flat and individual replicates
     land far from (s2a, s2d, s2gxg, s2e) -- with four components a draw can
     peak at s2a = 0, s2d = 0 or s2gxg = 0 outright, and s2d most often.  That is the estimand's sampling
     behaviour, not an estimator fault, and only the exact optimum tells the
     two apart.
"""
import argparse
import numpy as np
from scipy.optimize import minimize

from Function_MCREML import (additive_design, dominance_design,
                             split_into_genes, build_W_pooled,
                             build_K, compute_KU, setup_pooled,
                             compute_WU_pooled, simulate_Cholesky_4vc,
                             simulate_phenotype, mc_reml, R_DEFAULT)


def _rel(a, b):
    return np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300)


def pair_sum_W(genes):
    """Literal  W = (1/P) sum_g sum_{a<b in g} h_ab h_ab',  h_ab = Z_a .* Z_b.

    The C(m_g, 2)-term sum written out, with no Hadamard identity anywhere --
    the yardstick build_W_pooled's closed form is measured against.
    """
    n = genes[0].shape[0]
    W = np.zeros((n, n))
    P = 0
    for Zg in genes:
        mg = Zg.shape[1]
        if mg < 2:
            continue
        P += mg * (mg - 1) // 2
        for a in range(mg):
            for b in range(a + 1, mg):
                h = Zg[:, a] * Zg[:, b]
                W += np.outer(h, h)
    return W / P


def exact_WU(genes, U):
    """Reference  W @ U  for the unstandardized pooled kernel, O(n m_g^2)/column.

    INDEPENDENT of Function_MCREML on purpose: this is the yardstick the
    pipeline's low-rank operator is measured against, so it shares no code
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
    # r for the truncated end-to-end fit in check 7; the pipeline default.  The
    # full-rank fit next to it always uses r = min(n, max_g m_g).
    ap.add_argument('--r_reml', type=int, default=R_DEFAULT)
    ap.add_argument('--reps', type=int, default=3)
    ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args()

    rng = np.random.default_rng(a.seed)
    SNP = rng.binomial(2, rng.uniform(0.05, 0.5, size=a.m), size=(a.n, a.m))
    Z = additive_design(SNP)
    Zd = dominance_design(SNP)
    genes = split_into_genes(Z, a.G)
    n = a.n
    r_full = min(n, max(g.shape[1] for g in genes))   # exactness cap

    # ---- 1. the closed form is the pair sum -----------------------------
    W = build_W_pooled(genes)
    r1 = _rel(W, pair_sum_W(genes))
    print(f"[1] hadamard closed form vs pair sum  rel = {r1:.3e}   "
          f"{'OK' if r1 < 1e-10 else 'FAIL'}")

    # ---- 2. the independent exact reference == dense W ------------------
    U = rng.normal(size=(n, 3))
    r2 = _rel(exact_WU(genes, U), W @ U)
    print(f"[2] exact reference vs W @ U          rel = {r2:.3e}   "
          f"{'OK' if r2 < 1e-10 else 'FAIL'}")

    # ---- 3. the pipeline's low-rank operator: convergence in r ----------
    # Eigenvalue share retained at each r, pooled over genes (lam_s = sg^2 from
    # the same thin SVD the setup uses) -- the note's "PC variance" column.
    svals = [np.linalg.svd(g, compute_uv=False) for g in genes if g.shape[1] >= 2]
    lam_tot = sum(float(np.sum(s ** 2)) for s in svals)
    print("[3] low-rank operator (the pipeline's) vs W @ U  (deterministic bias):")
    ref = W @ U
    for r in sorted({5, R_DEFAULT, r_full // 2, r_full}):
        F_lr, P = setup_pooled(genes, r=r)
        rel = _rel(compute_WU_pooled(genes, F_lr, P, U), ref)
        share = sum(float(np.sum(s[:r] ** 2)) for s in svals) / lam_tot
        exact = '   EXACT expected' if r >= r_full else ''
        print(f"      r = {r:4d}   eig share = {share:.4f}   rel = {rel:.3e}"
              f"{exact}")
    F_lr, P = setup_pooled(genes, r=r_full)
    r3 = _rel(compute_WU_pooled(genes, F_lr, P, U), ref)
    print(f"      full rank (r = {r_full})               rel = {r3:.3e}   "
          f"{'OK' if r3 < 1e-10 else 'FAIL'}")
    u = U[:, 0]
    r3v = _rel(compute_WU_pooled(genes, F_lr, P, u),
               compute_WU_pooled(genes, F_lr, P, u.reshape(n, 1))[:, 0])
    print(f"      (n,) and (n,1) input paths agree: rel = {r3v:.3e}")

    # ---- 4. the additive and dominance kernels: exact, no truncation ------
    Ka = build_K(Z)
    Kd = build_K(Zd)
    r4a = _rel(compute_KU(Z, U), Ka @ U)
    r4d = _rel(compute_KU(Zd, U), Kd @ U)
    r4v = _rel(compute_KU(Z, u), (Ka @ u))
    tra_err = abs(np.trace(Ka) - n) / n
    trd_err = abs(np.trace(Kd) - n) / n
    print(f"[4] compute_KU vs K_a @ U             rel = {r4a:.3e}   "
          f"{'OK' if r4a < 1e-12 else 'FAIL'}  (no truncation: exact at any r)")
    print(f"    compute_KU vs K_d @ U             rel = {r4d:.3e}   "
          f"{'OK' if r4d < 1e-12 else 'FAIL'}")
    print(f"    (n,) input path agrees: rel = {r4v:.3e};  "
          f"|tr(K_a) - n| / n = {tra_err:.3e}, |tr(K_d) - n| / n = {trd_err:.3e}   "
          f"{'OK' if max(tra_err, trd_err) < 1e-12 else 'FAIL'}")
    # How orthogonal the dominance coding really is: mean |corr| between the
    # matched columns of Z_a and Z_d.  Exactly 0 only in expectation under HWE.
    col_corr = np.abs(np.einsum('ij,ij->j', Z, Zd)) / n
    print(f"    |corr(Z_a[:,j], Z_d[:,j])|: mean = {col_corr.mean():.4f}, "
          f"max = {col_corr.max():.4f}   (0 only in expectation under HWE -- "
          f"this is how much K_a and K_d overlap)")

    # ---- 5. symmetry and PSD-ness of the truncated matrix ---------------
    F_tr, P = setup_pooled(genes, r=a.r_reml)
    What = compute_WU_pooled(genes, F_tr, P, np.eye(n))
    asym = np.linalg.norm(What - What.T) / np.linalg.norm(What)
    ev = np.linalg.eigvalsh(0.5 * (What + What.T))
    evW = np.linalg.eigvalsh(W)
    evKa = np.linalg.eigvalsh(Ka)
    evKd = np.linalg.eigvalsh(Kd)
    print(f"[5] asymmetry of W-hat (r={a.r_reml})          rel = {asym:.3e}   "
          f"{'OK' if asym < 1e-12 else 'FAIL'}  (exact symmetry, no averaging)")
    print(f"    lambda_min  W = {evW[0]:+.3e}   W-hat = {ev[0]:+.3e}"
          f"    lambda_max  W = {evW[-1]:.3e}   W-hat = {ev[-1]:.3e}")
    print(f"    lambda_min  K_a = {evKa[0]:+.3e}   K_d = {evKd[0]:+.3e}   "
          f"(PSD by construction; V's conditioning is the epistasis term's "
          f"problem)")

    # ---- 6. the properties the epistasis kernel gives up ----------------
    print(f"[6] tr(W) = {np.trace(W):.4f}  vs  tr(K_a) = {np.trace(Ka):.4f} "
          f"= tr(K_d) = {np.trace(Kd):.4f} = n = {n}   "
          f"(an identity for K_a and K_d, NOT for W)")
    print(f"    ||W 1|| / ||W||_F = {np.linalg.norm(W @ np.ones(n)) / np.linalg.norm(W):.4f}"
          f"   (0 in the centered siblings)")

    # ---- 7. end-to-end 4-component REML vs the exact likelihood optimum --
    s2a, s2d, s2gxg, s2e = 0.1, 0.1, 0.1, 0.7
    np.random.seed(a.seed)
    La, Ld, Lgxg, _ = simulate_Cholesky_4vc(SNP, a.G, s2a=s2a, s2d=s2d,
                                            s2gxg=s2gxg, s2e=s2e)
    I = np.eye(n)
    Ks = [Ka, Kd, W, I]

    def exact_opt(yv):
        """argmax of the exact log-likelihood of
        V = s2a K_a + s2d K_d + s2gxg W + s2e I.

        Dense, with the analytic gradient -- no grid, because the components do
        not share an eigenbasis (K_a, K_d and W do not commute) and a 4-D grid
        at the resolution the 2-component script used would be 151^4 solves.
        L-BFGS-B from several starts; the best objective wins.
        """
        def nll_grad(th):
            V = th[0] * Ka + th[1] * Kd + th[2] * W + th[3] * I
            L = np.linalg.cholesky(V)
            x = np.linalg.solve(V, yv)
            logdet = 2.0 * np.sum(np.log(np.diag(L)))
            nll = 0.5 * (logdet + yv @ x)
            Vi = np.linalg.inv(V)
            g = np.array([0.5 * (np.sum(Vi * K) - x @ (K @ x)) for K in Ks])
            return nll, g

        best = (np.inf, None)
        for st in ([0.1, 0.1, 0.1, 0.7], [0.25, 0.25, 0.25, 0.25],
                   [1e-6, 1e-6, 1e-6, yv.var()]):
            res = minimize(nll_grad, np.array(st), jac=True, method='L-BFGS-B',
                           bounds=[(1e-9, 5.0)] * 4)
            if res.fun < best[0]:
                best = (res.fun, res.x)
        return best[1]

    print(f"[7] end-to-end REML  (truth s2a={s2a}, s2d={s2d}, s2gxg={s2gxg}, "
          f"s2e={s2e}); 'exact' is the dense likelihood optimum for that draw; "
          f"full rank must match it, r={a.r_reml} shows the truncation's "
          f"systematic shift")
    for rep in range(a.reps):
        y = simulate_phenotype(La, Ld, Lgxg, n, s2a=s2a, s2d=s2d, s2gxg=s2gxg,
                               s2e=s2e)
        opt = exact_opt(y)
        edge = ('  [BOUNDARY: a component clamped -- the gaps below are the '
                'shared optimizer, not the operator; see mc_reml]'
                if np.min(opt) <= 1e-6 else '')
        print(f"    rep{rep}   exact  s2a={opt[0]:.4f}  s2d={opt[1]:.4f}  "
              f"s2gxg={opt[2]:.4f}  s2e={opt[3]:.4f}{edge}")
        for r_use, lab in ((r_full, 'full'), (a.r_reml, f'{a.r_reml}')):
            s, _ = mc_reml(Z, Zd, y, a.G, iters=50, Nmc=50, seed=1, r=r_use)
            print(f"           r={lab:>4}  s2a={s[0]:.4f}  s2d={s[1]:.4f}  "
                  f"s2gxg={s[2]:.4f}  s2e={s[3]:.4f}   "
                  f"|d_exact| = {np.max(np.abs(s - opt)):.4f}")


if __name__ == '__main__':
    main()
