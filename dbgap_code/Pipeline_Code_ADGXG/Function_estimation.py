# -*- coding: utf-8 -*-
import numpy as np
from scipy.optimize import minimize
import scipy.sparse.linalg as spla
####################################################################
# Function_estimation.py -- MC AI-REML for variance components, matrix-free.
#
#     y = X b + g + e ,   V = sum_i s_i K_i + s_e I ,
#
# with the covariates X as FIXED effects.  Genetic components: additive (K_A) and
# dominance (K_D), from Function_input.genetic_kernels, and within-gene
# additive-by-additive epistasis (K_GxG), from Function_input.gxg_kernel.
#
# KERNELS ARE ONE FUNCTION  B -> [K_1 B, ..., K_K B]  ((n,) or (n, k) in and out),
# all kernels from one pass over the genotypes, so the fit never needs Z, D or
# any n x n matrix (int8 genotypes, standardized on the fly).
#
# COVARIATES BY PROJECTION.  M = I - X (X'X)^-1 X' (X includes the intercept).
# y and the Hutchinson probes are projected, every kernel is applied as M K_i M
# and the noise kernel is M.  This is exact REML with X as fixed effects: the
# fit lives in the range of M, and projected probes estimate the RESTRICTED
# traces tr(P K_i) that REML needs.
#
# THE OPTIMIZER is BOLT-REML's (Loh et al. 2015, Supplementary Note 3.2-3.5):
# Monte-Carlo AI-REML with fixed Hutchinson probes, a two-phase probe schedule,
# and each step the solution of a small QP -- the AI quadratic model maximized
# over the box s_i >= 0, s_e >= 1e-9 var(y), inside an adaptive trust region.
# A component whose likelihood peaks at 0 sits on its bound while the others
# converge; its SE is reported as NaN.
#
# The data term V^{-1} y is solved in the SAME batched CG as the probes,
# B = [y, U]: CG treats columns independently, so x is unchanged, and the
# merged solve is one wider apply instead of two.
####################################################################


# ------------------------------------------------------------------ covariates
def make_projector(X):
    """M = I - X (X'X)^-1 X' as a LinearOperator: M @ B for (n,) or (n, k) arrays.

    Uses an orthonormal basis Q of col(X) (QR), so M B = B - Q (Q' B): O(n p k).
    M is symmetric, so rmatvec = matvec.
    """
    X = np.asarray(X, dtype=np.float64)
    Q, _ = np.linalg.qr(X)
    f = lambda B: B - Q @ (Q.T @ B)
    n = X.shape[0]
    return spla.LinearOperator((n, n), matvec=f, matmat=f, rmatvec=f, rmatmat=f, dtype=np.float64)

# ------------------------------------------------------------------ linear algebra

def _cg_batched(V, Bmat, x0=None, tol=1e-6, maxiter=1000):
    """Conjugate gradient for the SPD system V X = Bmat.

    V      : (n, n) LinearOperator (or anything with V @ X), applied to (n, c) blocks.
    Bmat   : (n, c) right-hand sides, solved together with per-column scalars,
             so one V-pass advances every column.
    x0     : (n, c) warm start; None (or all zeros) starts from X = 0.

    Returns (X, ok).  ok is False, and X is NOT a solution, when
      - a column meets a search direction with p'Vp <= 0: V is not positive
        definite here and CG is not a valid solver, or
      - maxiter iterations pass without every column reaching tol.
    Costs (iterations taken) products V @ (n, c) block, plus one more for a
    nonzero x0 (the initial residual).
    """
    n, c = Bmat.shape
    X = np.zeros((n, c)) if x0 is None else x0.copy()
    R = Bmat.copy() if x0 is None or not np.any(x0) else Bmat - V @ X
    P = R.copy()
    rs_old = np.sum(R * R, axis=0)
    b_norm = np.sqrt(np.sum(Bmat * Bmat, axis=0))
    b_norm[b_norm == 0.0] = 1.0

    for _ in range(maxiter):
        VP = V @ P
        pVp = np.sum(P * VP, axis=0)
        if np.any((pVp <= 0.0) & (rs_old > 0.0)):
            return X, False
        alpha = rs_old / pVp
        X += alpha * P
        R -= alpha * VP
        rs_new = np.sum(R * R, axis=0)
        if np.max(np.sqrt(rs_new) / b_norm) < tol:
            break
        beta = rs_new / rs_old
        P = R + beta * P
        rs_old = rs_new
    else:
        return X, False                      
    return X, True

# ------------------------------------------------------------------ standard errors
def reml_se(AI, Nmc, at_bound=None):
    """SEs of the AI-REML estimate, with the Monte-Carlo inflation.

        SE_p = sqrt( [AI_free^{-1}]_pp ) * sqrt(1 + 1/Nmc)

    The second factor is what MC AI-REML adds (BOLT-REML note 2.3): the
    Hutchinson traces match the data's quadratics to an average over Nmc
    simulated references rather than to their expectations, inflating the
    variance by (1 + 1/Nmc) -- 0.5% on the SE at Nmc = 100.  Uses the FINE Nmc.

    AI is restricted to the components NOT at their bound.  A component at its
    bound gets NaN: analytic SEs are not valid on the boundary (BOLT-REML).
    """
    k = AI.shape[0]
    free = (np.ones(k, dtype=bool) if at_bound is None
            else ~np.asarray(at_bound, dtype=bool))
    se = np.full(k, np.nan)
    if free.any():
        sub = AI[np.ix_(free, free)]
        se[free] = np.sqrt(np.diag(np.linalg.inv(sub))) * np.sqrt(1.0 + 1.0 / Nmc)
    return se


def h2_se(s, AI, Nmc, at_bound=None):
    """h2_i = s_i / sum(s) for every genetic component, and its SE by the delta
    method on AI^{-1} (with the same (1 + 1/Nmc) inflation as reml_se).

    Returns (h2, se_h2), arrays of length len(s) - 1.

    As in reml_se, components at their bound are held FIXED there: AI is
    inverted on the free components only, and the bound ones get zero
    variance in the delta method.  A component at its bound gets an SE of NaN;
    the free components keep theirs.
    """
    s = np.asarray(s, dtype=float)
    k = len(s)
    total = s.sum()
    h2 = s[:-1] / total
    se = np.full(h2.shape, np.nan)
    free = (np.ones(k, dtype=bool) if at_bound is None
            else ~np.asarray(at_bound, dtype=bool))
    if free.any():
        cov = np.zeros((k, k))                            # bound components: fixed, no variance
        cov[np.ix_(free, free)] = (np.linalg.inv(AI[np.ix_(free, free)])
                                   * (1.0 + 1.0 / Nmc))
        for i in np.flatnonzero(free[:-1]):               # free genetic components only
            grad = -s[i] / total ** 2 * np.ones(k)        # d(s_i / total) / d s_j
            grad[i] += 1.0 / total
            se[i] = np.sqrt(grad @ cov @ grad)
    return h2, se


# ------------------------------------------------------------------ trust-region step
def _trust_region_step(g, A, lo, Delta, jitter=1e-10):
    """The BOLT-REML step subproblem (Supplementary Note 3.2.1, 3.4):

        maximize  g'p - 0.5 p'A p
        s.t.      p >= lo ,   ||diag(A) p|| <= Delta ,

    lo = (lower bounds) - theta, so lo <= 0 at a feasible theta.  The norm
    constraint is dropped while Delta is Inf.  A only gets a jitter of
    jitter * mean|diag A| so the QP is well posed.

    Returns (p, gain, how): gain is the objective at p -- dLL_pred, >= 0 by
    construction since p = 0 is feasible -- and how names the route taken:

      'newton'   the Newton step A^-1 g is feasible, hence the exact solution;
      'slsqp'    SLSQP (one of BOLT's three NLopt solvers), started from the
                 Newton step clipped into the box, then POLISHED on its active
                 set (exact to round-off when it satisfies the KKT conditions);
      'fallback' SLSQP failed: the Newton step clipped into the box and
                 shortened to the radius;
      'cauchy'   that fallback predicts a loss: the projected-gradient step,
                 whose gain is never negative.
    """
    k = g.shape[0]
    dvec = np.diag(A).copy()                 # the radius scaling, diag(A)
    Aj = A + jitter * (np.abs(dvec).mean() + 1e-300) * np.eye(k)
    finite = np.isfinite(Delta)

    def gain(p):
        return float(g @ p - 0.5 * p @ (Aj @ p))

    def rnorm(p):
        return float(np.linalg.norm(dvec * p))

    def clip_to_domain(p):
        """Clip into the box, then shorten onto the radius (stays in the box,
        since lo <= 0 and shrinking moves p toward 0)."""
        p = np.maximum(p, lo)
        if finite and rnorm(p) > Delta:
            p = p * (Delta / rnorm(p))
        return p

    p_newton = np.linalg.solve(Aj, g)
    if np.all(p_newton >= lo) and (not finite or rnorm(p_newton) <= Delta):
        return p_newton, gain(p_newton), 'newton'

    p0 = np.maximum(p_newton, lo)
    cons = []
    if finite:
        cons = [{'type': 'ineq',
                 'fun': lambda p: Delta ** 2 - np.sum((dvec * p) ** 2),
                 'jac': lambda p: -2.0 * dvec ** 2 * p}]
    try:
        res = minimize(lambda p: -gain(p), p0, jac=lambda p: -(g - Aj @ p),
                       method='SLSQP', bounds=[(l, None) for l in lo],
                       constraints=cons,
                       options={'ftol': 1e-12, 'maxiter': 200})
        ok = bool(res.success) and np.all(np.isfinite(res.x))
    except (ValueError, np.linalg.LinAlgError):
        ok = False
    if ok:
        p = clip_to_domain(res.x)               # remove round-off violations
        # Active-set polish: fix the components SLSQP left on their bound and
        # solve the rest exactly; keep it only if it is a KKT point and within
        # the radius.
        scale = np.abs(lo).max() + np.abs(p).max() + 1e-300
        act = (p - lo) <= 1e-7 * scale
        q = np.where(act, lo, 0.0)
        free = ~act
        if free.any():
            rhs = g[free] - Aj[np.ix_(free, act)] @ lo[act]
            q[free] = np.linalg.solve(Aj[np.ix_(free, free)], rhs)
        grad_q = g - Aj @ q                     # gradient of the objective at q
        gscale = np.abs(g).max() + 1e-300
        if (np.all(q[free] >= lo[free])
                and np.all(grad_q[act] <= 1e-10 * gscale)
                and (not finite or rnorm(q) <= Delta)
                and gain(q) >= gain(p) - 1e-12 * max(1.0, abs(gain(p)))):
            p = q
        if gain(p) >= 0.0:
            return p, gain(p), 'slsqp'

    p = clip_to_domain(p_newton)
    if gain(p) >= 0.0:
        return p, gain(p), 'fallback'

    # Projected gradient: drop the components pinned at their bound whose
    # gradient points out of the box, then the best point along it.
    d = g.copy()
    d[(lo >= 0.0) & (d < 0.0)] = 0.0
    if not np.any(d):
        return np.zeros(k), 0.0, 'cauchy'    # KKT point: nothing to gain
    curv = float(d @ (Aj @ d))
    t = float(g @ d) / curv if curv > 0.0 else np.inf
    neg = d < 0.0
    if neg.any():
        t = min(t, float(np.min(lo[neg] / d[neg])))
    if finite:
        t = min(t, Delta / rnorm(d))
    p = t * d
    return p, gain(p), 'cauchy'


# ------------------------------------------------------------------ MC AI-REML
def mc_reml(kernels, y, X, iters=30, Nmc=100, Nmc_coarse=15, cg_tol=1e-6,
            cg_maxiter=1000, jitter=1e-10, tol_ll=1e-4, tol_ll_coarse=1e-2,
            eta1=1e-4, eta2=0.99, alpha1=0.25, alpha2=3.5,
            s_init=None, seed=None, verbose=False):
    """Monte-Carlo AI-REML for  V = sum_i s_i K_i + s_e I  with fixed effects X,
    optimized as in BOLT-REML (Loh et al. 2015, Supplementary Note 3.2-3.5).

    kernels : ONE function B -> [K_1 B, ..., K_K B], (n,) or (n, k) in and out,
              computing all kernels in one call (one pass over the genotypes),
              e.g. Function_input.genetic_kernels(G, ('A', 'D')).  A list of
              functions B -> K_i B is also accepted.
    y       : (n,) phenotype.
    X       : (n, p) covariates, including the intercept column.

    Covariates (exact REML): y and the probes are projected by M = I - X(X'X)^-1 X',
    every kernel is applied as M K_i M, and the noise kernel is M.

    Parameter domain (Note 3.2.1): s_i >= 0, s_e >= 1e-9 var(y); no upper bound.
    s_init gives the starting point as FRACTIONS of var(y) (default: 0.1 split
    over the genetic components, 0.9 for noise).

    Score traces tr(P K_i) by HUTCHINSON with Nmc fixed Rademacher probes
    (projected), K_i U formed once.  TWO-PHASE SCHEDULE (Note 3.5): Nmc_coarse
    probes to tol_ll_coarse, then all Nmc to tol_ll; the probes are drawn ONCE
    and the coarse phase uses the first Nmc_coarse of them.

    Each iteration solves the QP  maximize g'p - 0.5 p'AI p  over the box, inside
    an adaptive trust region (_trust_region_step).  The actual gain is approximated
    by the TRAPEZOID rule 0.5 p'(g(s) + g(s+p)); rho = actual / predicted drives
    accept / reject and the radius.  Converged when dLL_pred < tol_ll in the fine
    phase.

    Returns
    -------
    s    : (K + 1,) estimates (s_1, ..., s_K, s_e).
    AI   : (K + 1, K + 1) average information at the FULL Nmc.  SEs via
           reml_se(AI, Nmc, info['at_bound']).
    info : dict -- converged, n_iters, n_reject, switch_it, at_bound.
    """
    M = make_projector(X)
    y = (M @ np.asarray(y, dtype=float).reshape(-1, 1))[:, 0]
    n = y.shape[0]
    if callable(kernels):
        kern_all = kernels
    else:
        kern_list = list(kernels)
        kern_all = lambda B: [f(B) for f in kern_list]

    def kall(B):
        """[M K_1 M B, ..., M K_K M B], all kernels from ONE call of kern_all."""
        return [M @ o for o in kern_all(M @ B)]

    rng = np.random.default_rng(seed)
    # Drawn ONCE at FULL width; the coarse phase reads U[:, :Nmc_coarse].
    U = M @ rng.choice([-1.0, 1.0], size=(n, Nmc))
    KU = kall(U)                                    # K_i U for the fixed probes, once
    K = len(KU)                                     # number of genetic components
    k = K + 1

    vary = float(y @ y) / (n - X.shape[1])          # residual variance of y given X
    lower = np.zeros(k)
    lower[-1] = 1e-9 * vary                         # the domain, Note 3.2.1
    if s_init is None:
        s_init = np.r_[np.full(K, 0.1 / K), 0.9]
    s = vary * np.asarray(s_init, dtype=float)
    if s.shape != (k,):
        raise ValueError(f"s_init must have {k} entries, got {s.shape}")
    if np.any(s < lower):
        raise ValueError(f"s_init={tuple(s_init)} violates the bounds")

    if Nmc_coarse is None or Nmc_coarse >= Nmc:
        n_probe, coarse = Nmc, False                # single-phase fit
    else:
        n_probe, coarse = Nmc_coarse, True
    switch_it = None

    yc = y.reshape(n, 1)
    xbuf = np.zeros((n, 1))                         # warm starts for the CG solves
    Pbuf = np.zeros((n, Nmc))
    Gbuf = None

    def vmat(s_at, B):
        """V B = s_e M B + sum_i s_i M K_i M B, all kernels from ONE call."""
        MB = M @ B
        outs = kern_all(MB)
        comb = s_at[0] * outs[0]
        for i in range(1, K):
            comb = comb + s_at[i] * outs[i]
        return s_at[-1] * MB + M @ comb

    def eval_grad_ai(s_at):
        """Score and average information at s_at, at the CURRENT n_probe.
        Returns (score, AI, ok, sols); ok is False when CG met negative curvature.
        sols are the CG solutions at s_at.  They become the warm starts only via
        commit(sols), i.e. once s_at is accepted: a rejected trial point (e.g. a
        near-singular V) must not leave its solutions as later starting points."""
        f = lambda B: vmat(s_at, B)
        V = spla.LinearOperator((n, n), matvec=f, matmat=f, rmatvec=f, rmatmat=f, dtype=np.float64)

        # x = V^{-1} y and the probe solves V^{-1} U: ONE batched, warm-started CG
        Uc_ = U[:, :n_probe]
        Sol, ok = _cg_batched(V, np.hstack([yc, Uc_]),
                              x0=np.hstack([xbuf, Pbuf[:, :n_probe]]),
                              tol=cg_tol, maxiter=cg_maxiter)
        if not ok:
            return None, None, False, None
        X_y, Psol = Sol[:, :1], Sol[:, 1:]
        x_ = X_y[:, 0]

        # score_i = 0.5 (x' K_i x - tr(P K_i)), traces by Hutchinson
        Kx = kall(x_) + [M @ x_]
        KUs = KU + [Uc_]
        sc = np.array([0.5 * (x_ @ Kx[i]
                              - np.mean(np.sum(Psol * KUs[i][:, :n_probe], axis=0)))
                       for i in range(k)])

        # average information: AI_ij = 0.5 (K_i x)' V^{-1} (K_j x)
        KX = np.column_stack(Kx)
        G_sol, ok = _cg_batched(V, KX, x0=Gbuf, tol=cg_tol, maxiter=cg_maxiter)
        if not ok:
            return None, None, False, None
        ai = 0.5 * (KX.T @ G_sol)
        return sc, 0.5 * (ai + ai.T), True, (X_y, Psol, G_sol)

    def commit(sols):
        """Make the CG solutions at the current point the warm starts."""
        nonlocal xbuf, Gbuf
        xbuf, Psol, Gbuf = sols
        Pbuf[:, :Psol.shape[1]] = Psol

    # Delta = Inf: the radius only exists once a step has been rejected
    Delta = np.inf
    score, AI, ok, sols = eval_grad_ai(s)
    if not ok:
        raise RuntimeError(f"V is not positive definite at the starting point s = {s}")
    commit(sols)
    n_reject = 0
    n_iters = 0
    converged = False

    for it in range(iters):
        n_iters = it + 1
        step, dLL_pred, how = _trust_region_step(score, AI, lower - s, Delta,
                                                 jitter=jitter)
        if verbose:
            print(f"iter {it:2d}  s/var(y)={np.round(s / vary, 5)}  dLL_pred={dLL_pred:.4e}"
                  f"  step={how}  phase={'coarse' if coarse else 'fine'}(S={n_probe})",
                  flush=True)
        # p = 0 is feasible, so the QP optimum is >= 0 by construction
        assert dLL_pred >= -1e-8 * max(1.0, abs(dLL_pred)), (
            f"dLL_pred={dLL_pred:.6e} < 0 at iter {it} ({how} step)")

        # Coarse phase: passing the loose tolerance switches to the full Nmc
        if coarse:
            if dLL_pred < tol_ll_coarse:
                coarse = False
                n_probe = Nmc
                switch_it = it
                score_f, AI_f, ok, sols = eval_grad_ai(s)    # re-evaluate at full Nmc
                if not ok:
                    break
                score, AI = score_f, AI_f
                commit(sols)
                continue
        elif dLL_pred < tol_ll:
            converged = True
            break

        # Trial point, trapezoid gain, accept / reject
        s_try = np.maximum(s + step, lower)
        taken = s_try - s
        score_try, AI_try, ok, sols_try = eval_grad_ai(s_try)
        if not ok:
            rho = -1.0                                 # V not PD at the trial point
        elif np.linalg.norm(score_try) > 2.0 * np.linalg.norm(score):
            rho = -1.0                                 # quadratic model broke down
        else:
            dLL_approx = float(taken @ (0.5 * (score + score_try)))
            rho = dLL_approx / dLL_pred if dLL_pred != 0.0 else -1.0

        taken_dnorm = np.linalg.norm(np.diag(AI) * taken)
        if rho > eta1:
            s, score, AI = s_try, score_try, AI_try     # accept
            commit(sols_try)                           # only accepted points seed later solves
        else:
            n_reject += 1
        if rho < eta1:
            Delta = alpha1 * taken_dnorm               # shrink onto the bad step
        elif rho >= eta2:
            Delta = max(Delta, alpha2 * taken_dnorm)   # model is good: expand

    at_bound = [bool(v) for v in (s - lower <= 1e-8 * vary)]
    info = {'converged': converged, 'n_iters': n_iters, 'n_reject': n_reject,
            'switch_it': switch_it, 'at_bound': at_bound}
    if verbose:
        print(f"converged={converged}; iters={n_iters}; rejected={n_reject}; "
              f"switch at iter {switch_it}; at_bound={at_bound}", flush=True)
    return s, AI, info