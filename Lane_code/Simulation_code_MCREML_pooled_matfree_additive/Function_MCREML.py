# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.linalg import cholesky
from scipy.optimize import minimize
import time

####################################################################
# TWO variance components: additive + noise.
#
#     y = g_a + e ,
#     V = s2a K_a + s2e I ,
#     K_a = Z_a Z_a'/m      (standardized design).
#
# This is the ADDITIVE-ONLY member of the Lowrank_Wu_std_4VC family: the same
# optimizer, probes, CG, counters and timers, with the dominance kernel K_d and
# the pooled within-gene epistasis kernel W removed.  So there are no genes, no
# G, no truncation level r, no c-hat and no stored-A apply; K_a is applied
# EXACTLY, and any error in s2a-hat is sampling error, never operator error.
# The design is column-standardized, so c_a = 1 exactly and s2a is on the
# realized-variance scale with no correction.
#
# EVERYTHING RUNS IN THE CONTRAST SPACE  P_c = I - 11'/n.  The phenotype, the
# component draws and the Hutchinson probes are centered.  Centered probes
# estimate the RESTRICTED trace tr(P_c V^-1 K_i), which supplies the -1/s2e
# correction on the residual component: proper REML with an intercept.  K_a
# needs no wrapping: its design is column-standardized, so K_a 1 = 0 already.
#
# THE OPTIMIZER is BOLT-REML's (Loh et al. 2015, Supplementary Note 3.2-3.5):
# Monte-Carlo AI-REML with fixed Hutchinson probes, a two-phase probe schedule,
# and each step the solution of a small QP -- the AI quadratic model maximized
# over the box s2a >= 0, s2e >= 1e-9 var(y), inside an adaptive trust region.
# A component whose likelihood peaks at 0 sits on its bound while the other
# converges; its SE is reported as NaN.  V is PD on the whole box here (K_a is
# exact and PSD, s2e > 0), so CG's p'Vp > 0 check (cg_negcurv) should never
# fire; it is kept as a guard.
#
# The data term V^{-1} y is solved in the SAME batched CG as the probes,
# B = [y, U]: CG treats columns independently, so x is unchanged, and the
# merged solve is one wider apply instead of two.
####################################################################


# --------------------------------------------- operator-apply counting
# The estimator's whole cost is APPLIES of K_a and I to blocks of vectors,
# almost all inside CG.  Unlike wall-clock, these counts are
# machine-independent.  Per operator: <op>_applies (calls) and <op>_columns
# (total n-vectors, the cost-bearing number).  V is the composite, K is
# compute_KU.  Also cg_solves, cg_iters, reml_iters and cg_negcurv (0 on a
# clean replicate).  GLOBAL and CUMULATIVE; mc_reml zeroes them, so one
# MC_REML call = one replicate.
#
# WALL-CLOCK TWINS.  _OP_TIMES holds seconds, keyed "<what>_sec", accumulated
# by ONE perf_counter pair per call at the same sites (V, K applies), plus the
# solve groups and setup pieces inside mc_reml.  These ARE machine-dependent:
# counts say how much work was asked for, these say how fast each kind of work
# ran.  Nesting: V_sec CONTAINS the K_sec spent inside V applies; every
# solve_*_sec contains its V applies.  So the pieces are quoted per column
# (driver), not summed.
_OP_COUNTS = {}
_OP_TIMES = {}


def reset_op_counts():
    """Zero every operator-apply counter AND timer.  mc_reml calls this on entry."""
    _OP_COUNTS.clear()
    _OP_TIMES.clear()


def get_op_counts():
    """Snapshot of counters and timers as ONE plain dict.  Absent keys mean zero.

    Timer keys end in "_sec"; everything else is a count.
    """
    out = dict(_OP_COUNTS)
    out.update(_OP_TIMES)
    return out


def _add_time(name, dt):
    """Accumulate dt seconds under timer `name` (the "_sec" suffix is added)."""
    key = name + '_sec'
    _OP_TIMES[key] = _OP_TIMES.get(key, 0.0) + dt


def _count_apply(name, ncols):
    """Record one apply of operator `name` on `ncols` n-vectors."""
    a, c = name + '_applies', name + '_columns'
    _OP_COUNTS[a] = _OP_COUNTS.get(a, 0) + 1
    _OP_COUNTS[c] = _OP_COUNTS.get(c, 0) + int(ncols)


def _count_event(name, k=1):
    """Record k occurrences of a plain counted event (no column dimension)."""
    _OP_COUNTS[name] = _OP_COUNTS.get(name, 0) + k


# ------------------------------------------------------------------ design
def _standardize_cols(M, stability_std=1e-12):
    """Column-standardize to mean 0, variance 1 (ddof=0).

    A near-constant column (std < stability_std) is left un-scaled.
    """
    M = np.asarray(M, dtype=float)
    mu = M.mean(axis=0)
    sd = M.std(axis=0)
    sd = np.where(sd < stability_std, 1.0, sd)
    return (M - mu) / sd


def _center_cols(B):
    """Project onto the contrast space: subtract the column mean.

    P_c = I - 11'/n applied column-wise.  K_a already annihilates 1 (its
    design is column-standardized), so the model lives entirely in this space
    and the 1-direction carries no information.
    """
    return B - B.mean(axis=0, keepdims=True)


def additive_design(real_data):
    """Additive design Z_a: column-standardized allele dosages.

    Feeds K_a = Z_a Z_a'/m, on both the simulation and the estimation side.
    """
    return _standardize_cols(real_data)


# --------------------------------------- standardized GRM (without GRM)
def compute_KU(Zi, U):
    """ K_i @ U  with  K_i = Z_i Z_i' / m.

    Zi is a standardized design; U is (n,) or (n, c) and the result matches.
    O(n m c) time, O(n c) storage -- K_i is never formed.  tr(K_i) = n exactly.
    """
    Zi = np.asarray(Zi, dtype=float)
    U = np.asarray(U, dtype=float)
    m = Zi.shape[1]
    single = (U.ndim == 1)
    if single:
        U = U.reshape(-1, 1)
    _count_apply('K', U.shape[1])
    t0 = time.perf_counter()
    out = (Zi @ (Zi.T @ U)) / m
    _add_time('K', time.perf_counter() - t0)
    return out[:, 0] if single else out


def build_K(Zi):
    """Explicit GRM  K_i = Z_i Z_i' / m.  SIMULATION only (n-by-n)."""
    Zi = np.asarray(Zi, dtype=float)
    return (Zi @ Zi.T) / Zi.shape[1]


# ------------------------------------------------------------- simulation
def simulate_Cholesky_additive(real_data, s2a=0.1, s2e=0.9, stability=1e-10):
    """Cholesky factor of the additive genetic covariance.

        La La' = s2a K_a .

    Returns (La, k_build_time); k_build_time covers the GRM gemm and its
    Cholesky, the O(n^2 m) + O(n^3) work that is this job.  The dense K_a is
    freed once its factor exists.
    """
    t_start = time.perf_counter()
    Za = additive_design(real_data)
    n = Za.shape[0]
    Ka = build_K(Za)
    La = cholesky(s2a * Ka + stability * np.eye(n), lower=True)
    del Ka
    k_build_time = time.perf_counter() - t_start
    return La, k_build_time


def _force_var(v, target, name):
    """Rescale v so its ddof=0 sample variance is exactly target.

    target = 0 means the component is absent; return zeros.  A nonzero target
    with a degenerate draw is a real fault and raises.
    """
    vv = float(np.var(v, ddof=0))
    if target <= 0.0:
        return np.zeros_like(v)
    if vv <= 0.0:
        raise ValueError(f"cannot force Var-hat({name}) to {target}: the draw "
                         f"has zero sample variance (is its Cholesky factor "
                         f"identically zero?)")
    return v * np.sqrt(target / vv)


def simulate_phenotype(La, n, s2a=0.1, s2e=0.9,
                       return_realized=False, force_realized=True):
    """Draw one phenotype  y = g_a + e,

        g_a = P_c La u1 ,   e = P_c sqrt(s2e) u2 ,

    the two draws independent and CENTERED.  Returns y, or (y, V_a, V_e) with
    return_realized=True; both Var-hat at ddof=0.

    Because La La' = (s2a/m) Z_a Z_a', g_a has exactly the law of Z_a beta with
    beta ~ N(0, (s2a/m) I_m), and E[Var-hat(g_a)] = s2a (c_a = 1 exactly).

    force_realized=True (THE DEFAULT, as in the 4VC parent) rescales each
    component by sqrt(target / Var-hat) so its realized variance hits the
    target exactly.  An estimate's deviation from the target is then the
    ESTIMATOR's error alone, but IT CHANGES THE GENERATING LAW: dividing by a
    random sample variance leaves y non-Gaussian, so REML fits a model the
    simulation does not obey (a small misspecification at these n).  A
    component with target 0 is left at zero.  force_realized=False restores
    the plain draw.

    NOTE the two do NOT sum to var(y) under either setting: the sample
    cross-product is not zero, and rescaling touches only the diagonal terms.
    """
    u1 = np.random.randn(n)
    u2 = np.random.randn(n)

    # CENTERED draws; the estimator works in the same contrast space.
    a = _center_cols(La @ u1)             # additive effect   ~ N(0, s2a K_a)
    e = _center_cols(np.sqrt(s2e) * u2)   # residual noise

    if force_realized:
        # ddof=0, matching how the realized columns report it.  See the
        # docstring: this changes the generating law.
        a = _force_var(a, s2a, 'a')
        e = _force_var(e, s2e, 'e')

    y = a + e
    if return_realized:
        return y, float(a.var()), float(e.var())
    return y


# ------------------------------------------------ Linear operator for the 2-VC system
def _v_matvec(Z, s2a, s2e, B):
    """Apply  V = s2a K_a + s2e I  to B, matrix-free.

    B is (n,) or (n, k) and the result matches.  No n-by-n GRM exists: K_a is
    one gemm pair against Z.  THE UNIT OF WORK: CG sees V and nothing else; a
    centered right-hand side stays centered, since both terms preserve it.
    """
    _count_apply('V', 1 if np.ndim(B) == 1 else np.shape(B)[1])
    t0 = time.perf_counter()
    out = s2a * compute_KU(Z, B) + s2e * B
    _add_time('V', time.perf_counter() - t0)    # includes its K apply
    return out


def _cg_batched(matvec, Bmat, x0=None, tol=1e-6, maxiter=1000):
    """Conjugate gradient for the SPD system V X = Bmat.

    matvec : X -> V X, (n, c) in and out.
    Bmat   : (n, c) right-hand sides, solved together with per-column scalars,
             so one V-pass advances every column.
    x0     : (n, c) warm start.

    Returns (X, ok).  ok is False when a column meets a search direction with
    p'Vp <= 0: V is not positive definite here, CG is not a valid solver, and
    X is NOT a solution.  Counted as cg_negcurv.  With V = s2a K_a + s2e I on
    the box this cannot happen short of round-off; the check is a guard.
    Costs 1 + (iterations taken) applies of matvec.
    """
    n, c = Bmat.shape
    _count_event('cg_solves')
    X = np.zeros((n, c)) if x0 is None else x0.copy()
    R = Bmat - matvec(X)
    P = R.copy()
    rs_old = np.sum(R * R, axis=0)
    b_norm = np.sqrt(np.sum(Bmat * Bmat, axis=0))
    b_norm[b_norm == 0.0] = 1.0

    for _ in range(maxiter):
        _count_event('cg_iters')
        VP = matvec(P)
        pVp = np.sum(P * VP, axis=0)
        if np.any((pVp <= 0.0) & (rs_old > 0.0)):
            _count_event('cg_negcurv')
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
    return X, True


def reml_se(AI, Nmc, at_bound=None):
    """SEs of the AI-REML estimate, with the Monte-Carlo inflation.

        SE_p = sqrt( [AI_free^{-1}]_pp ) * sqrt(1 + 1/Nmc)

    The second factor is what MC AI-REML adds (BOLT-REML note 2.3): the
    Hutchinson traces match the data's quadratics to an average over Nmc
    simulated references rather than to their expectations, inflating the
    variance by (1 + 1/Nmc) -- 0.5% on the SE at Nmc = 100.  Uses the FINE Nmc.

    AI is restricted to the components NOT at their bound.  A component at its
    bound gets NaN: analytic SEs are not valid on the boundary (BOLT-REML).

    The ONLY place SEs should come from.
    """
    k = AI.shape[0]
    free = (np.ones(k, dtype=bool) if at_bound is None
            else ~np.asarray(at_bound, dtype=bool))
    se = np.full(k, np.nan)
    if free.any():
        sub = AI[np.ix_(free, free)]
        se[free] = np.sqrt(np.diag(np.linalg.inv(sub))) * np.sqrt(1.0 + 1.0 / Nmc)
    return se


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
                 Newton step clipped into the box.  SLSQP stops at its
                 tolerance (~1e-8 in p), so its active set -- the components
                 it put on their bound -- is then POLISHED: the free
                 components are solved exactly with the active ones fixed,
                 and the result is used when it satisfies the KKT conditions
                 and the radius.  That makes the step exact to round-off,
                 and so a deterministic function of (g, A);
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
        # solve the rest exactly; keep it only if it is a KKT point (free
        # components inside the box, no active component that the gradient
        # would lift off its bound) and within the radius.
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


# ----------------------------------------------------------- MC AI-REML (k=2)
def mc_reml(Z, y, iters=30, Nmc=100, Nmc_coarse=15, cg_tol=1e-6,
            cg_maxiter=1000, jitter=1e-10, tol_ll=1e-4, tol_ll_coarse=1e-2,
            eta1=1e-4, eta2=0.99, alpha1=0.25, alpha2=3.5,
            s_init=(0.05, 0.95), seed=None, verbose=False):
    """Monte-Carlo AI-REML for  V = s2a K_a + s2e I, optimized as in BOLT-REML
    (Loh et al. 2015, Supplementary Note 3.2-3.5).

    The fit runs in the CONTRAST SPACE: y and the probes are centered, which
    makes this proper REML with an intercept.

    Z is the standardized genotype (K_a, all m SNPs); it is the only large
    array held.  K_a is applied EXACTLY, so any error in s2a-hat is sampling
    error, never operator error.

    Parameter domain (Note 3.2.1)
    -----------------------------
    s2a >= 0 and s2e >= 1e-9 var(y); no upper bound.  s_init is a fraction of
    var(y) and must satisfy the bounds.

    Score traces
    ------------
    tr(V^{-1} K_i) by HUTCHINSON: Nmc fixed centered Rademacher probes, K_a U
    formed once.  The data term x = V^{-1} y is solved in the SAME
    warm-started batched CG as the probes, on [y, U]: CG treats columns
    independently, so x is the same solution.  The AI solve stays separate,
    since its right-hand side needs x.  Fixed probes make the objective a
    deterministic function of s.

    TWO-PHASE SCHEDULE (Note 3.5): Nmc_coarse probes to tol_ll_coarse, then
    all Nmc to tol_ll.  The probes are drawn ONCE at full width and the coarse
    phase uses the first Nmc_coarse columns, so the switch EXTENDS the objective
    rather than replacing it -- the products are formed once, the fixed point
    moves continuously, and Delta carries over.  The step that triggers the
    switch is discarded and re-evaluated at full Nmc.  Nmc_coarse=None (or
    >= Nmc) gives the single-phase fit.

    Constrained trust-region step (Note 3.2.1, 3.4)
    -----------------------------------------------
    Each iteration solves the 2-dimensional QP

        maximize  g'p - 0.5 p'AI p
        s.t.      s + p >= lower bounds ,  ||diag(AI) p|| <= Delta

    (_trust_region_step: SLSQP, the norm constraint dropped while Delta is
    Inf).  A component whose likelihood peaks at 0 lands ON its bound.
    dLL_pred, the QP objective at the solution, is the predicted gain; the
    actual gain would need log|V|, which this construction never forms, so it
    is approximated by the TRAPEZOID rule 0.5 p'(g(s) + g(s+p)) -- the
    likelihood being the line integral of the score.  rho = actual/predicted
    drives the accept test and the radius update; a gradient that more than
    doubles means the quadratic model has broken down and sets rho = -1.
    Delta starts at Inf.  An accepted step's trial score is reused, so only
    REJECTED steps cost an extra score evaluation.  Converged when
    dLL_pred < tol_ll in the fine phase.

    Returns
    -------
    s    : (2,) estimated (s2a, s2e).
    AI   : (2, 2) average information at the FULL Nmc.  SEs via
           reml_se(AI, Nmc, info['at_bound']); a raw sqrt(diag(inv(AI)))
           understates them.
    info : dict -- converged (dLL_pred < tol_ll reached in the fine phase),
           n_iters, n_reject, switch_it, at_bound (bool per component).
    """
    reset_op_counts()                  # counts describe THIS replicate alone
    y = np.asarray(y, dtype=float).flatten()
    y = y - y.mean()                   # work in the contrast space throughout
    Z = np.asarray(Z, dtype=float)
    n = y.shape[0]
    k = 2

    # --- setup: genotype-only, done once ---
    # Only the fixed-probe product K_a U (setup_KU); setup_sec is the total.
    t_setup0 = time.perf_counter()
    Kaapply = lambda B: compute_KU(Z, B)

    vary = y.var()
    lower = np.array([0.0, 1e-9 * vary])   # the domain, Note 3.2.1

    rng = np.random.default_rng(seed)
    # Drawn ONCE at FULL width; the coarse phase reads U[:, :Nmc_coarse].  See
    # the docstring for why this is not a redraw at the switch.
    # CENTERED Rademacher probes: E[uu'] = P_c, so Hutchinson estimates
    # tr(P_c V^-1 K_i) -- the RESTRICTED trace REML wants with an intercept.
    U = _center_cols(rng.choice([-1.0, 1.0], size=(n, Nmc)))
    if Nmc_coarse is None or Nmc_coarse >= Nmc:
        n_probe, coarse = Nmc, False                # single-phase fit
    else:
        n_probe, coarse = Nmc_coarse, True
    switch_it = None                                # iteration the switch fired

    # --- starting point, as FRACTIONS OF var(y) -----------------------------
    # vary * s_init keeps the estimator scale-equivariant.
    s = vary * np.asarray(s_init, dtype=float)
    if s.shape != (k,):
        raise ValueError(f"s_init must have {k} entries, got {s.shape}")
    if np.any(s < lower):
        raise ValueError(f"s_init={tuple(s_init)} violates the bounds "
                         f"(genetic >= 0, s2e >= 1e-9 var(y)).")
    AI = np.eye(k)

    yc = y.reshape(n, 1)
    xbuf = np.zeros((n, 1))   # warm-start buffers for the CG solve groups
    # Held at the FULL Nmc so the fine phase inherits the coarse phase's warm
    # starts; the untouched columns are cold-started once, at the switch.
    Pbuf = np.zeros((n, Nmc))
    Gbuf = None

    t0 = time.perf_counter()
    KaU = Kaapply(U)  # K_a U for the fixed probes: genotype-only, formed once
    _add_time('setup_KU', time.perf_counter() - t0)
    _add_time('setup', time.perf_counter() - t_setup0)
    t_phase0 = time.perf_counter()      # coarse phase starts here

    def eval_grad_ai(s_at):
        """Score and average information at s_at, at the CURRENT n_probe.

        Returns (score, AI, ok).  ok is False when a CG solve met negative
        curvature (V not PD at s_at); score and AI are then None and the
        warm-start buffers are left untouched.  Lifted out of the loop because
        the trust region needs the score at a trial point too, and both must be
        the same function of s or rho compares two different objectives.
        """
        nonlocal xbuf, Gbuf
        s2a_, s2e_ = s_at
        mv = lambda B: _v_matvec(Z, s2a_, s2e_, B)

        # --- x = V^{-1} y and the probe solves V^{-1} U, ONE batched CG ------
        # on [y, U[:, :n_probe]] (n_probe + 1 columns), warm-started from the
        # previous solution of each.  CG treats the columns independently.
        Uc_ = U[:, :n_probe]
        t0 = time.perf_counter()
        Sol, ok = _cg_batched(mv, np.hstack([yc, Uc_]),
                              x0=np.hstack([xbuf, Pbuf[:, :n_probe]]),
                              tol=cg_tol, maxiter=cg_maxiter)
        _add_time('solve_probe_coarse' if coarse else 'solve_probe_fine',
                  time.perf_counter() - t0)
        if not ok:
            return None, None, False
        X_y, Psol = Sol[:, :1], Sol[:, 1:]
        x_ = X_y[:, 0]

        # data quadratic x'K_a x
        Kax_ = Kaapply(x_)

        # --- score traces tr(V^{-1} K_i): Hutchinson ------------------------
        sc = np.array([
            0.5 * (x_ @ Kax_ - np.mean(np.sum(Psol * KaU[:, :n_probe], axis=0))),
            0.5 * (x_ @ x_ - np.mean(np.sum(Psol * Uc_, axis=0)))])

        # --- average information: A_ij = 0.5 (K_i x)' V^{-1}(K_j x) ---
        KX = np.column_stack([Kax_, x_])
        t0 = time.perf_counter()
        G_sol, ok = _cg_batched(mv, KX, x0=Gbuf, tol=cg_tol, maxiter=cg_maxiter)
        _add_time('solve_ai', time.perf_counter() - t0)
        if not ok:
            return None, None, False
        xbuf, Gbuf = X_y, G_sol
        Pbuf[:, :n_probe] = Psol
        ai = 0.5 * (KX.T @ G_sol)
        return sc, 0.5 * (ai + ai.T), True

    # Delta = Inf: the radius only exists once a step has been rejected, so a
    # replicate that never overshoots follows the plain AI-Newton path.
    Delta = np.inf
    score, AI, ok = eval_grad_ai(s)      # the current gradient, carried forward
    if not ok:
        raise RuntimeError(f"V is not positive definite at the starting point "
                           f"s = {s}; raise s_init's s2e share.")
    n_reject = 0
    n_iters = 0
    converged = False

    for it in range(iters):
        _count_event('reml_iters')
        n_iters = it + 1

        # --- the constrained QP step inside the adaptive trust region ------
        step, dLL_pred, how = _trust_region_step(score, AI, lower - s, Delta,
                                                 jitter=jitter)

        # --- CONVERGENCE: BOLT-REML predicted log-likelihood gain ------------
        # From l(s + d) ~= l(s) + score'd - 0.5 d' AI d, dLL_pred is the model's
        # remaining height above the iterate within the feasible set, so
        # dLL_pred < tol_ll means within ~sqrt(tol_ll) SE of the constrained
        # optimum.  Tested BEFORE the step is taken, following BOLT.
        if verbose:
            print(f"iter {it:2d}  dLL_pred={dLL_pred:.6e}  step={how}"
                  f"  phase={'coarse' if coarse else 'fine'}(S={n_probe})"
                  f"  Delta={Delta:.4g}", flush=True)
        # p = 0 is feasible, so the QP optimum is >= 0 by construction; checked
        # because failing silently yields a plausible-looking number.
        assert dLL_pred >= -1e-8 * max(1.0, abs(dLL_pred)), (
            f"dLL_pred={dLL_pred:.6e} < 0 at iter {it} ({how} step); "
            f"eigs(AI)={np.linalg.eigvalsh(AI)}")

        # In the COARSE phase, passing the loose tolerance switches the schedule
        # rather than ending the fit.  The triggering step is discarded; Delta
        # carries over, the radius being independent of the probe count.
        if coarse:
            if dLL_pred < tol_ll_coarse:
                # Phase stamp: everything since setup was the coarse phase.
                t_now = time.perf_counter()
                _add_time('phase_coarse', t_now - t_phase0)
                t_phase0 = t_now
                coarse = False
                n_probe = Nmc
                switch_it = it
                if verbose:
                    cnt = get_op_counts()
                    print(f"iter {it:2d}  SWITCH coarse -> fine "
                          f"(S={Nmc_coarse} -> {Nmc}): cg_iters so far="
                          f"{cnt.get('cg_iters', 0)}, V_columns so far="
                          f"{cnt.get('V_columns', 0)}", flush=True)
                score_f, AI_f, ok = eval_grad_ai(s)   # re-evaluate at full Nmc
                if not ok:
                    # The extra probes found negative curvature at the accepted
                    # point itself: nowhere to step back to, so stop unconverged.
                    if verbose:
                        print(f"iter {it:2d}  V not PD at the full Nmc; "
                              f"stopping unconverged", flush=True)
                    break
                score, AI = score_f, AI_f
                continue
        elif dLL_pred < tol_ll:
            converged = True
            break

        # --- trial point, trapezoid gain, accept / reject -------------------
        # The QP keeps s + step inside the box; the maximum only removes
        # round-off below a bound.  rho uses the actual move s_try - s.
        s_try = np.maximum(s + step, lower)
        taken = s_try - s
        score_try, AI_try, ok = eval_grad_ai(s_try)

        if not ok:
            rho = -1.0                    # V not PD at the trial point: reject
        else:
            # Model breakdown: a gradient that GREW by more than 2x means the
            # quadratic model does not describe this region.  Reject outright.
            gnorm, gnorm_try = np.linalg.norm(score), np.linalg.norm(score_try)
            if gnorm_try > 2.0 * gnorm:
                rho = -1.0
            else:
                # TRAPEZOID rule: l(s+p) - l(s) is the line integral of the
                # score, the one computable thing (the likelihood needs log|V|).
                dLL_approx = float(taken @ (0.5 * (score + score_try)))
                rho = dLL_approx / dLL_pred if dLL_pred != 0.0 else -1.0

        taken_dnorm = np.linalg.norm(np.diag(AI) * taken)
        if rho > eta1:
            # Accept; reusing the trial score keeps the cost at one score/iter.
            s, score, AI = s_try, score_try, AI_try
        else:
            n_reject += 1
        if rho < eta1:
            Delta = alpha1 * taken_dnorm          # shrink onto the bad step
        elif rho >= eta2:
            Delta = max(Delta, alpha2 * taken_dnorm)   # model is good: expand

        if verbose:
            print(f"iter {it:2d}  s={s}  rho={rho:+.4f}"
                  f"  {'ACCEPT' if rho > eta1 else 'REJECT'}"
                  f"{'' if ok else ' (V not PD)'}"
                  f"  |D p|={taken_dnorm:.4g}  Delta->{Delta:.4g}", flush=True)

    # Phase stamp: from the switch (or from setup, single-phase) to here.
    _add_time('phase_fine', time.perf_counter() - t_phase0)
    at_bound = [bool(v) for v in (s - lower <= 1e-8 * vary)]
    info = {'converged': converged, 'n_iters': n_iters, 'n_reject': n_reject,
            'switch_it': switch_it, 'at_bound': at_bound}
    if verbose:
        print(f"switch at iter {switch_it}; rejected {n_reject} step(s); "
              f"converged={converged}; at_bound={at_bound}; "
              f"final Delta={Delta:.4g}; SE={reml_se(AI, Nmc, at_bound)}",
              flush=True)
    return s, AI, info


def MC_REML(Z, y, iters=30, Nmc=100, Nmc_coarse=15, cg_tol=1e-6,
            cg_maxiter=1000, seed=None, verbose=False):
    """Wrapper: returns (s2a_hat, s2e_hat, AI, info).

    SEs from reml_se(AI, Nmc, info['at_bound']).  Nmc_coarse=None gives the
    single-phase fit.  info is mc_reml's: converged, n_iters, n_reject,
    switch_it, at_bound.
    """
    s, AI, info = mc_reml(Z, y, iters=iters, Nmc=Nmc,
                          Nmc_coarse=Nmc_coarse, cg_tol=cg_tol,
                          cg_maxiter=cg_maxiter, seed=seed, verbose=verbose)
    s2a_hat, s2e_hat = s
    return s2a_hat, s2e_hat, AI, info
