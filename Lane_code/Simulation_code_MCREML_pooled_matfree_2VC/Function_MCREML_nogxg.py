# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.linalg import cholesky
import time

####################################################################
# TWO genetic variance components (hence "_2VC") + noise -- THREE parameters:
#
#     y = g_a + g_d + e ,
#     V = s2a K_a + s2d K_d + s2e I ,
#     K_a = Z_a Z_a'/m ,  K_d = Z_d Z_d'/m      (standardized designs).
#
# This is the _Lowrank_Wu_std_4VC pipeline with the pooled within-gene
# epistasis component REMOVED: no W, no genes, no c-hat, no rank-r truncation,
# no G.  What is left is the part of that pipeline that was already EXACT --
# K_a and K_d are applied matrix-free as one gemm pair each, and the estimator
# fits the same kernels the phenotype was drawn from, so any error in an
# estimate is sampling or optimizer error, never operator error.  It exists
# to give the 4VC runs a baseline: the same optimizer, the same probes, the
# same CG, the same timers, on the system without the expensive kernel.
#
# EVERYTHING RUNS IN THE CONTRAST SPACE  P_c = I - 11'/n.  The phenotype, the
# component draws and the Hutchinson probes are centered, so this is proper
# REML with an intercept (centered probes estimate the RESTRICTED trace
# tr(P_c V^-1 K_i), which supplies the -1/s2e correction on the residual
# component).  K_a and K_d need no wrapping: their designs are
# column-standardized, so K 1 = 0 already.  tr(K_a) = tr(K_d) = n exactly.
####################################################################


# --------------------------------------------- operator-apply counting
# The estimator's whole cost is APPLIES of K_a, K_d and I to blocks of vectors,
# almost all inside CG.  Unlike wall-clock, these counts are machine-
# independent.  Per operator: <op>_applies (calls) and <op>_columns (total
# n-vectors, the cost-bearing number).  V is the composite, K is compute_KU
# serving both designs.  Also cg_solves, cg_iters, reml_iters and
# lam_min_exact (0 on a clean replicate).  GLOBAL and CUMULATIVE; mc_reml
# zeroes them, so one MC_REML call = one replicate.
#
# WALL-CLOCK TWINS.  _OP_TIMES holds seconds, keyed "<what>_sec", accumulated
# by ONE perf_counter pair per call at the same sites (V and K applies), plus
# the solve groups and setup pieces inside mc_reml.  These ARE machine-
# dependent -- that is their point: counts say how much work was asked for,
# these say how fast each kind of work ran.  Nesting: V_sec CONTAINS the K_sec
# spent inside V applies; every solve_*_sec contains its V applies.  So the
# pieces are quoted per column (driver), not summed.
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


def _count_apply(name, ncols):
    """Record one apply of operator `name` on `ncols` n-vectors."""
    a, c = name + '_applies', name + '_columns'
    _OP_COUNTS[a] = _OP_COUNTS.get(a, 0) + 1
    _OP_COUNTS[c] = _OP_COUNTS.get(c, 0) + int(ncols)


def _count_event(name, k=1):
    """Record k occurrences of a plain counted event (no column dimension)."""
    _OP_COUNTS[name] = _OP_COUNTS.get(name, 0) + k


def _add_time(name, dt):
    """Accumulate dt seconds under timer `name` (the "_sec" suffix is added)."""
    key = name + '_sec'
    _OP_TIMES[key] = _OP_TIMES.get(key, 0.0) + dt


# ------------------------------------------------------------------ designs
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

    P_c = I - 11'/n applied column-wise.  Both kernels here already annihilate
    1 (their designs are column-standardized), so the model lives entirely in
    this space and the 1-direction carries no information.
    """
    return B - B.mean(axis=0, keepdims=True)


def additive_design(real_data):
    """Additive design Z_a: column-standardized allele dosages.  Feeds K_a."""
    return _standardize_cols(real_data)


def dominance_design(real_data, maf_floor=1e-6):
    """Dominance design Z_d: column-standardized GCTA dominance coding.

    With p = column mean / 2 and q = 1 - p:  0 -> -p/q, 1 -> 1, 2 -> -q/p, which
    has mean 0 under HWE -- the orthogonality to Z_a that makes s2a and s2d
    separately estimable, though only in expectation and only under HWE.
    Dosages are rounded first (the coding is defined on hard calls) and p is
    clipped to [maf_floor, 1 - maf_floor] against a monomorphic column.
    """
    X = np.rint(np.asarray(real_data, dtype=float))
    p = X.mean(axis=0) / 2.0
    p = np.clip(p, maf_floor, 1.0 - maf_floor)
    q = 1.0 - p
    D = (X == 0) * (-p / q) + (X == 1) * 1.0 + (X == 2) * (-q / p)
    return _standardize_cols(D)


# --------------------------------------- standardized GRMs (without GRM)
# One routine serves K_a and K_d: same structure, different design.
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
    _count_apply('K', U.shape[1])            # K_a and K_d share this routine
    t0 = time.perf_counter()
    out = (Zi @ (Zi.T @ U)) / m
    _add_time('K', time.perf_counter() - t0)
    return out[:, 0] if single else out


def build_K(Zi):
    """Explicit GRM  K_i = Z_i Z_i' / m.  SIMULATION only (n-by-n)."""
    Zi = np.asarray(Zi, dtype=float)
    return (Zi @ Zi.T) / Zi.shape[1]


# ------------------------------------------------------------- simulation
def simulate_Cholesky_2vc(real_data, s2a=0.1, s2d=0.1, s2e=0.7,
                          stability=1e-10):
    """Cholesky factors of the two genetic covariances.

        La La' = s2a K_a ,  Ld Ld' = s2d K_d .

    Returns (La, Ld, k_build_time); k_build_time covers both GRM builds and
    both factorizations -- SIMULATION-ONLY cost, the estimator never pays it.

    Built SEQUENTIALLY, each dense GRM freed once its factor exists, so the
    peak is ~3 n-by-n arrays.  Without the epistasis kernel this is O(n^2 m)
    for the two gemms and O(n^3) for the two Choleskys, and it is still the
    memory-critical job of the pipeline.
    """
    Za = additive_design(real_data)
    Zd = dominance_design(real_data)
    n, m = Za.shape

    t_start = time.perf_counter()

    # --- additive --------------------------------------------------------
    Ka = build_K(Za)
    La = cholesky(s2a * Ka + stability * np.eye(n), lower=True)
    del Ka

    # --- dominance -------------------------------------------------------
    Kd = build_K(Zd)
    Ld = cholesky(s2d * Kd + stability * np.eye(n), lower=True)
    del Kd

    k_build_time = time.perf_counter() - t_start
    return La, Ld, k_build_time


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


def simulate_phenotype(La, Ld, n, s2a=0.1, s2d=0.1, s2e=0.7,
                       return_realized=False, force_realized=True):
    """Draw one phenotype  y = g_a + g_d + e,

        g_a = P_c La u1 ,   g_d = P_c Ld u2 ,   e = P_c sqrt(s2e) u3 ,

    the three draws independent and CENTERED.  Returns y, or
    (y, V_a, V_d, V_e) with return_realized=True; all Var-hat at ddof=0.

    g_a = Z_a beta with beta ~ N(0, (s2a/m) I_m) in law, and likewise g_d.
    Both designs are column-standardized, so E[Var-hat(g_a)] = s2a and
    E[Var-hat(g_d)] = s2d with no realized-variance factor to correct.

    force_realized=True (THE DEFAULT, as in the 4VC parent) rescales each
    component by sqrt(target / Var-hat) so its realized variance hits the
    target exactly.  An estimate's deviation from the target is then the
    ESTIMATOR's error alone -- at the price that y no longer has the Gaussian
    law REML fits (dividing by a random sample variance): a MISSPECIFICATION,
    small at these n.  A component with target 0 is left at zero.
    force_realized=False restores the plain draw.

    NOTE the three do NOT sum to var(y) under either setting: the sample
    cross-products are not zero, and rescaling touches only the diagonal terms.
    """
    u1 = np.random.randn(n)
    u2 = np.random.randn(n)
    u3 = np.random.randn(n)

    a = _center_cols(La @ u1)             # additive effect   ~ N(0, s2a K_a)
    d = _center_cols(Ld @ u2)             # dominance effect  ~ N(0, s2d K_d)
    e = _center_cols(np.sqrt(s2e) * u3)   # residual noise

    if force_realized:
        a = _force_var(a, s2a, 'a')
        d = _force_var(d, s2d, 'd')
        e = _force_var(e, s2e, 'e')

    y = a + d + e
    if return_realized:
        return y, float(a.var()), float(d.var()), float(e.var())
    return y


# ------------------------------------------------ Linear operator for the 3-parameter system
def _v_matvec(Z, Zd, s2a, s2d, s2e, B):
    """Apply  V = s2a K_a + s2d K_d + s2e I  to B, matrix-free.

    B is (n,) or (n, k) and the result matches.  No n-by-n GRM exists: K_a and
    K_d are a gemm pair each.  THE UNIT OF WORK: CG sees V and nothing else; a
    centered right-hand side stays centered, since every term preserves it.
    """
    _count_apply('V', 1 if np.ndim(B) == 1 else np.shape(B)[1])
    t0 = time.perf_counter()
    out = (s2a * compute_KU(Z, B)
           + s2d * compute_KU(Zd, B)
           + s2e * B)
    _add_time('V', time.perf_counter() - t0)    # includes its two K applies
    return out


def _cg_batched(matvec, Bmat, x0=None, tol=1e-6, maxiter=1000):
    """Conjugate gradient for the SPD system V X = Bmat.

    matvec : X -> V X, (n, c) in and out.
    Bmat   : (n, c) right-hand sides, solved together with per-column scalars,
             so one V-pass advances every column.
    x0     : (n, c) warm start.

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
        alpha = rs_old / np.sum(P * VP, axis=0)
        X += alpha * P
        R -= alpha * VP
        rs_new = np.sum(R * R, axis=0)
        if np.max(np.sqrt(rs_new) / b_norm) < tol:
            break
        beta = rs_new / rs_old
        P = R + beta * P
        rs_old = rs_new
    return X


def _spectral_range(apply, n, iters=80, tol=1e-7, seed=0):
    """(lam_min, lam_max) of a SYMMETRIC operator given only its apply.

    Two shifted power iterations on ONE column: the dominant-MAGNITUDE
    eigenvalue mu1 as a SIGNED Rayleigh quotient, then the dominant of
    A - mu1 I, which is the opposite extreme.  Deterministic (fixed seed).

    Both kernels here are PSD, so only the upper end is read from it in setup;
    it is kept general because feasible() also runs it on the composite V,
    which CAN be indefinite when a genetic component goes negative.
    """
    rng = np.random.default_rng(seed)

    def _dominant(ap):
        v = rng.standard_normal((n, 1))
        v /= np.linalg.norm(v)
        lam = 0.0
        for _ in range(iters):
            w = ap(v)
            nw = np.linalg.norm(w)
            if nw <= 0.0:
                return 0.0                       # the zero operator
            v = w / nw
            lam_new = float(v.T @ ap(v))
            if abs(lam_new - lam) <= tol * max(1.0, abs(lam_new)):
                return lam_new
            lam = lam_new
        return lam

    mu1 = _dominant(apply)
    mu2 = _dominant(lambda B: apply(B) - mu1 * B) + mu1
    return (min(mu1, mu2), max(mu1, mu2))


def reml_se(AI, Nmc):
    """SEs of the AI-REML estimate, with the Monte-Carlo inflation.

        SE_p = sqrt( [AI^{-1}]_pp ) * sqrt(1 + 1/Nmc)

    The second factor is what MC AI-REML adds (BOLT-REML note 2.3): the
    Hutchinson traces match the data's quadratics to an average over Nmc
    simulated references rather than to their expectations, inflating the
    variance by (1 + 1/Nmc) -- 0.5% on the SE at Nmc = 100.  Uses the FINE Nmc.

    The ONLY place SEs should come from.
    """
    return np.sqrt(np.diag(np.linalg.inv(AI))) * np.sqrt(1.0 + 1.0 / Nmc)


# ----------------------------------------------------------- MC AI-REML (k=3)
def mc_reml(Z, Zd, y, iters=30, Nmc=100, Nmc_coarse=15, cg_tol=1e-6,
            cg_maxiter=1000, jitter=1e-8, tol_ll=1e-4, tol_ll_coarse=1e-2,
            lm=1e-3, eta1=1e-4, eta2=0.99, alpha1=0.25, alpha2=3.5,
            upper_mult=1.5, s_init=(0.05, 0.05, 0.9),
            seed=None, verbose=False):
    """Monte-Carlo AI-REML for  V = s2a K_a + s2d K_d + s2e I.

    The fit runs in the CONTRAST SPACE: y and the probes are centered, which
    makes this proper REML with an intercept.

    Z is the standardized genotype (K_a); Zd is the dominance design (K_d).
    Both are n-by-m and are the only large arrays held.  Both kernels are
    EXACT, so any error in an estimate is sampling or optimizer error.

    Score traces
    ------------
    tr(V^{-1} K_i) by HUTCHINSON: Nmc fixed Rademacher probes, K_a U / K_d U
    formed once, Nmc warm-started CG solves per iteration.  Fixed probes make
    the objective a deterministic function of s.

    TWO-PHASE SCHEDULE (BOLT-REML note 3.5): Nmc_coarse probes to
    tol_ll_coarse, then all Nmc to tol_ll.  The probes are drawn ONCE at full
    width and the coarse phase uses the first Nmc_coarse columns, so the switch
    EXTENDS the objective rather than replacing it.  Nmc_coarse=None gives the
    single-phase fit.

    OPTIMIZER.  The block below is the 4VC parent's, with k = 3 instead of 4
    and the epistasis kernel removed from V and from the feasibility bound --
    nothing else changed, so a run here differs from a 4VC run in the KERNELS
    alone.  It carries the parent's KNOWN DEFECT, inherited deliberately: if a
    replicate's likelihood peaks at a component = 0 (s2e has a floor; the two
    genetic components have none, and are held in the PD cone by feasible()),
    that component sticks and the others converge to the WRONG values -- the
    Newton step is never re-projected onto the free subspace.  calc_stats.py
    counts the affected replicates.

    Returns
    -------
    s  : (3,) estimated (s2a, s2d, s2e).
    AI : (3, 3) average information at the FULL Nmc.  SEs via reml_se(AI, Nmc);
         a raw sqrt(diag(inv(AI))) understates them.
    """
    reset_op_counts()                  # counts describe THIS replicate alone
    y = np.asarray(y, dtype=float).flatten()
    y = y - y.mean()                   # work in the contrast space throughout
    Z = np.asarray(Z, dtype=float)
    Zd = np.asarray(Zd, dtype=float)
    n = y.shape[0]
    k = 3

    # --- setup: genotype-only, done once ---
    # Timed in two pieces (spectral ranges, the fixed-probe products K_i U)
    # plus their total, setup_sec.  There is no SVD piece here: nothing is
    # truncated.
    t_setup0 = time.perf_counter()
    Kaapply = lambda B: compute_KU(Z, B)
    Kdapply = lambda B: compute_KU(Zd, B)

    # Spectral ranges for the feasibility bound, genotype-only, computed ONCE.
    # K_a and K_d are Gram matrices, so their lower end is written as 0 rather
    # than estimated.
    t0 = time.perf_counter()
    _, lmax_a = _spectral_range(Kaapply, n)
    _, lmax_d = _spectral_range(Kdapply, n)
    _add_time('setup_spectral', time.perf_counter() - t0)
    lmin_k = np.array([0.0, 0.0])               # K_a, K_d PSD exactly
    lmax_k = np.array([lmax_a, lmax_d])

    def lam_min_V_bound(s_vec):
        """Worst-case lower bound on lam_min(V).  Positive => V is PD.

        s_i lam_min(K_i) when s_i > 0, s_i lam_max(K_i) when s_i < 0.  With
        both kernels PSD this is s2e + sum over NEGATIVE genetic components of
        s_i lam_max(K_i): a positive component can never take V indefinite,
        a negative one can.
        """
        g = np.asarray(s_vec[:2], dtype=float)
        worst = np.where(g > 0.0, g * lmin_k, g * lmax_k)
        return float(s_vec[2] + worst.sum())

    def feasible(s_vec):
        """Is V at s_vec positive definite?  Cheap bound first, exact if it fails.

        lam_min_V_bound is WORST-CASE, so with lam_max(K) >> 1 it refuses
        negative components that are perfectly feasible.  So it is used only
        as a fast accept, and when it fails the true lam_min(V) is measured by
        power iteration on the composite apply.  Costs a few hundred
        single-column V applies, and only at points the bound could not
        certify -- lam_min_exact counts them, 0 on a clean replicate.

        An under-converged Rayleigh quotient OVERestimates lam_min, hence the
        extra iterations and the strictly positive margin.
        """
        margin = 1e-4 * max(s_vec[2], 1e-12)
        if lam_min_V_bound(s_vec) > margin:
            return True
        _count_event('lam_min_exact')
        s2a_, s2d_, s2e_ = s_vec
        t0 = time.perf_counter()
        lo, _ = _spectral_range(
            lambda B: _v_matvec(Z, Zd, s2a_, s2d_, s2e_, B),
            n, iters=200)
        _add_time('lam_min_exact', time.perf_counter() - t0)
        return lo > margin

    vary = y.var()
    # UPPER bound at 1.5 * var(y): the three components decompose the
    # phenotypic variance, so one approaching var(y) has diverged.
    s_upper = upper_mult * vary

    # --- LOWER bound: NO box on the two genetic components ------------------
    # Only s2e keeps a floor.  Any floor folds the sampling distribution back
    # onto itself, biasing a near-zero component UP; a negative estimate is read
    # as "0, plus the noise around it".  Nothing is un-guarded by this: what
    # holds V in the PD cone is feasible() below, not the box.  s2e keeps its
    # floor because both kernels are PSD and singular, so nothing else does.
    s_lower = np.array([-np.inf, -np.inf, 1e-9])

    rng = np.random.default_rng(seed)
    # Drawn ONCE at FULL width; the coarse phase reads U[:, :Nmc_coarse].
    # CENTERED Rademacher probes: E[uu'] = P_c, so Hutchinson estimates
    # tr(P_c V^-1 K_i) -- the RESTRICTED trace REML wants with an intercept.
    U = _center_cols(rng.choice([-1.0, 1.0], size=(n, Nmc)))
    if Nmc_coarse is None or Nmc_coarse >= Nmc:
        n_probe, coarse = Nmc, False                # single-phase fit
    else:
        n_probe, coarse = Nmc_coarse, True
    switch_it = None                                # iteration the switch fired

    # --- starting point, as FRACTIONS OF var(y) -----------------------------
    # vary * s_init keeps the estimator scale-equivariant.  The genetic
    # components start SMALL, as in the 4VC parent (0.05 each; the residual
    # takes the rest, 0.9 here against 0.85 there because there is one
    # component fewer).
    s = vary * np.asarray(s_init, dtype=float)
    if s.shape != (k,):
        raise ValueError(f"s_init must have {k} entries, got {s.shape}")
    AI = np.eye(k)

    yc = y.reshape(n, 1)
    xbuf = None       # warm-start buffers for the CG solve groups
    # Held at the FULL Nmc so the fine phase inherits the coarse phase's warm
    # starts; the untouched columns are cold-started once, at the switch.
    Pbuf = np.zeros((n, Nmc))
    Gbuf = None

    t0 = time.perf_counter()
    KaU = Kaapply(U)  # K_a U for the fixed probes: genotype-only, formed once
    KdU = Kdapply(U)  # K_d U   "        "
    _add_time('setup_KU', time.perf_counter() - t0)
    _add_time('setup', time.perf_counter() - t_setup0)
    t_phase0 = time.perf_counter()      # coarse phase starts here

    def eval_grad_ai(s_at):
        """Score and average information at s_at, at the CURRENT n_probe.

        Lifted out of the loop because the trust region needs the score at a
        trial point too, and both must be the same function of s or rho compares
        two different objectives.  Updates the warm-start buffers in place.
        """
        nonlocal xbuf, Gbuf
        s2a_, s2d_, s2e_ = s_at
        mv = lambda B: _v_matvec(Z, Zd, s2a_, s2d_, s2e_, B)

        # --- x = V^{-1} y ---
        t0 = time.perf_counter()
        xbuf = _cg_batched(mv, yc, x0=xbuf, tol=cg_tol, maxiter=cg_maxiter)
        _add_time('solve_y', time.perf_counter() - t0)
        x_ = xbuf[:, 0]

        # data quadratics x'K_i x
        Kax_, Kdx_ = Kaapply(x_), Kdapply(x_)

        # --- score traces tr(V^{-1} K_i): Hutchinson, exact probe solves ---
        Uc_ = U[:, :n_probe]
        t0 = time.perf_counter()
        Psol = _cg_batched(mv, Uc_, x0=Pbuf[:, :n_probe], tol=cg_tol,
                           maxiter=cg_maxiter)
        _add_time('solve_probe_coarse' if coarse else 'solve_probe_fine',
                  time.perf_counter() - t0)
        Pbuf[:, :n_probe] = Psol
        sc = np.array([
            0.5 * (x_ @ Kax_ - np.mean(np.sum(Psol * KaU[:, :n_probe], axis=0))),
            0.5 * (x_ @ Kdx_ - np.mean(np.sum(Psol * KdU[:, :n_probe], axis=0))),
            0.5 * (x_ @ x_ - np.mean(np.sum(Psol * Uc_, axis=0)))])

        # --- average information: A_ij = 0.5 (K_i x)' V^{-1}(K_j x) ---
        KX = np.column_stack([Kax_, Kdx_, x_])
        t0 = time.perf_counter()
        Gbuf = _cg_batched(mv, KX, x0=Gbuf, tol=cg_tol, maxiter=cg_maxiter)
        _add_time('solve_ai', time.perf_counter() - t0)
        ai = 0.5 * (KX.T @ Gbuf)
        return sc, 0.5 * (ai + ai.T)

    # Delta = Inf: the radius only exists once a step has been rejected, so a
    # replicate that never overshoots follows the plain AI-Newton path.
    Delta = np.inf
    score, AI = eval_grad_ai(s)          # the current gradient, carried forward
    n_reject = 0

    for it in range(iters):
        _count_event('reml_iters')

        # --- AI-Newton step inside the adaptive trust region ----------------
        # AI can be near-singular (K_a with K_d, K_d with I), so the LM ridge
        # handles the singularity and the trust region the overshoot.
        dA = np.abs(np.diag(AI))
        ridge = lm * (dA.mean() + 1e-12)
        # The LM ridge alone cannot make AI + ridge I strictly PD (it is scaled
        # to mean|diag(AI)|).  Under the feasibility guard this adds nothing;
        # kept so a marginal case degrades into a small step, not a wrong one.
        lam_lo = float(np.linalg.eigvalsh(AI)[0])
        if lam_lo < 0.0:
            ridge += abs(lam_lo) * (1.0 + 1e-6)
        step = np.linalg.solve(AI + (ridge + jitter) * np.eye(k), score)

        # Radius constraint on ||diag(AI) * p||: shorten the Newton step along
        # its own direction.  Inactive while Delta is Inf.
        dnorm = np.linalg.norm(np.diag(AI) * step)
        if dnorm > Delta and dnorm > 0.0:
            step *= Delta / dnorm
            dnorm = Delta

        # --- CONVERGENCE: BOLT-REML predicted log-likelihood gain ------------
        # From l(s + d) ~= l(s) + score'd - 0.5 d' AI d, dLL_pred is the model's
        # remaining height above the iterate, i.e. 0.5 sum_p ((s_opt-s)_p/SE_p)^2
        # -- so dLL_pred < tol_ll means within ~sqrt(tol_ll) SE of the optimum, a
        # statement about the ESTIMATE, not about step size.  Tested BEFORE the
        # step is taken, on the radius-constrained step, following BOLT.
        dLL_pred = float(score @ step - 0.5 * step @ (AI @ step))
        if verbose:
            print(f"iter {it:2d}  dLL_pred={dLL_pred:.6e}"
                  f"  phase={'coarse' if coarse else 'fine'}(S={n_probe})"
                  f"  Delta={Delta:.4g}", flush=True)
        # An INVARIANT given the feasibility guard (V PD => AI PSD), not a hope.
        # Checked because failing silently yields a plausible-looking number.
        assert dLL_pred > -1e-8 * max(1.0, abs(dLL_pred)), (
            f"dLL_pred={dLL_pred:.6e} < 0 at iter {it}: AI + ridge is not "
            f"positive definite despite the feasibility guard "
            f"(eigs(AI)={np.linalg.eigvalsh(AI)}, ridge={ridge:.6g}, "
            f"lam_min(V) bound={lam_min_V_bound(s):.6g})")

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
                score, AI = eval_grad_ai(s)   # re-evaluate at the full Nmc
                continue
        elif dLL_pred < tol_ll:
            break

        # --- trial point, trapezoid gain, accept / reject -------------------
        # The clamp applies to the TRIAL point, so rho uses the actual move
        # s_try - s.  FEASIBILITY: the genetic components have no lower bound, so
        # the step is backtracked until feasible() certifies it -- the only thing
        # limiting a negative component now.
        s_try = np.clip(s + step, s_lower, s_upper)
        shrink = 0
        ok = feasible(s_try)
        while not ok and shrink < 40:
            step = 0.5 * step
            s_try = np.clip(s + step, s_lower, s_upper)
            shrink += 1
            ok = feasible(s_try)
        if not ok:
            # Even a vanishing step is infeasible: the CURRENT iterate is on the
            # boundary.  Reject and shrink the radius hard.
            n_reject += 1
            Delta = alpha1 * np.linalg.norm(np.diag(AI) * step)
            if verbose:
                print(f"iter {it:2d}  INFEASIBLE even at 2^-{shrink} of the "
                      f"step (bound and exact lam_min(V) both <= 0); rejected, "
                      f"Delta->{Delta:.4g}", flush=True)
            continue
        if shrink and verbose:
            print(f"iter {it:2d}  step backtracked 2^-{shrink} to keep V "
                  f"positive definite (lam_min(V) bound="
                  f"{lam_min_V_bound(s_try):.4g})", flush=True)

        taken = s_try - s
        score_try, AI_try = eval_grad_ai(s_try)

        # Model breakdown: a gradient that GREW by more than 2x means the
        # quadratic model does not describe this region.  Reject outright.
        gnorm, gnorm_try = np.linalg.norm(score), np.linalg.norm(score_try)
        if gnorm_try > 2.0 * gnorm:
            rho = -1.0
        else:
            # TRAPEZOID rule: l(s+p) - l(s) is the line integral of the score,
            # the one computable thing (the likelihood itself needs log|V|).
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
                  f"  |D p|={taken_dnorm:.4g}  Delta->{Delta:.4g}", flush=True)

    # Phase stamp: from the switch (or from setup, single-phase) to here.
    _add_time('phase_fine', time.perf_counter() - t_phase0)
    if verbose:
        print(f"switch at iter {switch_it}; rejected {n_reject} step(s); "
              f"final Delta={Delta:.4g}; SE={reml_se(AI, Nmc)}", flush=True)
    return s, AI


def MC_REML(Z, Zd, y, iters=30, Nmc=100, Nmc_coarse=15, cg_tol=1e-6,
            cg_maxiter=1000, seed=None, verbose=False):
    """Wrapper: returns (s2a_hat, s2d_hat, s2e_hat, AI).

    SEs from reml_se(AI, Nmc).  Nmc_coarse=None gives the single-phase fit.
    """
    s, AI = mc_reml(Z, Zd, y, iters=iters, Nmc=Nmc, Nmc_coarse=Nmc_coarse,
                    cg_tol=cg_tol, cg_maxiter=cg_maxiter, seed=seed,
                    verbose=verbose)
    s2a_hat, s2d_hat, s2e_hat = s
    return s2a_hat, s2d_hat, s2e_hat, AI
