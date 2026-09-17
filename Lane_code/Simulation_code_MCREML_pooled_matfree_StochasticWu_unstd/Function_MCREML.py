# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.linalg import cholesky, eigh_tridiagonal
import time

####################################################################
# POOLED within-gene pairwise-epistasis phenotype simulation + MC AI-REML.
# UNSTANDARDIZED-KERNEL variant of Simulation_code_MCREML_pooled_matfree_
# StochasticWu.  MATRIX-FREE ESTIMATION (simulation still uses a dense W).
#
# WHAT CHANGES HERE, AND WHY IT IS THE POINT OF THIS DIRECTORY.
# The StochasticWu pipeline had a kernel MISMATCH it could not remove: its
# simulation built the STANDARDIZED kernel, h_ab = std(Z_a .* Z_b), while its
# O(n m Nw) stochastic operator can only reach the CENTERED-but-UNSCALED one,
# because the 1/sigma_ab weight is a full-rank m-by-m matrix sitting inside the
# pair sum and is exactly what forces O(n m^2).  The estimator there therefore
# carried a deterministic kernel bias on top of its Monte-Carlo error.
#
# This directory drops the standardization from BOTH SIDES.  The kernel is the
# raw within-gene interaction Gram,
#
#     h_ab = Z_a .* Z_b            (NO centering, NO 1/sigma_ab scaling)
#     W = (1/P) sum_{g=1}^G sum_{a<b in g} h_ab h_ab' ,
#     P = sum_{g=1}^G C(m_g, 2)    (total within-gene pairs, equal weight),
#
# and both build_W_pooled (simulation) and compute_WU_pooled (estimation) target
# THAT SAME W, so the ONLY error the estimator carries is Monte Carlo in the Nw
# frozen probes.  There is no kernel bias to disentangle any more, which is what
# makes a run here interpretable.
#
# THE OPERATOR (the note's "unstandardized W").  With K_w = Z Z' the additive
# GRM (n-by-n, NEVER formed) and D = Z .* Z (n-by-m),
#
#     W u = 1/(2p) [ (K_w .* K_w) u - D (D' u) ]
#         ~ 1/(2p) [ 1/Nmc ((Z(Z'(u .* (Z(Z'V))))) .* V) 1_Nmc - D (D' u) ] ,
#
# V = [v_1, ..., v_Nmc] in R^{n x Nmc} with i.i.d. Rademacher entries.  The
# SECOND line is what this module implements, and the only thing it implements:
# the stochastic operator with the diagonal estimator, one configuration, Nw the
# only knob.  The exact first line is not a code path here -- verify_unstd.py
# carries its own independent reference implementation of it, which is where a
# check belongs (a check that reuses the code under test proves nothing).
#
# WHY THE TWO LINES AGREE.  Summing over UNORDERED pairs and doubling,
#     2 sum_{a<b} (Z_a .* Z_b)(Z_a .* Z_b)'  =  sum_{a,b} (...) - sum_{a=b} (...)
#                                           =  (K_w .* K_w) - D D' ,
# since sum_{a,b} Z_ta Z_tb Z_sa Z_sb = (sum_a Z_ta Z_sa)^2 = K_ts^2.  The
# stochastic line is then the standard diagonal estimator applied to
# A = K diag(u) K:  E[v .* (A v)] = diag(A), and diag(K diag(u) K)_t
# = sum_s K_ts^2 u_s = ((K .* K) u)_t.  Exact in expectation, no bias.
#
# WHAT IS *GONE* RELATIVE TO StochasticWu.  The centering vector v_R and scalar
# s_T, and with them the whole "-s_u v_R + (s_u s_T - u'v_R) 1" tail.  They
# existed only to mean-centre the interaction columns.  Consequently:
#
#   * setup_pooled returns (V_list, P), NOT (V_list, v_R, s_T, P), and
#     compute_WU_pooled takes (Z_list, V_list, P, U).  This is a deliberate
#     signature break from the sibling pipelines -- keeping dead v_R / s_T
#     arguments that are structurally meaningless here would be worse.
#   * W 1 != 0 any more.  lambda_min(W) = 0 EXACTLY is therefore no longer
#     available, and every place that leaned on it (the SLQ conditioning guard)
#     is re-derived below from PSD-ness alone: W is a sum of outer products, so
#     lambda_min(W) >= 0 and lambda_min(V) >= s2e still bounds the conditioning
#     -- an inequality where the sibling had an identity.  The guard stays
#     conservative, which is the direction that matters.
#
# COST.  Identical to StochasticWu, since only the per-gene weight matrices were
# removed:  8 n m_g Nw per right-hand side per gene (four gemms for the
# symmetrized quartic), against 2 n m_g^2 for an exact apply.  So the operator
# is cheaper than exact when 4 Nw < m_g, its relative error is 0.24 sqrt(n/Nw)
# and FLAT in m, and it needs no m-by-m storage at all:
#
#     setup   O(n m Nw) time,  O(n m + G n Nw) storage
#
# against O(n^2) storage for the dense pre-computed W of the _preW pipeline.
# Because the error is flat in m, Nw must scale with n, not with m, to hold
# accuracy -- that is the trade the whole directory exists to measure.
#
# FROZEN PROBES.  Drawn ONCE from a seeded generator in setup_pooled and reused
# for the whole run.  Redrawing inside an apply would make W-hat a different
# matrix on every call: CG would descend on a moving target and Lanczos would
# build a basis for no operator at all.
#
# SYMMETRY, AND THE PRICE OF LOSING W 1 = 0.  At finite Nw the raw estimator's
# matrix is  A-hat_ts = (1/Nw) sum_i V_i(t) K_ts (K V_i)_s , symmetric only in
# expectation, and CG and Lanczos both require symmetry -- so the apply is
# averaged with its transpose (the SAME expression with V and Y swapped), which
# doubles the work and changes nothing else.  It is PSD only in expectation
# too, and here there is no exact null-vector anchor to certify it against, so
# slq_setup reports the smallest Ritz value and the guard uses it (see
# slq_reliable).  Raising Nw is the only remedy; there is no cheap exact
# fallback inside the REML loop.
#
# SIMULATION (build_W_pooled + simulate_Cholesky_gxg) builds the dense W ONCE to
# form the Cholesky factor Lgxg and draw a correctly-correlated epistasis
# effect.  The O(n^2) W here is intentional and unavoidable: a Cholesky factor
# needs the explicit matrix.  W is NOT cached for estimation.  Dropping the
# standardization makes this build CHEAPER than the sibling's by a factor m_g:
# the same Hadamard identity used by the operator gives W in O(n^2 m) instead of
# the O(n^2 m^2) pair loop (see build_W_pooled).
#
# THE SCORE TRACES (trace_method).  The AI-REML score needs tr(V^{-1} K_i) for
# K = (W, I).  Two estimators are provided:
#
#   'hutchinson' : the classical probe estimator, tr(V^{-1}K) ~ mean_l
#                  (V^{-1}u_l)'(K u_l), with V^{-1}U obtained by Nmc batched CG
#                  solves EVERY REML iteration.  This is the original path.
#
#   'slq'        : stochastic Lanczos quadrature (default).  For each probe, k
#                  steps of Lanczos give a Gauss rule
#                      u' f(V) u ~ ||u||^2 sum_j tau_j f(theta_j)
#                  whose nodes theta_j (Ritz values) and weights tau_j (squared
#                  first components of the tridiagonal eigenvectors) summarise
#                  the whole spectral measure of V at u.
#
# WHY SLQ IS RUN ON W, NOT ON V.  V = s2gxg W + s2e I is AFFINE in W, so both
# generate the SAME Krylov space from the same probe and the same orthonormal
# basis Q; hence exactly
#
#     T_k(V) = s2gxg T_k(W) + s2e I ,
#
# i.e. the Ritz values transform as theta_j = s2gxg mu_j + s2e and the weights
# tau_j are UNCHANGED.  Running Lanczos on W therefore yields nodes/weights that
# do not depend on the variance components at all: they are computed ONCE, before
# the REML loop, and every later iteration evaluates
#
#     tr(V^{-1} I) ~ mean_l ||u_l||^2 sum_j tau_j / (s2gxg mu_j + s2e) ,
#     tr(V^{-1} W) ~ mean_l ||u_l||^2 sum_j tau_j mu_j / (s2gxg mu_j + s2e) ,
#
# in O(Nmc k) arithmetic with NO matrix apply.  This is the real saving: the
# Hutchinson path spends Nmc CG right-hand sides per iteration (~97% of the
# solve work), the SLQ path spends k W-applies once and nothing thereafter.
#
# Note this is NOT the algebraically equivalent shortcut
# tr(V^{-1}W) = (n - s2e tr(V^{-1}))/s2gxg (exact, since s2gxg W = V - s2e I):
# that form cancels catastrophically as s2gxg -> 0, where both sides -> 0/0.
# Reading tr(V^{-1}W) off the same quadrature with an extra mu_j factor is the
# same identity evaluated stably, node by node, with no subtraction.
#
# ACCURACY.  The k-node Gauss rule is exact for polynomials of degree <= 2k-1.
# For f(x) = 1/x, f^(2k) > 0 on x > 0, so the rule is a strict LOWER bound that
# increases monotonically with k -- raise slq_k until the traces stop moving.
# Lanczos is run with full reorthogonalization, O(n k^2) per probe -- negligible
# beside the W applies.
#
# ONE SUBSTANTIVE STATISTICAL DIFFERENCE.  The unstandardized kernel weights a
# pair by its own sigma_ab^2 = Var(Z_a .* Z_b) instead of normalising it away,
# so rare-SNP pairs (large sigma_ab) dominate W and its spectrum is far more
# spread out than the standardized kernel's.  That moves W AWAY from I, which
# is the direction that helps: the near-collinearity of W and I is precisely
# what makes the 2-component AI matrix near-singular and h^2 weakly identified
# at large m.  Expect the same fit to be better conditioned here than in the
# standardized siblings at the same (n, m, G).
####################################################################


# ------------------------------------------------------------------ designs
def _standardize_cols(M, stability_std=1e-12):
    """Column-standardize to mean 0, variance 1 (ddof=0).

    A near-constant column (std < stability_std) is left un-scaled to avoid a
    divide-by-zero; such columns contribute ~0 to the GRM anyway.

    NOTE this standardizes the SNPs, which every variant does.  It is the
    standardization of the interaction columns h_ab that this directory drops.
    """
    M = np.asarray(M, dtype=float)
    mu = M.mean(axis=0)
    sd = M.std(axis=0)
    sd = np.where(sd < stability_std, 1.0, sd)
    return (M - mu) / sd


def additive_design(real_data):
    """Additive design Z_a: column-standardized allele dosages."""
    return _standardize_cols(real_data)


def split_into_genes(Z, G):
    """Split the m SNP columns of Z into G contiguous gene blocks.

    Returns a list [Z_1, ..., Z_G] of column sub-matrices -- the SEVERAL Z, one
    per gene, that the pooled kernel consumes.  np.array_split handles m not
    divisible by G by making the first (m mod G) genes one SNP larger, so the
    genes may be unequal -- exactly the case the Pooled Model is built for.
    """
    m = Z.shape[1]
    return [Z[:, cols] for cols in np.array_split(np.arange(m), G)]


# ------------------------------------- epistasis GRM (dense, SIMULATION only)
def _within_gene_sum(Zg, method='hadamard', pair_batch_size=5000):
    """UN-normalized within-gene epistasis GRM  S = sum_{a<b} h_ab h_ab',
    with the UNSTANDARDIZED interaction column  h_ab = Z_a .* Z_b.

    There is NO centering and NO 1/sigma_ab scaling -- that is the whole
    difference from the sibling pipelines -- and no 1/p_g division either: the
    Pooled Model divides ONCE by the global pair total P (see build_W_pooled),
    giving every pair equal weight.  Returns (S, p_g).

    method
    ------
    'hadamard' (default)
        The closed form, from  2 sum_{a<b} h_ab h_ab' = (K .* K) - D D'  with
        K = Z_g Z_g' and D = Z_g .* Z_g:

            S = 0.5 [ (K .* K) - D D' ] ,     O(n^2 m_g) time, O(n^2) storage.

        This is the SAME identity the matrix-free operator uses, so simulation
        and estimation are algebraically the same object by construction.
    'pairs'
        The literal C(m_g, 2)-term sum, built in pair-batches to bound memory,
        at O(n^2 m_g^2).  Kept as the independent check that 'hadamard' is the
        pair sum and not merely something close to it (verify_unstd.py asserts
        they agree to machine precision); never used by the pipeline.

    The standardized sibling has no such shortcut: its 1/sigma_ab weight does
    not factor through K, which is why it must loop over pairs.
    """
    n, mg = Zg.shape
    pg = mg * (mg - 1) // 2

    if method == 'hadamard':
        K = Zg @ Zg.T
        D = Zg * Zg
        S = 0.5 * (K * K - D @ D.T)
        return S, pg

    if method != 'pairs':
        raise ValueError(f"method must be 'hadamard' or 'pairs'; got {method!r}.")

    idx_i, idx_j = np.triu_indices(mg, k=1)
    S = np.zeros((n, n))
    for start in range(0, pg, pair_batch_size):
        end = min(start + pair_batch_size, pg)
        H = Zg[:, idx_i[start:end]] * Zg[:, idx_j[start:end]]   # raw product
        S += H @ H.T
    return S, pg


def build_W_pooled(Z_list, method='hadamard', pair_batch_size=5000):
    """Pooled WITHIN-gene pairwise-epistasis GRM, UNSTANDARDIZED interactions.

        W = (1/P) sum_{g=1}^G sum_{a<b in g} h_ab h_ab' ,   h_ab = Z_a .* Z_b ,
        P = sum_{g=1}^G C(m_g, 2) .

    A 1-SNP gene contributes no pair and is skipped.  n-by-n storage; O(n^2 m)
    time with the default 'hadamard' route.  USED ONLY BY THE SIMULATION (the
    Cholesky factor needs an explicit matrix); estimation applies the SAME W
    matrix-free via compute_WU_pooled.

    tr(W) is NOT n here.  With standardized SNPs, W_tt = (1/2P) sum_g
    [(sum_a Z_ta^2)^2 - sum_a Z_ta^4], so tr(W) = n only in expectation under
    independence and fluctuates with the genotype -- the standardized sibling's
    tr(W) = n identity does not survive dropping the 1/sigma_ab weight.
    """
    n = Z_list[0].shape[0]
    W = np.zeros((n, n))
    P = 0
    for Zg in Z_list:
        if Zg.shape[1] < 2:                  # a 1-SNP gene has no within-gene pair
            continue
        S, pg = _within_gene_sum(Zg, method=method, pair_batch_size=pair_batch_size)
        W += S
        P += pg
    if P == 0:
        raise ValueError("No gene block has >= 2 SNPs; increase m/G.")
    W /= P
    return W


# ------------------------------------- matrix-free stochastic pooled W apply
NW_DEFAULT = 200        # frozen operator probes (the note's Nmc for W u)
W_SEED_DEFAULT = 0      # seed for that draw; W-hat, hence REML, is fixed by it


def setup_pooled(Z_list, Nw=NW_DEFAULT, seed=W_SEED_DEFAULT, symmetrize=True):
    """SETUP phase: everything that depends on the genotype alone.

    Returns (V_list, P):

        V_list : per-gene state, one dict per gene (None for a 1-SNP gene)
        P      : sum_g C(m_g, 2)  -- the global pair total

    The sibling pipelines also return the pooled centering pair (v_R, s_T);
    the unstandardized kernel has NO centering term, so there is nothing to
    return and the argument is not carried downstream.  See the header.

    A gene keeps Z_g, D_g = Z_g .* Z_g, the FROZEN probes V_g, and the
    u-independent Y_g = Z_g(Z_g'V_g) = K_g V_g -- which is why each apply is
    only gemms.  Nothing m-by-m is formed, so unlike an exact route the setup
    has no m^2 term at all:

        O(n m Nw) time,  O(n m + G n Nw) storage.
    """
    rng = np.random.default_rng(seed)

    V_list = []
    P = 0
    for Zg in Z_list:
        mg = Zg.shape[1]
        if mg < 2:                               # no within-gene pair
            V_list.append(None)
            continue
        P += mg * (mg - 1) // 2

        Vp = rng.choice([-1.0, 1.0], size=(Zg.shape[0], Nw))   # frozen probes
        V_list.append({'Z': Zg, 'D': Zg * Zg, 'V': Vp,
                       'Y': Zg @ (Zg.T @ Vp),    # K_g V,  u-independent: once
                       'sym': symmetrize})

    if P == 0:
        raise ValueError("No gene block has >= 2 SNPs; increase m/G.")
    return V_list, P


def _gene_WU(gene, U, buf_elems):
    """ONE gene's un-normalized contribution,  2 S_g U, by the probe estimator.

        (K .* K) u ~ (1/Nw) sum_i V_i .* (Z(Z'(u .* Y_i))) ,   Y = K V ,

    averaged with its transpose (V and Y swapped) when symmetrizing, then minus
    D(D'u).  Columns run in blocks of w, lifting (Y .* u) into an (n, Nw*w)
    buffer so each block is two gemms (four symmetrized) independent of w:
    O(n m_g Nw) per column, nothing m-by-m.

    ON THE BLOCK SIZE.  Bigger is better, measured on the sibling pipeline at
    n=1000/m_g=100/Nw=200 over a width-100 apply:

        w =     1     2     5    10    20    25    50   100
        s   4.14  8.30  6.04  4.30  3.67  3.51  3.41  3.32

    Column-at-a-time was tried on the theory that the 32 MB buffer thrashes
    cache; it is the SECOND WORST option.  The two gemms dominate and they want
    a large right-hand side, so w is capped only by memory (buf_elems doubles,
    64 MB by default -> w = 40 here, within 6% of the unblocked optimum).
    """
    Zg, Dg, Vp, Y = gene['Z'], gene['D'], gene['V'], gene['Y']
    sym = gene['sym']
    n, c = U.shape
    Nw = Vp.shape[1]
    out = np.empty((n, c))

    cb = max(1, min(c, buf_elems // max(1, n * Nw)))
    for s in range(0, c, cb):
        e = min(s + cb, c)
        w = e - s
        Ub = U[:, s:e]

        Tb = (Y[:, :, None] * Ub[:, None, :]).reshape(n, Nw * w)   # u .* Y
        Ob = (Zg @ (Zg.T @ Tb)).reshape(n, Nw, w)                  # Z(Z'(.))
        acc = np.einsum('ni,niw->nw', Vp, Ob)                      # (.) .* V, 1'

        if sym:                          # the transpose: V and Y swapped
            Tb = (Vp[:, :, None] * Ub[:, None, :]).reshape(n, Nw * w)
            Ob = (Zg @ (Zg.T @ Tb)).reshape(n, Nw, w)
            acc = 0.5 * (acc + np.einsum('ni,niw->nw', Y, Ob))

        out[:, s:e] = acc / Nw

    out -= Dg @ (Dg.T @ U)                   # the a = b terms the pair sum drops
    return out


def compute_WU_pooled(Z_list, V_list, P, U, buf_elems=8_000_000):
    """Matrix-free product  W @ U  for the pooled UNSTANDARDIZED epistasis GRM.

    U may be (n,) or (n, c); the result matches its shape.  Never forms W
    (n-by-n), any interaction block H_g (n-by-p_g), or anything m-by-m.
    Implements

        W U = 1/(2P) sum_g 2 S_g U ,   S_g = sum_{a<b in g} h_ab h_ab' ,

    each gene's complete 2 S_g U coming from _gene_WU.  P is the one genuinely
    global quantity: the Pooled Model divides ONCE by the total pair count,
    giving every pair equal weight.  buf_elems caps the per-gene column-block
    buffer.

    W-hat is the SAME matrix on every call (the probes were frozen at setup),
    which is what lets CG and Lanczos run on it at all.

    COST vs the m-by-m-weighted exact route of the standardized siblings,
    measured at n=1000, m=1000, G=10, Nw=200 -- it depends entirely on the
    WIDTH, and the two routes scale oppositely:

        width 1     exact 0.613 s     stochastic 0.055 s     11x FASTER
        width 2     exact 0.611 s     stochastic 0.087 s      7x FASTER
        width 100   exact 1.002 s     stochastic 3.32  s      3.3x slower

    The exact route's per-apply cost is nearly FLAT in c because its m-by-m
    contraction dominates, so it is very inefficient at width 1; this operator
    has no such fixed cost and scales linearly in c.  The REML loop spends its
    CG solves at widths 1 and 2 (where the stochastic route wins) and the
    one-off SLQ Lanczos at width Nmc (where it loses), so which is faster
    overall is set by Nmc and slq_k, not by the per-apply figure alone.
    """
    n = Z_list[0].shape[0]
    U = np.asarray(U, dtype=float)
    single = (U.ndim == 1)
    if single:
        U = U.reshape(n, 1)

    T1 = np.zeros_like(U)
    for gene in V_list:
        if gene is None:                     # 1-SNP gene: no pair, no work
            continue
        T1 += _gene_WU(gene, U, buf_elems)

    out = T1 / (2.0 * P)
    return out[:, 0] if single else out


# ------------------------------------------------------------- simulation
def simulate_Cholesky_gxg(real_data, G, s2gxg=0.5, s2e=0.5, stability=1e-10):
    """Cholesky factor of the pooled epistasis covariance.

    The m SNPs of real_data are column-standardized and split into G contiguous
    genes (split_into_genes -> several Z), and the pooled within-gene kernel is
    built densely by build_W_pooled -- the UNSTANDARDIZED one, the same object
    the estimator applies.  Returns (Lgxg, w_build_time) with

        Lgxg Lgxg' = s2gxg W .

    The dense W is used here and then discarded: unlike the _preW pipeline it is
    NOT cached, because estimation rebuilds its action matrix-free from the
    genotype.  w_build_time is the wall-clock seconds spent building W only.

    W is PSD by construction but singular in general (P < n gives it rank < n,
    and it is a Gram matrix of correlated columns regardless), so `stability`
    jitter on the diagonal is what makes the factorisation succeed; it is 1e-10
    against a kernel whose diagonal is O(1), i.e. numerically invisible.
    """
    Za = additive_design(real_data)
    n, m = Za.shape

    genes = split_into_genes(Za, G)          # several Z, one per gene

    t_start = time.perf_counter()
    W = build_W_pooled(genes)
    w_build_time = time.perf_counter() - t_start

    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)
    return Lgxg, w_build_time


def simulate_remove_sampling_err(Lgxg, n, s2gxg=0.5, s2e=0.5):
    """Phenotype y = g_gxg + e with per-component sampling-error removal.

    g_gxg = Lgxg u1  (~ N(0, s2gxg W)),  e white noise.  Each component is
    rescaled to its exact target variance and y is mean-centred, so the fitted
    model carries no fixed effect.  Returns y (n,).

    ON THE MEAN-CENTRING, which is not innocuous here.  The centered siblings
    had W 1 = 0, so 1 was a null direction of the epistasis component and
    removing y's mean removed none of the signal.  The unstandardized kernel has
    W 1 != 0 (measured ||W 1||/sqrt(n) ~ 1.2 at n=300), so centring DOES discard
    a genuinely informative direction, and Var(y_centred) = M V M with
    M = I - 11'/n is not exactly the V that mc_reml fits.  The discrepancy is
    one direction out of n and is invisible at pipeline sizes -- verified by
    grid search on the exact likelihood, where mc_reml matches the true optimum
    to the Monte-Carlo trace noise (verify_unstd.py, check 6).  The centring is
    kept because it is what every sibling pipeline does, so a phenotype drawn
    here is comparable to one drawn there; dropping it would be the cleaner
    model and a break in comparability.
    """
    u1 = np.random.randn(n)
    u2 = np.random.randn(n)

    gxg = Lgxg @ u1                       # epistasis effect ~ N(0, s2gxg W)
    e = np.sqrt(s2e) * u2                 # residual noise

    # Eliminate sampling variances: rescale each component to its exact target.
    gxg *= np.sqrt(s2gxg / np.var(gxg, ddof=0))
    e *= np.sqrt(s2e / np.var(e, ddof=0))

    y = gxg + e
    y -= y.mean()
    return y


# ------------------------------------------------ matrix-free linear algebra
def _v_matvec(Z_list, V_list, P, s2gxg, s2e, B):
    """Apply V = s2gxg W + s2e I to B, with W applied matrix-free.

        V B = s2gxg (W B) + s2e B .

    B may be (n,) or (n, c); the result matches its shape.  The epistasis term
    is rebuilt from the genotype on every call -- no n-by-n W exists.
    """
    return s2gxg * compute_WU_pooled(Z_list, V_list, P, B) + s2e * B


def _cg_batched(matvec, Bmat, x0=None, tol=1e-6, maxiter=1000):
    """Conjugate gradient for the SPD system V X = Bmat.

    matvec : callable  X -> V X   (accepts / returns (n, c) arrays).
    Bmat   : (n, c) right-hand sides -- all c columns solved together with
             per-column CG scalars, so one V-pass advances every column.
    x0     : (n, c) warm start (e.g. the previous REML iteration's solution).
    """
    n, c = Bmat.shape
    X = np.zeros((n, c)) if x0 is None else x0.copy()
    R = Bmat - matvec(X)
    P = R.copy()
    rs_old = np.sum(R * R, axis=0)
    b_norm = np.sqrt(np.sum(Bmat * Bmat, axis=0))
    b_norm[b_norm == 0.0] = 1.0

    for _ in range(maxiter):
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


# --------------------------------------- stochastic Lanczos quadrature (SLQ)
def _lanczos_W_batched(Wapply, U, k):
    """k steps of Lanczos on W, run for every column of U at once.

    All c probes share ONE W-apply per step (the probes differ only in their
    scalar recurrence coefficients), so the whole run costs k applies of width c
    -- the same batching the CG solver uses.

    Reorthogonalization is FULL and applied twice ("twice is enough"): in
    floating point the Lanczos basis loses orthogonality as soon as a Ritz value
    converges, which duplicates nodes and corrupts the quadrature weights.  At
    O(n k^2) per probe this is negligible beside the W applies.

    A probe whose Krylov space is exhausted (beta ~ 0) is frozen: its basis
    vector is set to 0, and since W is linear every later alpha / beta for that
    probe is 0 too.  kmax records where that happened so the caller can
    eigendecompose only the meaningful leading block.

    Returns (alpha, beta, kmax, sq_norms) with alpha, beta of shape (k, c),
    kmax (c,) ints, sq_norms (c,) = ||u_l||^2.
    """
    n, c = U.shape
    sq_norms = np.sum(U * U, axis=0)
    nrm = np.sqrt(sq_norms)

    Q = np.zeros((k + 1, n, c))
    Q[0] = U / np.where(nrm > 0.0, nrm, 1.0)
    alpha = np.zeros((k, c))
    beta = np.zeros((k, c))
    kmax = np.full(c, k, dtype=int)

    for j in range(k):
        z = Wapply(Q[j])
        a = np.sum(Q[j] * z, axis=0)
        alpha[j] = a

        z -= a * Q[j]
        if j > 0:
            z -= beta[j - 1] * Q[j - 1]
        for _ in range(2):                       # full reorthogonalization
            coef = np.einsum('lnc,nc->lc', Q[:j + 1], z)
            z -= np.einsum('lnc,lc->nc', Q[:j + 1], coef)

        b = np.sqrt(np.sum(z * z, axis=0))
        thr = 1e-12 * np.maximum(np.abs(alpha[:j + 1]).max(axis=0), 1.0)
        broke = b <= thr                         # Krylov space exhausted
        kmax = np.where(broke & (kmax == k), j + 1, kmax)
        beta[j] = np.where(broke, 0.0, b)
        Q[j + 1] = np.where(broke, 0.0, z / np.where(broke, 1.0, b))

    return alpha, beta, kmax, sq_norms


def slq_setup(Wapply, U, k=25, buf_elems=16_000_000):
    """SLQ SETUP: Gauss-quadrature nodes and weights of W at each probe.

    Runs Lanczos on W -- NOT on V -- so the result is free of the variance
    components: V = s2gxg W + s2e I is affine in W, so the two share a Krylov
    basis and T_k(V) = s2gxg T_k(W) + s2e I exactly.  The nodes below are
    therefore Ritz values mu_j of W, mapped to V's nodes by the caller, and the
    weights are identical for both.  Done ONCE, outside the REML iteration.

    Returns (nodes, weights, sq_norms, mu_min), nodes/weights of shape (Nmc, k).
    Unused slots of a probe that broke down early carry weight 0, so they drop
    out of every quadrature sum without special-casing.

    mu_min is the smallest RAW Ritz value seen across all probes, recorded
    BEFORE the nodes are clipped at 0.  The standardized siblings did not need
    it: there W 1 = 0 held exactly on every draw, pinning lambda_min(W) = 0.
    Here the operator has no exact null vector, and the stochastic W-hat is PSD
    only in expectation, so a genuinely negative eigenvalue is possible at small
    Nw -- and it would make V = s2gxg W-hat + s2e I indefinite, which breaks CG
    silently.  mu_min is the cheap early warning (Lanczos converges to the
    extremes of the spectrum first), used by slq_reliable and reported by
    mc_reml.

    Probes are processed in chunks sized so the Lanczos basis stays near
    buf_elems doubles: peak extra memory is O(buf_elems), not O(n k Nmc).
    """
    n, c = U.shape
    nodes = np.zeros((c, k))
    weights = np.zeros((c, k))
    sq_norms = np.zeros(c)

    chunk = max(1, min(c, buf_elems // max(1, n * (k + 1))))
    for s in range(0, c, chunk):
        e = min(s + chunk, c)
        alpha, beta, kmax, sqn = _lanczos_W_batched(Wapply, U[:, s:e], k)
        sq_norms[s:e] = sqn
        for i in range(e - s):
            ki = int(kmax[i])
            if ki == 1:
                th = alpha[:1, i]
                tau = np.ones(1)
            else:
                th, Y = eigh_tridiagonal(alpha[:ki, i], beta[:ki - 1, i])
                tau = Y[0] ** 2               # squared first components
            nodes[s + i, :ki] = th
            weights[s + i, :ki] = tau

    mu_min = float(np.min(nodes[weights > 0.0])) if np.any(weights > 0.0) else 0.0

    # W is PSD by construction (and W-hat in expectation); clip Ritz values that
    # rounding or Monte-Carlo noise pushed below 0 so the denominator
    # s2gxg*mu + s2e can never be driven non-positive inside slq_traces.
    np.maximum(nodes, 0.0, out=nodes)
    return nodes, weights, sq_norms, mu_min


def slq_traces(nodes, weights, sq_norms, s2gxg, s2e):
    """Both score traces for V = s2gxg W + s2e I from one cached quadrature.

        tr(V^{-1} I) ~ mean_l ||u_l||^2 sum_j tau_j       / (s2gxg mu_j + s2e)
        tr(V^{-1} W) ~ mean_l ||u_l||^2 sum_j tau_j mu_j  / (s2gxg mu_j + s2e)

    Pure O(Nmc k) arithmetic -- no matrix apply, no linear solve.  The W trace
    carries the node mu_j in the numerator rather than being recovered from the
    identity tr(V^{-1}W) = (n - s2e tr(V^{-1}))/s2gxg, which is algebraically
    the same but cancels catastrophically as s2gxg -> 0.
    """
    denom = s2gxg * nodes + s2e
    wd = weights / denom
    trV1I = np.mean(sq_norms * np.sum(wd, axis=1))
    trV1W = np.mean(sq_norms * np.sum(wd * nodes, axis=1))
    return trV1W, trV1I


def slq_reliable(nodes, mu_min, s2gxg, s2e, k, tol=1e-3):
    """Is the CACHED k-node rule still accurate at this (s2gxg, s2e)?

    One Lanczos run has to serve every variance setting the optimizer visits,
    and a fixed k stops being enough once V becomes ill-conditioned.  That
    happens exactly as s2e -> 0 -- the boundary the box clamp permits and that a
    weakly-identified fit is known to drift onto -- so the cached quadrature
    CANNOT be trusted unconditionally.  Measured relative error of tr(V^{-1}W)
    at k=25 on the sibling pipeline: ~1e-16 at cond(V)=3, ~2e-6 at cond=130, but
    ~7e-4 at cond=16000, which is the scale of the Monte-Carlo noise itself.

    THE CONDITIONING, AND WHERE THIS DIFFERS FROM THE SIBLINGS.  There, every
    interaction column was mean-centred, so W 1 = 0 gave lambda_min(W) = 0 and
    hence lambda_min(V) = s2e EXACTLY.  The unstandardized kernel has no such
    null vector.  What survives is the inequality: W is a sum of outer products
    h_ab h_ab', so it is PSD, lambda_min(W) >= 0 and

        lambda_min(V) >= s2e ,   lambda_max(V) ~ s2gxg max_j mu_j + s2e ,

    the latter read off the cached Ritz values, which converge to the extremes
    of the spectrum first.  Using s2e as lambda_min OVERSTATES kappa whenever
    lambda_min(W) > 0, i.e. the guard is if anything more conservative here than
    in the siblings -- the safe direction.  The one case it would understate is
    a stochastic W-hat with a genuinely negative eigenvalue, so mu_min (the
    smallest raw Ritz value, pre-clipping) is folded in: a negative mu_min
    lowers the assumed lambda_min(V) and tightens the guard.

    The classical Gauss/CG bound then predicts a relative error ~ 2 rho^(2k)
    with rho = (sqrt(kappa) - 1)/(sqrt(kappa) + 1).  The bound is conservative
    by ~2 orders of magnitude, which is what you want in a guard: it errs toward
    falling back.

    Returns True when the cached rule is safe to use.
    """
    lo = max(s2gxg * min(mu_min, 0.0) + s2e, 1e-300)
    kappa = (s2gxg * float(nodes.max()) + s2e) / lo
    if not np.isfinite(kappa) or kappa <= 1.0:
        return True
    rho = (np.sqrt(kappa) - 1.0) / (np.sqrt(kappa) + 1.0)
    return 2.0 * rho ** (2 * k) <= tol


# ----------------------------------------------------------- MC AI-REML (k=2)
def mc_reml(Z, y, G, iters=30, Nmc=50, cg_tol=1e-6, cg_maxiter=1000,
            jitter=1e-8, tol=1e-8, lm=1e-3, step_frac=0.5, upper_mult=5.0,
            seed=None, verbose=False, trace_method='slq', slq_k=25,
            slq_guard_tol=1e-3, Nw=NW_DEFAULT, w_seed=W_SEED_DEFAULT):
    """Monte-Carlo AI-REML for V = s2gxg W + s2e I (pooled UNSTANDARDIZED W).

    Z is the column-standardized genotype; it is split into G contiguous genes
    and the pooled setup (per-gene state, pair total P) is built ONCE here,
    outside the iteration.

    W u is applied by the note's O(n m Nw) operator and by nothing else: the
    quartic (K .* K) u is estimated from Nw frozen Rademacher probes, drawn once
    at w_seed and reused for the whole run so W-hat is a FIXED matrix.  Relative
    error 0.24 sqrt(n/Nw), flat in m, so Nw is the accuracy knob and it tracks
    n.  verify_unstd.py holds an independent exact implementation to check this
    one against.

    trace_method
    ------------
    'slq' (default)
        Stochastic Lanczos quadrature.  slq_k Lanczos steps on W are run ONCE
        before the loop; because V is affine in W the resulting nodes/weights
        serve every variance setting, so each iteration costs 3 CG right-hand
        sides (1 for x, 2 for the AI matrix) and NO probe solves at all.
        An iteration whose (s2gxg, s2e) is too ill-conditioned for the cached
        rule (see slq_reliable) falls back to the exact probe solves, so the
        speedup is given up only where it would otherwise cost accuracy.
    'hutchinson'
        The original estimator: Nmc CG solves per iteration for V^{-1}U, with
        W U for the fixed probes formed once up front.  Kept for verification.

    Both use the SAME fixed Rademacher probes, so the objective stays a
    deterministic function of s and the iteration converges to a fixed point.
    Note SLQ replaces how each probe's quadratic form is EVALUATED, not the
    probe sampling: the Monte-Carlo error is the Hutchinson one either way, and
    slq_k controls only how exactly that same target is reproduced.

    A KNOWN DEFECT OF THE SHARED OPTIMIZER, inherited deliberately.  When a
    replicate's likelihood peaks at s2gxg = 0, s2gxg sticks on the lower box
    clamp and s2e then converges to the WRONG value: the AI-Newton step is
    solved jointly and never re-projected onto the free subspace, so the blocked
    gxg direction keeps driving s2e through the coupling term.  It settles where
    (AI^{-1} score)_e = 0 rather than where score_e = 0, i.e. above var(y).
    Measured on pure-noise phenotypes (n=200, m=80, G=4), s2e overshoots var(y)
    by up to +0.35 here and up to +0.26 in the standardized sibling, whose
    optimizer block is byte-identical -- this is NOT a property of the
    unstandardized kernel, and it does not touch s2gxg, which is 0 either way.
    It is left uncorrected so that a run here differs from a sibling run in the
    W apply ALONE; an active-set projection would fix it but would change
    exactly the replicates that make the two comparable.  verify_unstd.py check
    6 flags boundary replicates so their large grid gap is not misread as an
    operator fault.

    Returns
    -------
    s  : (2,) estimated (s2gxg, s2e).
    AI : (2, 2) final average-information matrix.
    """
    if trace_method not in ('slq', 'hutchinson'):
        raise ValueError(f"trace_method must be 'slq' or 'hutchinson'; "
                         f"got {trace_method!r}.")
    y = np.asarray(y, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float)
    n = y.shape[0]
    k = 2

    # --- setup: genotype-only, done once ---
    genes = split_into_genes(Z, G)
    V_list, P = setup_pooled(genes, Nw=Nw, seed=w_seed)
    Wapply = lambda B: compute_WU_pooled(genes, V_list, P, B)

    vary = y.var()
    s_upper = upper_mult * vary          # no component can exceed ~total var

    rng = np.random.default_rng(seed)
    U = rng.choice([-1.0, 1.0], size=(n, Nmc))      # Rademacher probes in {+-1}

    s = np.full(k, vary / k)
    AI = np.eye(k)

    yc = y.reshape(n, 1)
    xbuf = None       # warm-start buffers for the CG solve groups
    Pbuf = None
    Gbuf = None

    WU = None         # W U: needed only on the Hutchinson path, built on demand
    n_fallback = 0
    if trace_method == 'slq':
        # Lanczos on W: nodes/weights are parameter-free, so this is the ONLY
        # spectral work the whole run does.  slq_k applies of width Nmc, once.
        nodes, weights, sq_norms, mu_min = slq_setup(Wapply, U, k=slq_k)
        if verbose and mu_min < -1e-8 * max(float(nodes.max()), 1.0):
            # No W 1 = 0 anchor here, so this is the only PSD check available.
            print(f"WARNING: smallest Ritz value of W is {mu_min:.3e} < 0 -- the "
                  f"operator is not PSD at Nw={Nw}; CG may stall.  Raise Nw.")
    else:
        WU = Wapply(U)                    # W U for the fixed probes: once

    for it in range(iters):
        s2gxg, s2e = s
        matvec = lambda B: _v_matvec(genes, V_list, P, s2gxg, s2e, B)

        # --- x = V^{-1} y ---
        xbuf = _cg_batched(matvec, yc, x0=xbuf, tol=cg_tol, maxiter=cg_maxiter)
        x = xbuf[:, 0]

        # data quadratics x'K_i x
        Wx = Wapply(x)
        xWx = x @ Wx
        xIx = x @ x

        # --- score traces tr(V^{-1} K_i) ---
        use_slq = (trace_method == 'slq'
                   and slq_reliable(nodes, mu_min, s2gxg, s2e, slq_k,
                                    slq_guard_tol))
        if use_slq:
            # cached quadrature: O(Nmc slq_k) arithmetic, no solve, no apply
            trV1W, trV1I = slq_traces(nodes, weights, sq_norms, s2gxg, s2e)
        else:
            # V too ill-conditioned for the cached rule (s2e near 0): pay for
            # the exact probe solves this iteration rather than bias the score.
            if trace_method == 'slq':
                n_fallback += 1
            if WU is None:
                WU = Wapply(U)
            Pbuf = _cg_batched(matvec, U, x0=Pbuf, tol=cg_tol,
                               maxiter=cg_maxiter)
            trV1W = np.mean(np.sum(Pbuf * WU, axis=0))
            trV1I = np.mean(np.sum(Pbuf * U, axis=0))

        score = np.array([0.5 * (xWx - trV1W),
                          0.5 * (xIx - trV1I)])

        # --- average information: A_ij = 0.5 (K_i x)' V^{-1}(K_j x) ---
        KX = np.column_stack([Wx, x])
        Gbuf = _cg_batched(matvec, KX, x0=Gbuf, tol=cg_tol, maxiter=cg_maxiter)
        AI = 0.5 * (KX.T @ Gbuf)
        AI = 0.5 * (AI + AI.T)

        # --- damped, bounded AI-Newton step ---------------------------------
        # W can still be nearly collinear with I, so the AI matrix can be
        # near-singular and an undamped step explodes.  Levenberg-Marquardt
        # ridge (scaled to AI) + trust region on the step + a box clamp keep the
        # path stable without perturbing well-identified cases.  The
        # unstandardized kernel has a MORE spread spectrum than the standardized
        # one, so this damping should bind less often here -- but it is kept
        # identical so a run is comparable to the siblings' line for line.
        dA = np.abs(np.diag(AI))
        ridge = lm * (dA.mean() + 1e-12)
        step = np.linalg.solve(AI + (ridge + jitter) * np.eye(k), score)

        mx = np.abs(step).max()
        max_step = step_frac * vary                  # trust region
        if mx > max_step:
            step *= max_step / mx
        s = np.clip(s + step, 1e-9, s_upper)         # box clamp

        if verbose:
            src = 'slq' if use_slq else 'cg '
            print(f"iter {it:2d}  [{src}]  s={s}  "
                  f"max|step|={np.abs(step).max():.3e}")
        if np.abs(step).max() < tol:
            break

    if verbose and n_fallback:
        print(f"SLQ guard: {n_fallback}/{it + 1} iterations fell back to exact "
              f"probe solves (V ill-conditioned; raise slq_k to avoid).")

    return s, AI


def MC_REML(Z, y, G, iters=30, Nmc=50, cg_tol=1e-6, cg_maxiter=1000, seed=None,
            trace_method='slq', slq_k=25, Nw=NW_DEFAULT, verbose=False):
    """Wrapper: returns (s2gxg_hat, s2e_hat, AI)."""
    s, AI = mc_reml(Z, y, G, iters=iters, Nmc=Nmc, cg_tol=cg_tol,
                    cg_maxiter=cg_maxiter, seed=seed, verbose=verbose,
                    trace_method=trace_method, slq_k=slq_k, Nw=Nw)
    s2gxg_hat, s2e_hat = s
    return s2gxg_hat, s2e_hat, AI
