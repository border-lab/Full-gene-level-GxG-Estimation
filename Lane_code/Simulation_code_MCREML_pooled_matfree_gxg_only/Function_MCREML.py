# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.linalg import cholesky
from scipy.sparse.linalg import svds
import time

####################################################################
# ONE genetic variance component -- the pooled within-gene EPISTASIS kernel --
# plus noise (hence "_gxg_only").  TWO parameters:
#
#     y = g_gxg + e ,
#     V = s2gxg W + s2e I ,
#     W   = W_raw / c-hat ,
#     W_raw = (1/P) sum_g sum_{a<b in g} h_ab h_ab' ,  h_ab = Z_a .* Z_b ,
#     P = sum_g C(m_g, 2)      (UNSTANDARDIZED within-gene interactions).
#
# This is the _Lowrank_Wu_std_4VC pipeline with the ADDITIVE and DOMINANCE
# components REMOVED: no K_a, no K_d, no Z_d.  Z_a survives only as the
# standardized genotype the interactions h_ab = Z_a .* Z_b are built from.
# Everything about W is the parent's, verbatim: the c-hat normalization
# (pooled_c, C_METHOD), the gene split, the rank-r truncation of the apply,
# and the P_c wrapping.  It is the complement of the _additive sibling: that
# one keeps the cheap exact kernel and drops W, this one keeps W and drops
# the rest, so between them the cost and the estimation error of the
# epistasis component are isolated from everything else.
#
# Dividing W by c-hat puts s2gxg on the REALIZED-variance scale, so the
# estimate compares directly with its target and NO post-fit correction is
# applied anywhere.  c-hat comes from pooled_c() (C_METHOD) and both
# simulation and estimator call it on the same genotype, so the two sides
# differ by the rank-r truncation alone.
#
# EVERYTHING RUNS IN THE CONTRAST SPACE  P_c = I - 11'/n.  The phenotype, the
# component draws and the Hutchinson probes are centered, and the W apply is
# wrapped in P_c on both sides, so the epistasis kernel is effectively
# P_c W P_c.  Two payoffs: it is proper REML with an intercept (centered probes
# estimate the RESTRICTED trace tr(P_c V^-1 K_i), which supplies the -1/s2e
# correction on the residual component), and it removes the all-ones
# direction, which carries ~85% of lam_max(W) and was pinning s2gxg against
# the positive-definiteness boundary under a null.  c is defined through the
# CENTERED per-pair variance, so the normalization is unchanged.
#
# tr(W) != n: dividing by a scalar fixes only the scale.
####################################################################


# --------------------------------------------- operator-apply counting
# The estimator's whole cost is APPLIES of W-hat and I to blocks of vectors,
# almost all inside CG.  Unlike wall-clock, these counts are machine-
# independent.  Per operator: <op>_applies (calls) and <op>_columns (total
# n-vectors, the cost-bearing number).  V is the composite, W is
# compute_WU_pooled.  Also cg_solves, cg_iters, reml_iters and lam_min_exact
# (0 on a clean replicate).  GLOBAL and CUMULATIVE; mc_reml zeroes them, so
# one MC_REML call = one replicate.
#
# WALL-CLOCK TWINS.  _OP_TIMES holds seconds, keyed "<what>_sec", accumulated
# by ONE perf_counter pair per call at the same sites (V and W applies), plus
# the solve groups and setup pieces inside mc_reml and the four pieces of a
# _gene_WU call.  These ARE machine-dependent -- that is their point: counts
# say how much work was asked for, these say how fast each kind of work ran.
# Nesting: V_sec CONTAINS the W_sec spent inside V applies; W_sec contains the
# W_* pieces; every solve_*_sec contains its V applies.  So the pieces are
# quoted per column (driver), not summed.
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

    P_c = I - 11'/n applied column-wise.  The W apply is wrapped in this on
    both sides, so the model lives entirely in this space and the 1-direction
    carries no information.
    """
    return B - B.mean(axis=0, keepdims=True)


def additive_design(real_data):
    """Standardized genotype Z_a: column-standardized allele dosages.

    Here it feeds ONLY the interactions h_ab = Z_a .* Z_b; there is no K_a.
    """
    return _standardize_cols(real_data)


def split_into_genes(Z, G):
    """Split the m SNP columns of Z into G contiguous gene blocks.

    np.array_split makes the first (m mod G) genes one SNP larger.
    """
    m = Z.shape[1]
    return [Z[:, cols] for cols in np.array_split(np.arange(m), G)]


# ------------------------------------- epistasis GRM (SIMULATION only)
def _within_gene_sum(Zg):
    """UN-normalized within-gene epistasis GRM  S = sum_{a<b} h_ab h_ab', with
    the UNSTANDARDIZED  h_ab = Z_a .* Z_b.  Returns (S, p_g).

    Closed form, O(n^2 m_g) time and O(n^2) storage:

        S = 0.5 [ (K .* K) - D D' ] ,   K = Z_g Z_g' ,  D = Z_g .* Z_g ,

    the same Hadamard square the low-rank operator truncates.  No centering, no
    1/sigma_ab and no 1/p_g -- build_W_pooled divides once by the global P.
    """
    n, mg = Zg.shape
    pg = mg * (mg - 1) // 2
    K = Zg @ Zg.T
    D = Zg * Zg
    S = 0.5 * (K * K - D @ D.T)
    return S, pg


def pooled_c_exact(Z_list):
    """EXACT realized-variance factor c -- the yardstick for the O(nm) plug-ins.

        c = (1/P) sum_g sum_{a<b in g} Var-hat(Z_a .* Z_b) ,  P = sum_g C(m_g, 2)

    at ddof=0, so E[Var-hat(H gamma)] = c * s2gxg.  NOT the divisor: the kernels
    use pooled_c() (C_METHOD), and c_exact/c-hat is a run's residual scale error.

    Takes the COLUMN-STANDARDIZED gene blocks.  Formed without the n-by-P design
    H, one gemm per gene: O(n sum_g m_g^2) time, O(max_g m_g^2) storage --
    QUADRATIC in gene size, hence a diagnostic only.  1-SNP genes are skipped.
    """
    total = 0.0
    P = 0
    for Zg in Z_list:
        mg = Zg.shape[1]
        if mg < 2:                               # no within-gene pair
            continue
        n = Zg.shape[0]
        P += mg * (mg - 1) // 2
        Dg = Zg * Zg
        second = 0.5 * (np.sum(Dg.sum(axis=1) ** 2) - np.sum(Dg * Dg)) / n
        Gm = Zg.T @ Zg
        mean_sq = 0.5 * (np.sum(Gm * Gm) - np.sum(np.diag(Gm) ** 2)) / n ** 2
        total += second - mean_sq                # sum_{a<b} Var-hat(Z_a Z_b)
    if P == 0:
        raise ValueError("No gene block has >= 2 SNPs; increase m/G.")
    c = total / P
    if not (c > 0.0):
        raise ValueError(f"Pooled realized-variance factor c = {c!r} is not "
                         f"positive; the genotype has no varying interaction "
                         f"column. Check the MAF filter.")
    return c


# ---------------------------------------- c: the O(nm) plug-in estimators
# Both plug-ins are one closed form fed a different per-SNP skewness s.  From
# Var(Z_a Z_b) = 1 + r_ab s_a s_b, with R_g = Z_g'Z_g/n having unit diagonal,
#
#     c-hat = 1 + (1 / 2P) sum_g ( ||Z_g s_g||^2 / n - ||s_g||^2 ) ,
#
# i.e. ONE mat-vec per gene, O(n m_g) time, nothing m-by-m or n-by-P.  The
# routes differ only in s: third_moment_skewness reads it off the genotype,
# hwe_skewness predicts it from the allele frequency.
C_METHOD = 'moment'     # what pooled_c(), hence every kernel, uses.
                        # 'moment' | 'hwe' | 'exact'


def third_moment_skewness(Z):
    """Per-SNP skewness  s-hat_i = (1/n) sum_t Z_ti^3, with NO HWE assumption.

    On a column-standardized design the raw third moment IS the skewness.
    O(nm) time, O(m) storage.  THE PIPELINE'S ROUTE (C_METHOD).  A monomorphic
    column gets 0 and then biases c-hat, so use MAF-filtered genotypes.
    """
    Z = np.asarray(Z, dtype=float)
    return np.einsum('ij,ij,ij->j', Z, Z, Z) / Z.shape[0]


def hwe_skewness(real_data):
    """Skewness under HWE, from the allele frequency alone:

        s_i = (1 - 2 p_i) / sqrt(2 p_i (1 - p_i)) ,   p_i = mean(X_i) / 2 .

    Takes the RAW dosages.  O(nm).  DIAGNOSTIC here: its gap against
    third_moment_skewness is the HWE departure.  A monomorphic SNP gets 0.
    """
    X = np.asarray(real_data, dtype=float)
    p = X.mean(axis=0) / 2.0
    denom = 2.0 * p * (1.0 - p)
    ok = denom > 0.0
    s = np.zeros(p.shape[0])
    s[ok] = (1.0 - 2.0 * p[ok]) / np.sqrt(denom[ok])
    return s


def _split_like_genes(v, Z_list):
    """Cut a length-m per-SNP vector into the gene blocks of Z_list."""
    v = np.asarray(v, dtype=float)
    total = sum(Zg.shape[1] for Zg in Z_list)
    if v.shape[0] != total:
        raise ValueError(f"per-SNP vector has {v.shape[0]} entries but the "
                         f"gene blocks hold {total} columns.")
    out, k = [], 0
    for Zg in Z_list:
        mg = Zg.shape[1]
        out.append(v[k:k + mg])
        k += mg
    return out


def _pooled_c_plugin(Z_list, s_list, method):
    """The O(nm) closed form above, given a per-SNP skewness for each gene.

        c-hat = 1 + (1 / 2P) sum_g ( ||Z_g s_g||^2 / n - ||s_g||^2 ) .

    1-SNP genes are skipped, as everywhere else.
    """
    P = 0
    cross = 0.0
    for Zg, sg in zip(Z_list, s_list):
        mg = Zg.shape[1]
        if mg < 2:                               # no within-gene pair
            continue
        n = Zg.shape[0]
        P += mg * (mg - 1) // 2
        Zs = Zg @ sg                             # the ONE mat-vec per gene
        cross += (Zs @ Zs) / n - sg @ sg         # s' R_g s - ||s||^2
    if P == 0:
        raise ValueError("No gene block has >= 2 SNPs; increase m/G.")
    c = 1.0 + cross / (2.0 * P)
    if not (c > 0.0):
        raise ValueError(f"Plug-in realized-variance factor c-hat = {c!r} "
                         f"(method={method!r}) is not positive; the skewness "
                         f"term overwhelms the leading 1.  Check the MAF "
                         f"filter, or fall back to method='exact'.")
    return c


def pooled_c_moment(Z_list):
    """c-hat by the THIRD-MOMENT plug-in -- the divisor this pipeline uses."""
    return _pooled_c_plugin(Z_list,
                            [third_moment_skewness(Zg) for Zg in Z_list],
                            'moment')


def pooled_c_hwe(Z_list, real_data):
    """c-hat by the HWE closed form -- DIAGNOSTIC here.  Needs RAW dosages."""
    s = hwe_skewness(real_data)
    return _pooled_c_plugin(Z_list, _split_like_genes(s, Z_list), 'hwe')


def pooled_c(Z_list, method=C_METHOD, real_data=None):
    """THE DIVISOR.  Every kernel in this pipeline gets its c from here.

    Z_list is the COLUMN-STANDARDIZED gene blocks; method is

        'moment'  third-moment plug-in, O(nm)     -- the default (C_METHOD)
        'hwe'     HWE closed form,      O(nm)     -- needs real_data
        'exact'   mean per-pair sample variance, O(n sum_g m_g^2)

    build_W_pooled and setup_pooled both call this with the default on blocks
    from the same genotype, so the two sides cannot disagree.  There is no CLI
    switch for method on purpose: editing C_METHOD changes both, or neither.
    """
    if method == 'moment':
        return pooled_c_moment(Z_list)
    if method == 'exact':
        return pooled_c_exact(Z_list)
    if method == 'hwe':
        if real_data is None:
            raise ValueError("method='hwe' needs the RAW 0/1/2 dosages "
                             "(real_data=...) for the allele frequencies; the "
                             "standardized blocks no longer carry them.")
        return pooled_c_hwe(Z_list, real_data)
    raise ValueError(f"method must be 'moment', 'hwe' or 'exact'; "
                     f"got {method!r}.")


def build_W_pooled(Z_list, return_c=False):
    """Pooled within-gene epistasis GRM, UNSTANDARDIZED interactions,
    C-NORMALIZED.

        W_raw = (1/P) sum_g sum_{a<b in g} h_ab h_ab' ,  h_ab = Z_a .* Z_b ,
        W     = W_raw / c-hat ,   c-hat = pooled_c(Z_list) .

    SIMULATION ONLY (n-by-n, O(n^2 m)); the estimator applies a rank-r
    truncation of the same W via compute_WU_pooled, dividing by the same c-hat.
    c-hat is a positive scalar, so W stays symmetric PSD.  return_c=True returns
    (W, c).  1-SNP genes are skipped.  tr(W) != n and W 1 != 0.
    """
    n = Z_list[0].shape[0]
    W = np.zeros((n, n))
    P = 0
    for Zg in Z_list:
        if Zg.shape[1] < 2:                  # a 1-SNP gene has no within-gene pair
            continue
        S, pg = _within_gene_sum(Zg)
        W += S
        P += pg
    if P == 0:
        raise ValueError("No gene block has >= 2 SNPs; increase m/G.")
    c = pooled_c(Z_list)                     # C_METHOD: the third-moment c-hat
    W /= (P * c)                             # the c-normalization
    return (W, c) if return_c else W


# ---------------------------------------  low-rank pooled W apply
R_DEFAULT = 20          # truncation level r (the note's default; raise until
                        # the estimates stop moving -- LD-dependent, see header)


def setup_pooled(Z_list, r=R_DEFAULT):
    """SETUP phase: everything that depends on the genotype alone.

    Returns (F_list, P, c) -- per-gene state (None for a 1-SNP gene), the global
    pair total, and the divisor c-hat from pooled_c.  c belongs to the EXACT
    kernel, not the truncation: the phenotype came from the exact W_raw/c-hat.

    Leading eigenpairs of K_g = Z_g Z_g' come from a truncated SVD of Z_g, so
    K_g is never formed.  ARPACK needs k < min(n, m_g), so a gene where r
    reaches full rank is factorized exactly -- and the operator is then EXACT.
    """
    F_list = []
    P = 0
    for Zg in Z_list:
        mg = Zg.shape[1]
        if mg < 2:                               # no within-gene pair
            F_list.append(None)
            continue
        P += mg * (mg - 1) // 2

        full = min(Zg.shape)
        if r >= full - 1:                        # ARPACK needs k < full
            Ug, sg, _ = np.linalg.svd(Zg, full_matrices=False)
            rg = min(r, full)
        else:
            Ug, sg, _ = svds(Zg, k=r, v0=np.ones(full))
            idx = np.argsort(-sg)                # svds returns ascending
            Ug, sg = Ug[:, idx], sg[idx]
            rg = r

        F_list.append({'Z': Zg, 'D': Zg * Zg,
                       'Q': Ug[:, :rg],          # leading left singular vectors
                       'lam': sg[:rg] ** 2})     # eigenvalues of K_g

    if P == 0:
        raise ValueError("No gene block has >= 2 SNPs; increase m/G.")
    return F_list, P, pooled_c(Z_list)       # C_METHOD, the same c-hat the
                                             # simulation divided W by


def _gene_WU(gene, U, buf_elems):
    """ONE gene's un-normalized contribution  2 S_g U, by the rank-r truncation.

        (K .* K) u ~ sum_{s=1}^r lam_s q_s .* (Z(Z'(q_s .* u))) .
    """
    Zg, Dg, Q, lam = gene['Z'], gene['D'], gene['Q'], gene['lam']
    n, c = U.shape
    r = Q.shape[1]
    out = np.empty((n, c))
    Ql = Q * lam                                 # fold lam into the left factor

    # Four timers, one perf_counter stamp per boundary and ONE dict update per
    # piece per call: the gemm pair (Zg.T @ Tb, Zg @ .) against the two
    # elementwise pieces (the Q*U broadcast, the einsum + the D correction).
    t_bc = t_g1 = t_g2 = t_ew = 0.0
    cb = max(1, min(c, buf_elems // max(1, n * r)))
    for s in range(0, c, cb):
        e = min(s + cb, c)
        w = e - s
        Ub = U[:, s:e]

        t0 = time.perf_counter()
        Tb = (Q[:, :, None] * Ub[:, None, :]).reshape(n, r * w)   # q_s .* u
        t1 = time.perf_counter()
        Yb = Zg.T @ Tb                                            # (m_g, r w)
        t2 = time.perf_counter()
        Ob = (Zg @ Yb).reshape(n, r, w)                           # K(q_s .* u)
        t3 = time.perf_counter()
        out[:, s:e] = np.einsum('ns,nsw->nw', Ql, Ob)             # lam q_s .* (.)
        t4 = time.perf_counter()
        t_bc += t1 - t0
        t_g1 += t2 - t1
        t_g2 += t3 - t2
        t_ew += t4 - t3

    t0 = time.perf_counter()
    out -= Dg @ (Dg.T @ U)                   # the a = b terms the pair sum drops
    t_ew += time.perf_counter() - t0
    _add_time('W_bcast', t_bc)               # Q * U broadcast (elementwise)
    _add_time('W_gemm1', t_g1)               # Zg.T @ Tb          (gemm)
    _add_time('W_gemm2', t_g2)               # Zg @ (Zg.T @ Tb)   (gemm)
    _add_time('W_einsumD', t_ew)             # einsum + D term    (elementwise
                                             #   plus the small D gemm pair)
    return out


def compute_WU_pooled(Z_list, F_list, P, c, U, buf_elems=8_000_000):
    """Matrix-free  W-hat @ U, W-hat the rank-r truncation of the C-NORMALIZED
    pooled epistasis GRM.

        W-hat U = 1/(2 P c) sum_g 2 S-hat_g U .

    U is (n,) or (n, k) and the result matches.  Forms nothing n-by-n, n-by-p_g
    or m-by-m.  c is the same divisor the simulation used, so the two sides
    differ by the truncation alone.  Pass (F_list, P, c) from one setup_pooled.
    """
    n = Z_list[0].shape[0]
    U = np.asarray(U, dtype=float)
    single = (U.ndim == 1)
    if single:
        U = U.reshape(n, 1)
    _count_apply('W', U.shape[1])
    t0 = time.perf_counter()

    T1 = np.zeros_like(U)
    for gene in F_list:
        if gene is None:                     # 1-SNP gene: no pair, no work
            continue
        T1 += _gene_WU(gene, U, buf_elems)

    out = T1 / (2.0 * P * c)
    _add_time('W', time.perf_counter() - t0)
    return out[:, 0] if single else out


# ------------------------------------------------------------- simulation
def simulate_Cholesky_gxg(real_data, G, s2gxg=0.1, s2e=0.9, stability=1e-10):
    """Cholesky factor of the epistasis covariance.

        Lgxg Lgxg' = s2gxg W_raw/c .

    Returns (Lgxg, w_build_time, c); w_build_time covers the kernel build (the
    O(n^2 m) Hadamard squares plus the O(nm) c-hat), not the Cholesky.  The
    dense W is freed once its factor exists.  This is the memory-critical job
    of the pipeline.  s2e is accepted for a uniform signature with the sibling
    pipelines; it plays no part here.
    """
    Za = additive_design(real_data)
    genes = split_into_genes(Za, G)          # several Z, one per gene

    # The c-normalization is inside build_W_pooled and inside the timed block;
    # it is one gemm per gene, so the timing barely moves.
    t_start = time.perf_counter()
    W, c = build_W_pooled(genes, return_c=True)
    w_build_time = time.perf_counter() - t_start

    n = Za.shape[0]
    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)
    del W
    return Lgxg, w_build_time, c


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


def simulate_phenotype(Lgxg, n, s2gxg=0.1, s2e=0.9,
                       return_realized=False, force_realized=True):
    """Draw one phenotype  y = g_gxg + e,

        g_gxg = P_c Lgxg u1 ,    e = P_c sqrt(s2e) u2 ,

    the two draws independent and CENTERED.  Centering a draw from W is, in
    distribution, a draw from P_c W P_c -- the kernel the estimator applies --
    so the Cholesky factor is unchanged.  Returns y, or (y, V_ell, V_e) with
    return_realized=True; all Var-hat at ddof=0.

    Because Lgxg Lgxg' = (s2gxg/(P c)) H H', g_gxg has exactly the law of
    H gamma with gamma ~ N(0, (s2gxg/(P c)) I_P) without ever forming H.

    Realized variances are RANDOM.  For the raw kernel E[Var-hat(g_gxg)] =
    c * s2gxg; here W carries 1/c-hat, so it is (c/c-hat) s2gxg -- i.e. s2gxg to
    the plug-in's accuracy.

    force_realized=True (THE DEFAULT, as in the parents) rescales each
    component by sqrt(target / Var-hat) so its realized variance hits the
    target exactly.  An estimate's deviation from the target is then the
    ESTIMATOR's error alone and the residual c/c-hat factor is absorbed -- at
    the price that y no longer has the Gaussian law REML fits (dividing by a
    random sample variance): a MISSPECIFICATION, small at these n.  A
    component with target 0 is left at zero.  force_realized=False restores
    the plain draw.

    NOTE the two do NOT sum to var(y) under either setting: the sample
    cross-product is not zero, and rescaling touches only the diagonal terms.
    """
    u1 = np.random.randn(n)
    u2 = np.random.randn(n)

    gxg = _center_cols(Lgxg @ u1)         # epistasis  ~ N(0, s2gxg P_c W P_c)
    e = _center_cols(np.sqrt(s2e) * u2)   # residual noise

    if force_realized:
        gxg = _force_var(gxg, s2gxg, 'gxg')
        e = _force_var(e, s2e, 'e')

    y = gxg + e
    if return_realized:
        return y, float(gxg.var()), float(e.var())
    return y


# ------------------------------------ realized-variance scale factor  c
def compute_c_pooled(real_data, G, method=C_METHOD):
    """Realized-variance factor c of the RAW pooled kernel, from raw dosages.

    pooled_c() with the standardize-and-split done for you: builds Z_a, cuts it
    into the same G genes, and dispatches on method ('moment' the default,
    'hwe', 'exact').  Takes RAW dosages because the HWE route needs the allele
    frequencies.  The Cholesky job calls all three and writes them to
    result/c_<FILENAME>.txt.

    NOT where the kernel gets its divisor -- that calls pooled_c directly.
    """
    Z = additive_design(real_data)
    genes = split_into_genes(Z, G)
    return pooled_c(genes, method=method, real_data=real_data)


# ------------------------------------------------ Linear operator for the 2-parameter system
def _v_matvec(Z_list, F_list, P, c, s2gxg, s2e, B):
    """Apply  V = s2gxg W-hat + s2e I  to B, matrix-free.

    B is (n,) or (n, k) and the result matches.  No n-by-n GRM exists: the
    epistasis term is rebuilt from the cached low-rank factors and wrapped in
    P_c on both sides.  THE UNIT OF WORK: CG sees V and nothing else; a
    centered right-hand side stays centered, since every term preserves it.
    """
    _count_apply('V', 1 if np.ndim(B) == 1 else np.shape(B)[1])
    t0 = time.perf_counter()
    out = (s2gxg * _center_cols(
               compute_WU_pooled(Z_list, F_list, P, c, _center_cols(B)))
           + s2e * B)
    _add_time('V', time.perf_counter() - t0)    # includes its W apply
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

    NOT a PSD routine, deliberately: the rank-r W-hat is symmetric but can be
    indefinite, and under linkage equilibrium at small r its negative end
    dominates (r = 20, m = 1000, G = 10: -5.73 against +0.51).

    Two shifted power iterations on ONE column: the dominant-MAGNITUDE
    eigenvalue mu1 as a SIGNED Rayleigh quotient, then the dominant of
    A - mu1 I, which is the opposite extreme.  Deterministic (fixed seed).
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


# ----------------------------------------------------------- MC AI-REML (k=2)
def mc_reml(Z, y, G, iters=30, Nmc=100, Nmc_coarse=15, cg_tol=1e-6,
            cg_maxiter=1000, jitter=1e-8, tol_ll=1e-4, tol_ll_coarse=1e-2,
            lm=1e-3, eta1=1e-4, eta2=0.99, alpha1=0.25, alpha2=3.5,
            upper_mult=1.5, s_init=(0.05, 0.95),
            seed=None, verbose=False, r=R_DEFAULT):
    """Monte-Carlo AI-REML for  V = s2gxg W + s2e I, with W = W_raw/c applied
    by its rank-r truncation.

    The fit runs in the CONTRAST SPACE: y and the probes are centered and the
    W apply is wrapped in P_c, so the epistasis kernel is P_c W P_c.  That makes
    this proper REML with an intercept and keeps lam_max(W) off the all-ones
    direction, where most of it otherwise sits.

    Z is the standardized genotype, split into G genes for the epistasis
    setup; it is the only large array held.  Because the kernel carries
    1/c-hat, s2gxg is the REALIZED epistasis variance and needs no post-fit
    correction.  W-hat is deterministic in (genotype, r) and its error is a
    truncation BIAS -- small under LD, large under linkage equilibrium,
    reduced only by raising r.  With no other kernel in the model, that bias
    lands on s2gxg-hat and s2e-hat alone, uncoupled from anything else.

    Score traces
    ------------
    tr(V^{-1} W) by HUTCHINSON: Nmc fixed Rademacher probes, W U formed once,
    Nmc warm-started CG solves per iteration.  Fixed probes make the objective
    a deterministic function of s.

    TWO-PHASE SCHEDULE (BOLT-REML note 3.5): Nmc_coarse probes to
    tol_ll_coarse, then all Nmc to tol_ll.  The probes are drawn ONCE at full
    width and the coarse phase uses the first Nmc_coarse columns, so the switch
    EXTENDS the objective rather than replacing it.  Nmc_coarse=None gives the
    single-phase fit.

    Positive definiteness
    ---------------------
    The genetic component is UNBOUNDED BELOW (only s2e keeps a floor), so the
    iterate can leave the PD cone -- and not only through a negative
    component: W-hat is not PSD, so a POSITIVE s2gxg can do it too.
    feasible() therefore certifies every trial point, cheap Weyl bound first
    and true lam_min(V) when that fails; an uncertifiable step is backtracked,
    then rejected.

    OPTIMIZER.  The block below is the 4VC parent's, with k = 2 instead of 4
    and the additive and dominance kernels removed from V and from the
    feasibility bound -- nothing else changed, so a run here differs from a
    4VC run in the KERNELS alone.  It carries the parent's KNOWN DEFECT,
    inherited deliberately: if a replicate's likelihood peaks at a component
    = 0 (s2e has a floor), that component sticks and the other converges to
    the WRONG value -- the Newton step is never re-projected onto the free
    subspace.  calc_stats.py counts the affected replicates.

    Returns
    -------
    s  : (2,) estimated (s2gxg, s2e).
    AI : (2, 2) average information at the FULL Nmc.  SEs via reml_se(AI, Nmc);
         a raw sqrt(diag(inv(AI))) understates them.
    """
    reset_op_counts()                  # counts describe THIS replicate alone
    y = np.asarray(y, dtype=float).flatten()
    y = y - y.mean()                   # work in the contrast space throughout
    Z = np.asarray(Z, dtype=float)
    n = y.shape[0]
    k = 2

    # --- setup: genotype-only, done once ---
    # Timed in three pieces (svd factors incl. c-hat, spectral range, the
    # fixed-probe product W U) plus their total, setup_sec.
    t_setup0 = time.perf_counter()
    genes = split_into_genes(Z, G)
    F_list, P, c = setup_pooled(genes, r=r)      # c: the normalization divisor
    _add_time('setup_svd', time.perf_counter() - t_setup0)
    Wapply = lambda B: _center_cols(
        compute_WU_pooled(genes, F_list, P, c, _center_cols(B)))

    # Spectral range for the feasibility bound, genotype-only, computed ONCE.
    # W-hat is NOT PSD, so both its ends are estimated and the lower one may
    # be negative.
    t0 = time.perf_counter()
    lmin_w, lmax_w = _spectral_range(Wapply, n)
    _add_time('setup_spectral', time.perf_counter() - t0)
    lmin_k = np.array([lmin_w])
    lmax_k = np.array([lmax_w])

    def lam_min_V_bound(s_vec):
        """Worst-case lower bound on lam_min(V).  Positive => V is PD.

        s_i lam_min(W-hat) when s_i > 0, s_i lam_max(W-hat) when s_i < 0.
        Using lam_max for both signs would assume W-hat PSD and silently
        un-guard the direction where a POSITIVE s2gxg times a negative
        lam_min(W-hat) takes V indefinite.
        """
        g = np.asarray(s_vec[:1], dtype=float)
        worst = np.where(g > 0.0, g * lmin_k, g * lmax_k)
        return float(s_vec[1] + worst.sum())

    def feasible(s_vec):
        """Is V at s_vec positive definite?  Cheap bound first, exact if it fails.

        lam_min_V_bound is WORST-CASE, so with lam_max(W-hat) >> 1 it refuses
        negative components that are perfectly feasible.  So it is used only
        as a fast accept, and when it fails the true lam_min(V) is measured by
        power iteration on the composite apply.  Costs a few hundred
        single-column V applies, and only at points the bound could not
        certify -- lam_min_exact counts them, 0 on a clean replicate.

        An under-converged Rayleigh quotient OVERestimates lam_min, hence the
        extra iterations and the strictly positive margin.
        """
        margin = 1e-4 * max(s_vec[1], 1e-12)
        if lam_min_V_bound(s_vec) > margin:
            return True
        _count_event('lam_min_exact')
        s2gxg_, s2e_ = s_vec
        t0 = time.perf_counter()
        lo, _ = _spectral_range(
            lambda B: _v_matvec(genes, F_list, P, c, s2gxg_, s2e_, B),
            n, iters=200)
        _add_time('lam_min_exact', time.perf_counter() - t0)
        return lo > margin

    vary = y.var()
    # UPPER bound at 1.5 * var(y): the two components decompose the
    # phenotypic variance, so one approaching var(y) has diverged.
    s_upper = upper_mult * vary

    # --- LOWER bound: NO box on the genetic component ---------------------
    # Only s2e keeps a floor.  Any floor folds the sampling distribution back
    # onto itself, biasing a near-zero component UP; a negative estimate is read
    # as "0, plus the noise around it".  Nothing is un-guarded by this: what
    # holds V in the PD cone is feasible() below, not the box.  s2e keeps its
    # floor because W is singular, so nothing else does.
    s_lower = np.array([-np.inf, 1e-9])

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
    # component starts SMALL, as in the parents (0.05; the residual takes the
    # rest, 0.95 here against 0.85 in 4VC because there are fewer components).
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
    WU = Wapply(U)    # W U for the fixed probes: genotype-only, formed once
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
        s2gxg_, s2e_ = s_at
        mv = lambda B: _v_matvec(genes, F_list, P, c, s2gxg_, s2e_, B)

        # --- x = V^{-1} y ---
        t0 = time.perf_counter()
        xbuf = _cg_batched(mv, yc, x0=xbuf, tol=cg_tol, maxiter=cg_maxiter)
        _add_time('solve_y', time.perf_counter() - t0)
        x_ = xbuf[:, 0]

        # data quadratic x'W x
        Wx_ = Wapply(x_)

        # --- score traces tr(V^{-1} K_i): Hutchinson, exact probe solves ---
        Uc_ = U[:, :n_probe]
        t0 = time.perf_counter()
        Psol = _cg_batched(mv, Uc_, x0=Pbuf[:, :n_probe], tol=cg_tol,
                           maxiter=cg_maxiter)
        _add_time('solve_probe_coarse' if coarse else 'solve_probe_fine',
                  time.perf_counter() - t0)
        Pbuf[:, :n_probe] = Psol
        sc = np.array([
            0.5 * (x_ @ Wx_ - np.mean(np.sum(Psol * WU[:, :n_probe], axis=0))),
            0.5 * (x_ @ x_ - np.mean(np.sum(Psol * Uc_, axis=0)))])

        # --- average information: A_ij = 0.5 (K_i x)' V^{-1}(K_j x) ---
        KX = np.column_stack([Wx_, x_])
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
        # AI can be near-singular (W with I), so the LM ridge handles the
        # singularity and the trust region the overshoot.
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
        # s_try - s.  FEASIBILITY: the genetic component has no lower bound, so
        # the step is backtracked until feasible() certifies it -- the only
        # thing limiting a negative component now.
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


def MC_REML(Z, y, G, iters=30, Nmc=100, Nmc_coarse=15, cg_tol=1e-6,
            cg_maxiter=1000, seed=None, r=R_DEFAULT, verbose=False):
    """Wrapper: returns (s2gxg_hat, s2e_hat, AI).

    SEs from reml_se(AI, Nmc).  Nmc_coarse=None gives the single-phase fit.
    """
    s, AI = mc_reml(Z, y, G, iters=iters, Nmc=Nmc, Nmc_coarse=Nmc_coarse,
                    cg_tol=cg_tol, cg_maxiter=cg_maxiter, seed=seed,
                    verbose=verbose, r=r)
    s2gxg_hat, s2e_hat = s
    return s2gxg_hat, s2e_hat, AI
