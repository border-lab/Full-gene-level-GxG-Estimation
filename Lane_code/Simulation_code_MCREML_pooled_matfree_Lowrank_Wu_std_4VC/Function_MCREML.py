# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.linalg import cholesky
from scipy.sparse.linalg import svds
from scipy.optimize import minimize
import time

####################################################################
# FOUR variance components: additive + dominance + pooled within-gene
# epistasis + noise.
#
#     y = g_a + g_d + g_gxg + e ,
#     V = s2a K_a + s2d K_d + s2gxg W + s2e I ,
#     K_a = Z_a Z_a'/m ,  K_d = Z_d Z_d'/m      (standardized designs),
#     W   = W_raw / c-hat ,
#     W_raw = (1/P) sum_g sum_{a<b in g} h_ab h_ab' ,  h_ab = Z_a .* Z_b ,
#     P = sum_g C(m_g, 2)      (UNSTANDARDIZED within-gene interactions).
#
# Dividing W by c-hat is what separates this pipeline from _unstd_4VC, which
# corrects after the fit instead.  It puts s2gxg on the REALIZED-variance
# scale, so all four estimates compare directly with their targets and NO
# post-fit correction is applied anywhere.  c-hat comes from pooled_c() (the
# third-moment plug-in) and both simulation and estimator call it on the same
# genotype, so the two sides differ by the rank-r truncation alone.
#
# EVERYTHING RUNS IN THE CONTRAST SPACE  P_c = I - 11'/n.  The phenotype, the
# component draws and the Hutchinson probes are centered, and the W apply is
# wrapped in P_c on both sides, so the epistasis kernel is effectively
# P_c W P_c.  Two payoffs: it is proper REML with an intercept (centered probes
# estimate the RESTRICTED trace tr(P_c V^-1 K_i), which supplies the -1/s2e
# correction on the residual component), and it removes the all-ones direction,
# which carries ~85% of lam_max(W).  K_a and K_d need no wrapping:
# their designs are column-standardized, so K 1 = 0 already.  c is defined
# through the CENTERED per-pair variance, so the normalization is unchanged.
#
# THE OPTIMIZER is BOLT-REML's (Loh et al. 2015, Supplementary Note 3.2-3.5):
# Monte-Carlo AI-REML with fixed Hutchinson probes, a two-phase probe schedule,
# and each step the solution of a small QP -- the AI quadratic model maximized
# over the box s2a, s2d, s2gxg >= 0, s2e >= 1e-9 var(y), inside an adaptive
# trust region.  A component whose likelihood peaks at 0 sits on its bound
# while the others converge; its SE is reported as NaN.  W-hat is a truncation
# and not exactly PSD, so CG checks p'Vp > 0 and a trial point where it fails
# is rejected (cg_negcurv).
#
# W-HAT HAS TWO APPLY ROUTES (mc_reml's w_route), the same truncation and the
# same arithmetic to round-off.  STORED-A (default) builds At = [sqrt(lam_s)
# q_s .* Z_g] once and applies two gemm pairs: faster on wide applies, but it
# costs n * r * m_W memory (m_W = SNPs in genes) and is sensitive to memory
# bandwidth.  BROADCAST works gene by gene from the SVD factors, scaling U by
# q_s: no extra memory, and more robust on a crowded node.
#
# tr(W) != n: dividing by a scalar fixes only the scale.
# h_ab pairs ADDITIVE columns; dominance enters through K_d alone.
####################################################################


# --------------------------------------------- operator-apply counting
# The estimator's whole cost is APPLIES of K_a, K_d, W-hat and I to blocks of
# vectors, almost all inside CG.  Unlike wall-clock, these counts are
# machine-independent.  Per operator: <op>_applies (calls) and <op>_columns
# (total n-vectors, the cost-bearing number).  V is the composite, K is
# compute_KU serving both designs, W is compute_WU_storedA or
# compute_WU_bcast (whichever route the fit uses).  Also cg_solves,
# cg_iters, reml_iters and cg_negcurv (0 on a clean replicate).  GLOBAL and
# CUMULATIVE; mc_reml zeroes them, so one MC_REML call = one replicate.
#
# WALL-CLOCK TWINS.  _OP_TIMES holds seconds, keyed "<what>_sec", accumulated
# by ONE perf_counter pair per call at the same sites (V, K, W applies), plus
# the solve groups and setup pieces inside mc_reml and the pieces of the W
# apply: W_gemm_A / W_gemm_D on the stored-A route, W_bcast / W_gemm1 /
# W_gemm2 / W_einsumD on the broadcast route.  These ARE machine-dependent: counts
# say how much work was asked for, these say how fast each kind of work ran.
# Nesting: V_sec CONTAINS the K_sec and W_sec spent inside V applies; W_sec
# contains the W_* pieces; every solve_*_sec contains its V applies.  So the
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

    P_c = I - 11'/n applied column-wise.  Every kernel here already annihilates
    1 (K_a and K_d because their designs are column-standardized, W because the
    apply is wrapped), so the model lives entirely in this space and the
    1-direction carries no information.
    """
    return B - B.mean(axis=0, keepdims=True)


def additive_design(real_data):
    """Additive design Z_a: column-standardized allele dosages.

    Feeds both K_a = Z_a Z_a'/m and the interactions h_ab = Z_a .* Z_b.
    """
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


def split_into_genes(Z, G):
    """Split the m SNP columns of Z into G contiguous gene blocks.

    np.array_split makes the first (m mod G) genes one SNP larger.  K_a does not
    use this split; only the epistasis kernel does.
    """
    m = Z.shape[1]
    return [Z[:, cols] for cols in np.array_split(np.arange(m), G)]


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


# ---------------------------------------- c-hat: the third-moment plug-in
# c-hat puts s2gxg on the realized-variance scale.  From
# Var(Z_a Z_b) = 1 + r_ab s_a s_b, with R_g = Z_g'Z_g/n having unit diagonal,
#
#     c-hat = 1 + (1 / 2P) sum_g ( ||Z_g s_g||^2 / n - ||s_g||^2 ) ,
#
# i.e. ONE mat-vec per gene, O(n m_g) time, nothing m-by-m or n-by-P, with the
# skewness s_g read off the genotype by third_moment_skewness.

def third_moment_skewness(Z):
    """Per-SNP skewness  s-hat_i = (1/n) sum_t Z_ti^3, with NO HWE assumption.

    On a column-standardized design the raw third moment IS the skewness.
    O(nm) time, O(m) storage.  A monomorphic column gets 0 and then biases
    c-hat, so use MAF-filtered genotypes.
    """
    Z = np.asarray(Z, dtype=float)
    return np.einsum('ij,ij,ij->j', Z, Z, Z) / Z.shape[0]


def pooled_c(Z_list):
    """THE DIVISOR c-hat, by the third-moment plug-in.  Every kernel gets it here.

    Takes the COLUMN-STANDARDIZED gene blocks, so the divisor is a deterministic
    function of exactly what is being normalized.  build_W_pooled and
    setup_pooled both call this on blocks from the same genotype, so the two
    sides cannot disagree.  1-SNP genes are skipped, as everywhere else.
    """
    P = 0
    cross = 0.0
    for Zg in Z_list:
        mg = Zg.shape[1]
        if mg < 2:                               # no within-gene pair
            continue
        n = Zg.shape[0]
        P += mg * (mg - 1) // 2
        sg = third_moment_skewness(Zg)
        Zs = Zg @ sg                             # the ONE mat-vec per gene
        cross += (Zs @ Zs) / n - sg @ sg         # s' R_g s - ||s||^2
    if P == 0:
        raise ValueError("No gene block has >= 2 SNPs; increase m/G.")
    c = 1.0 + cross / (2.0 * P)
    if not (c > 0.0):
        # c-hat is an estimate, not an average of variances, so it CAN go
        # non-positive; dividing by it would flip the epistasis component's sign.
        raise ValueError(f"Plug-in realized-variance factor c-hat = {c!r} is "
                         f"not positive; the skewness term overwhelms the "
                         f"leading 1.  Check the MAF filter.")
    return c


def build_W_pooled(Z_list, return_c=False):
    """Pooled within-gene epistasis GRM, UNSTANDARDIZED interactions,
    C-NORMALIZED.

        W_raw = (1/P) sum_g sum_{a<b in g} h_ab h_ab' ,  h_ab = Z_a .* Z_b ,
        W     = W_raw / c-hat ,   c-hat = pooled_c(Z_list) .

    SIMULATION ONLY (n-by-n, O(n^2 m)); the estimator applies a rank-r
    truncation of the same W via compute_WU_storedA, dividing by the same c-hat.
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
    c = pooled_c(Z_list)                     # the third-moment c-hat
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
    return F_list, P, pooled_c(Z_list)       # the same c-hat the
                                             # simulation divided W by


# --- THE W-HAT APPLY: STORED A -----------------------------------------------
# Per gene, the rank-r truncation is
#
#     2 S-hat_g U = sum_s lam_s A_s (A_s' U) - D_g (D_g' U) ,   A_s = q_s .* Z_g .
#
# Built ONCE with sqrt(lam) folded in (lam_s = s_s^2 >= 0),
#
#     At_g = [A_1 | ... | A_r] Lam_g^{1/2} ,   Lam_g = blockdiag(lam_s I_{m_g}) ,
#
# gives 2 S-hat_g U = At_g (At_g' U) - D_g (D_g' U), and stacking every gene
# horizontally, At = [At_1 | ... | At_G], D = [D_1 | ... | D_G], turns the sum
# over genes into the gemm's own inner sum:
#
#     W-hat U = ( At (At' U) - D (D' U) ) / (2 P c) .
#
# Two gemm pairs, no per-gene loop, no per-apply build -- paid for in memory:
# At is (n, sum_g r_g m_g), i.e. r times Z when every SNP is in a gene.
def setup_storedA(F_list, dtype=np.float64, max_gb=16.0):
    """Build (At, D) for compute_WU_storedA from one setup_pooled's F_list.

    At : (n, sum_g r_g m_g), gene g's block ordered s-major, sqrt(lam_s) q_s .* Z_g
    D  : (n, sum_g m_g), the D_g = Z_g .* Z_g side by side
    Written block by block into preallocated arrays (no stacking copies).  1-SNP
    genes (None) carry no pair and are skipped.  The
    size of At + D is printed and checked against max_gb BEFORE allocating.
    """
    t0 = time.perf_counter()
    dtype = np.dtype(dtype)
    genes = [g for g in F_list if g is not None]
    n = genes[0]['Z'].shape[0]
    ncol_A = sum(g['Q'].shape[1] * g['Z'].shape[1] for g in genes)
    ncol_D = sum(g['Z'].shape[1] for g in genes)
    gb = n * (ncol_A + ncol_D) * dtype.itemsize / 1e9
    print(f"stored-A: At ({n} x {ncol_A}) + D ({n} x {ncol_D}) in {dtype.name} "
          f"= {gb:.3f} GB (limit {max_gb:g} GB)", flush=True)
    if gb > max_gb:
        raise MemoryError(
            f"the stored-A W apply needs {gb:.3f} GB for At + D, over the limit of "
            f"{max_gb:g} GB (--max_A_gb).  Use --A_dtype float32, a "
            f"smaller r, or raise --max_A_gb.")

    At = np.empty((n, ncol_A), dtype=dtype)
    D = np.empty((n, ncol_D), dtype=dtype)
    ja = jd = 0
    for g in genes:
        Zg, Q, lam = g['Z'], g['Q'], g['lam']
        mg = Zg.shape[1]
        sq = np.sqrt(lam)
        for s in range(Q.shape[1]):              # one (n, m_g) block per s,
            np.multiply((sq[s] * Q[:, s])[:, None], Zg,   # straight into At
                        out=At[:, ja:ja + mg])
            ja += mg
        D[:, jd:jd + mg] = g['D']
        jd += mg
    _add_time('setup_A', time.perf_counter() - t0)
    return At, D


def compute_WU_storedA(At, D, P, c, U):
    """Matrix-free  W-hat @ U  from the stored (At, D) of setup_storedA.

        W-hat U = ( At (At' U) - D (D' U) ) / (2 P c) .

    U is (n,) or (n, k) and the result
    matches, always float64.  With float32 At / D, U is cast to float32 for the
    two gemm pairs (else numpy would upcast At on every apply) and the result
    cast back.
    """
    n = At.shape[0]
    U = np.asarray(U, dtype=float)
    single = (U.ndim == 1)
    if single:
        U = U.reshape(n, 1)
    _count_apply('W', U.shape[1])
    t0 = time.perf_counter()
    Uc = U.astype(At.dtype, copy=False)
    t1 = time.perf_counter()
    out = At @ (At.T @ Uc)                       # sum_g 2 x pair part
    t2 = time.perf_counter()
    out -= D @ (D.T @ Uc)                        # the a = b terms
    t3 = time.perf_counter()
    out = out.astype(np.float64, copy=False) / (2.0 * P * c)
    _add_time('W_gemm_A', t2 - t1)
    _add_time('W_gemm_D', t3 - t2)
    _add_time('W', time.perf_counter() - t0)
    return out[:, 0] if single else out


# --- THE W-HAT APPLY: BROADCAST ----------------------------------------------
# The same rank-r truncation, gene by gene, with the r diagonal scalings q_s
# applied to U instead of Z_g:
#
#     2 S-hat_g U = sum_s lam_s q_s .* (Z_g (Z_g' (q_s .* U))) - D_g (D_g' U) .
#
# Nothing is stored beyond setup_pooled's per-gene factors; the n-by-(r w)
# intermediate q_s .* U is rebuilt on every apply, blocked over columns so it
# stays under buf_elems.
def _gene_WU_bcast(gene, U, buf_elems):
    """ONE gene's un-normalized contribution  2 S_g U, by the rank-r truncation,
    scaling U (the BROADCAST route).

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


def compute_WU_bcast(F_list, P, c, U, buf_elems=8_000_000):
    """Matrix-free  W-hat @ U  by the broadcast route, from setup_pooled's F_list.

        W-hat U = 1/(2 P c) sum_g 2 S-hat_g U .

    The same W-hat as compute_WU_storedA.  U is (n,) or (n, k) and the result
    matches.  1-SNP genes (None) carry no pair and are skipped.
    """
    genes = [g for g in F_list if g is not None]
    n = genes[0]['Z'].shape[0]
    U = np.asarray(U, dtype=float)
    single = (U.ndim == 1)
    if single:
        U = U.reshape(n, 1)
    _count_apply('W', U.shape[1])
    t0 = time.perf_counter()
    out = np.zeros_like(U)
    for gene in genes:
        out += _gene_WU_bcast(gene, U, buf_elems)
    out /= (2.0 * P * c)
    _add_time('W', time.perf_counter() - t0)
    return out[:, 0] if single else out


W_ROUTES = ('storedA', 'bcast')


# ------------------------------------------------------------- simulation
def simulate_Cholesky_4vc(real_data, G, s2a=0.1, s2d=0.1, s2gxg=0.1, s2e=0.7,
                          stability=1e-10):
    """Cholesky factors of the three genetic covariances.

        La La' = s2a K_a ,  Ld Ld' = s2d K_d ,  Lgxg Lgxg' = s2gxg W_raw/c .

    Returns (La, Ld, Lgxg, w_build_time, c); w_build_time covers the EPISTASIS
    kernel alone, the O(n^2 m) term that dominates this job.

    Built SEQUENTIALLY, each dense GRM freed once its factor exists, so the peak
    is ~4 n-by-n arrays rather than 7 (~8 GB at n = 16000).  This is the
    memory-critical job of the pipeline.
    """
    Za = additive_design(real_data)
    Zd = dominance_design(real_data)
    n, m = Za.shape

    # --- epistasis first (the expensive one), then free W -------------------
    genes = split_into_genes(Za, G)          # several Z, one per gene

    # The c-normalization is inside build_W_pooled and inside the timed block;
    # it is one gemm per gene, so the timing barely moves.
    t_start = time.perf_counter()
    W, c = build_W_pooled(genes, return_c=True)
    w_build_time = time.perf_counter() - t_start

    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)
    del W

    # --- additive --------------------------------------------------------
    Ka = build_K(Za)
    La = cholesky(s2a * Ka + stability * np.eye(n), lower=True)
    del Ka

    # --- dominance -------------------------------------------------------
    Kd = build_K(Zd)
    Ld = cholesky(s2d * Kd + stability * np.eye(n), lower=True)
    del Kd

    return La, Ld, Lgxg, w_build_time, c


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


def simulate_phenotype(La, Ld, Lgxg, n, s2a=0.1, s2d=0.1, s2gxg=0.1, s2e=0.7,
                       return_realized=False, force_realized=True):
    """Draw one phenotype  y = g_a + g_d + g_gxg + e,

        g_a = P_c La u1 ,        g_d = P_c Ld u2 ,
        g_gxg = P_c Lgxg u3 ,    e = P_c sqrt(s2e) u4 ,

    the four draws independent and CENTERED.  Centering a draw from W is, in
    distribution, a draw from P_c W P_c -- the kernel the estimator applies --
    so the Cholesky factors are unchanged.  Returns y, or
    (y, V_ell, V_a, V_d, V_e) with return_realized=True; all Var-hat at ddof=0.

    Because Lgxg Lgxg' = (s2gxg/(P c)) H H', g_gxg has exactly the law of
    H gamma with gamma ~ N(0, (s2gxg/(P c)) I_P) without ever forming H; the
    same argument gives g_a = Z_a beta and g_d = Z_d delta.

    Realized variances are RANDOM.  For the raw kernel E[Var-hat(g_gxg)] =
    c * s2gxg; here W carries 1/c-hat, so it is (c/c-hat) s2gxg -- i.e. s2gxg to
    the plug-in's accuracy.  The other three need no correction.

    force_realized=True (THE DEFAULT) rescales each component by
    sqrt(target / Var-hat) so its realized variance hits the target exactly.
    Three consequences:

      - an estimate's deviation from the target becomes the ESTIMATOR's error
        alone, and the residual c/c-hat factor is absorbed;
      - IT CHANGES THE GENERATING LAW.  Dividing by a random sample variance
        leaves y non-Gaussian, so REML fits a model the simulation does not
        obey: a MISSPECIFICATION, not a variance reduction, though small at
        these n;
      - the realized-variance columns become constants, so calc_stats.py's
        paired diagnostic collapses onto the unpaired one.

    A component with target 0 is left at zero.  force_realized=False restores
    the plain draw, which is what the sibling pipelines do.

    NOTE the four do NOT sum to var(y) under either setting: the sample
    cross-products are not zero, and rescaling touches only the diagonal terms.
    """
    u1 = np.random.randn(n)
    u2 = np.random.randn(n)
    u3 = np.random.randn(n)
    u4 = np.random.randn(n)

    # CENTERED draws.  Centering a draw from W is, in distribution, a draw from
    # P_c W P_c, so the Cholesky factors are unchanged and only the projection
    # is new.  The estimator works in the same contrast space.
    a = _center_cols(La @ u1)             # additive effect   ~ N(0, s2a K_a)
    d = _center_cols(Ld @ u2)             # dominance effect  ~ N(0, s2d K_d)
    gxg = _center_cols(Lgxg @ u3)         # epistasis  ~ N(0, s2gxg P_c W P_c)
    e = _center_cols(np.sqrt(s2e) * u4)   # residual noise

    if force_realized:
        # ddof=0, matching how the realized columns report it.  See the
        # docstring: this changes the generating law.
        a = _force_var(a, s2a, 'a')
        d = _force_var(d, s2d, 'd')
        gxg = _force_var(gxg, s2gxg, 'gxg')
        e = _force_var(e, s2e, 'e')

    y = a + d + gxg + e
    if return_realized:
        # realized Var-hat of each drawn component, ddof=0
        return (y, float(gxg.var()), float(a.var()), float(d.var()),
                float(e.var()))
    return y


# ------------------------------------ realized-variance scale factor  c
def compute_c_pooled(real_data, G):
    """c-hat from raw dosages: builds Z_a, cuts it into G genes, returns pooled_c.

    NOT where the kernels get their divisor -- they call pooled_c directly.
    Here the kernel already carries 1/c-hat, so no post-fit correction is
    applied.  K_a and K_d need no analogue: c_a = c_d = 1 exactly, their
    columns being standardized.
    """
    Z = additive_design(real_data)
    genes = split_into_genes(Z, G)
    return pooled_c(genes)


# ------------------------------------------------ Linear operator for the 4-VC system
def _v_matvec(Z, Zd, wu, s2a, s2d, s2gxg, s2e, B):
    """Apply  V = s2a K_a + s2d K_d + s2gxg W-hat + s2e I  to B, matrix-free.

    B is (n,) or (n, k) and the result matches.  No n-by-n GRM exists: K_a and
    K_d are a gemm pair each, the epistasis term is wu (the un-wrapped W-hat
    apply B -> W-hat B, compute_WU_storedA) wrapped in P_c on both sides (K_a and K_d need no
    wrapping -- K 1 = 0 already).  THE UNIT OF WORK: CG sees V and nothing else;
    a centered right-hand side stays centered, since every term preserves it.

    """
    _count_apply('V', 1 if np.ndim(B) == 1 else np.shape(B)[1])
    t0 = time.perf_counter()
    Bc = _center_cols(B)
    WB = wu(Bc)
    out = (s2a * compute_KU(Z, B)
           + s2d * compute_KU(Zd, B)
           + s2gxg * _center_cols(WB)
           + s2e * B)
    _add_time('V', time.perf_counter() - t0)    # includes its K and W applies
    return out


def _cg_batched(matvec, Bmat, x0=None, tol=1e-6, maxiter=1000):
    """Conjugate gradient for the SPD system V X = Bmat.

    matvec : X -> V X, (n, c) in and out.
    Bmat   : (n, c) right-hand sides, solved together with per-column scalars,
             so one V-pass advances every column.
    x0     : (n, c) warm start.

    Returns (X, ok).  ok is False when a column meets a search direction with
    p'Vp <= 0: V is not positive definite here (W-hat is a truncation, so V can
    be indefinite even with every component >= 0), CG is not a valid solver,
    and X is NOT a solution.  Counted as cg_negcurv.  Costs 1 + (iterations
    taken) applies of matvec.
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
                 and so a deterministic function of (g, A), whichever W
                 route produced them;
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


# ----------------------------------------------------------- MC AI-REML (k=4)
def mc_reml(Z, Zd, y, G, iters=30, Nmc=100, Nmc_coarse=15, cg_tol=1e-6,
            cg_maxiter=1000, jitter=1e-10, tol_ll=1e-4, tol_ll_coarse=1e-2,
            eta1=1e-4, eta2=0.99, alpha1=0.25, alpha2=3.5,
            s_init=(0.05, 0.05, 0.05, 0.85),
            seed=None, verbose=False, r=R_DEFAULT, w_route='storedA',
            A_dtype='float64', max_A_gb=16.0):
    """Monte-Carlo AI-REML for  V = s2a K_a + s2d K_d + s2gxg W + s2e I, with
    W = W_raw/c applied by its rank-r truncation, optimized as in BOLT-REML
    (Loh et al. 2015, Supplementary Note 3.2-3.5).

    The fit runs in the CONTRAST SPACE: y and the probes are centered and the
    W apply is wrapped in P_c, so the epistasis kernel is P_c W P_c.  That makes
    this proper REML with an intercept and keeps lam_max(W) off the all-ones
    direction, where most of it otherwise sits.

    Z is the standardized genotype (K_a, and split into G genes for the
    epistasis setup); Zd is the dominance design (K_d only).  Both are n-by-m
    and are the only large arrays held.  Because the kernel carries 1/c-hat,
    s2gxg is the REALIZED epistasis variance and needs no post-fit correction.
    K_a and K_d are EXACT; only W is truncated, so any error in s2a-hat or
    s2d-hat is sampling or coupling, never operator error.  W-hat is
    deterministic in (genotype, r) and its error is a truncation BIAS -- small
    under LD, large under linkage equilibrium, reduced only by raising r.

    W apply route (w_route)
    -----------------------
    Every W apply in the fit -- the P_c-wrapped Wapply, the V applies inside
    CG, the fixed-probe W U -- goes through ONE apply chosen here, so the two
    routes differ only in how W-hat U is computed; the arithmetic is the same
    truncation and agrees to round-off.

      'storedA'  (default) builds (At, D) once via setup_storedA, in A_dtype
                 and refused over max_A_gb, and applies compute_WU_storedA: two
                 gemm pairs per apply.  Fastest on wide applies (the probe
                 solves), but At costs n * r * m_W memory (m_W = SNPs in genes)
                 and the gemms stream all of it, so it is sensitive to memory
                 bandwidth.
      'bcast'    keeps setup_pooled's per-gene factors and applies
                 compute_WU_bcast: per gene, q_s .* U is broadcast and pushed
                 through Z_g.  No extra memory beyond the factors, and more
                 robust on a crowded node.  A_dtype / max_A_gb are ignored.

    Parameter domain (Note 3.2.1)
    -----------------------------
    s2a, s2d, s2gxg >= 0 and s2e >= 1e-9 var(y); no upper bound.  s_init is a
    fraction of var(y) and must satisfy the bounds.

    Score traces
    ------------
    tr(V^{-1} K_i) by HUTCHINSON: Nmc fixed centered Rademacher probes (BOLT's
    simulated y_V needs a PSD factor of every kernel, which W-hat does not
    have), K_a U / K_d U / W U formed once, Nmc warm-started CG solves per
    iteration.  Fixed probes make the objective a deterministic function of s.

    TWO-PHASE SCHEDULE (Note 3.5): Nmc_coarse probes to tol_ll_coarse, then
    all Nmc to tol_ll.  The probes are drawn ONCE at full width and the coarse
    phase uses the first Nmc_coarse columns, so the switch EXTENDS the objective
    rather than replacing it -- the products are formed once, the fixed point
    moves continuously, and Delta carries over.  The step that triggers the
    switch is discarded and re-evaluated at full Nmc.  Nmc_coarse=None (or
    >= Nmc) gives the single-phase fit.

    Constrained trust-region step (Note 3.2.1, 3.4)
    -----------------------------------------------
    Each iteration solves the 4-dimensional QP

        maximize  g'p - 0.5 p'AI p
        s.t.      s + p >= lower bounds ,  ||diag(AI) p|| <= Delta

    (_trust_region_step: SLSQP, the norm constraint dropped while Delta is
    Inf).  A component whose likelihood peaks at 0 lands ON its bound and the
    others keep converging in the free subspace.  dLL_pred, the QP objective at
    the solution, is the predicted gain; the actual gain would need log|V|,
    which this construction never forms, so it is approximated by the
    TRAPEZOID rule 0.5 p'(g(s) + g(s+p)) -- the likelihood being the line
    integral of the score.  rho = actual/predicted drives the accept test and
    the radius update; a gradient that more than doubles means the quadratic
    model has broken down and sets rho = -1.  Delta starts at Inf.  An accepted
    step's trial score is reused, so only REJECTED steps cost an extra score
    evaluation.  Converged when dLL_pred < tol_ll in the fine phase.

    The bounds keep every component >= 0, but W-hat is a truncation and not
    exactly PSD, so V can still be indefinite at a trial point.  CG detects
    that (p'Vp <= 0, counted as cg_negcurv) and the trial point is rejected
    exactly like a rho = -1 step.

    Returns
    -------
    s    : (4,) estimated (s2a, s2d, s2gxg, s2e).
    AI   : (4, 4) average information at the FULL Nmc.  SEs via
           reml_se(AI, Nmc, info['at_bound']); a raw sqrt(diag(inv(AI)))
           understates them.
    info : dict -- converged (dLL_pred < tol_ll reached in the fine phase),
           n_iters, n_reject, switch_it, at_bound (bool per component).
    """
    reset_op_counts()                  # counts describe THIS replicate alone
    if w_route not in W_ROUTES:
        raise ValueError(f"w_route must be one of {W_ROUTES}, got {w_route!r}")
    y = np.asarray(y, dtype=float).flatten()
    y = y - y.mean()                   # work in the contrast space throughout
    Z = np.asarray(Z, dtype=float)
    Zd = np.asarray(Zd, dtype=float)
    n = y.shape[0]
    k = 4

    # --- setup: genotype-only, done once ---
    # Timed in two pieces (svd factors incl. c-hat, the fixed-probe products
    # K_i U / W U) plus their total, setup_sec; the At build (stored-A route
    # only) is setup_A_sec.
    t_setup0 = time.perf_counter()
    genes = split_into_genes(Z, G)
    F_list, P, c = setup_pooled(genes, r=r)      # c: the normalization divisor
    _add_time('setup_svd', time.perf_counter() - t_setup0)
    # The ONE un-wrapped W-hat apply every W apply below goes through.
    if w_route == 'storedA':
        At, Dst = setup_storedA(F_list, dtype=A_dtype, max_gb=max_A_gb)
        # At/Dst hold everything the apply reads, so drop the per-gene copies.
        # split_into_genes COPIES (fancy indexing), so this frees the Z_g copies.
        F_list = genes = None
        wu = lambda B: compute_WU_storedA(At, Dst, P, c, B)
    else:
        # F_list holds each gene's Z_g, D_g, Q and lam -- all the apply reads.
        genes = None
        wu = lambda B: compute_WU_bcast(F_list, P, c, B)
    Kaapply = lambda B: compute_KU(Z, B)
    Kdapply = lambda B: compute_KU(Zd, B)
    Wapply = lambda B: _center_cols(wu(_center_cols(B)))

    vary = y.var()
    lower = np.array([0.0, 0.0, 0.0, 1e-9 * vary])   # the domain, Note 3.2.1

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
    xbuf = None       # warm-start buffers for the CG solve groups
    # Held at the FULL Nmc so the fine phase inherits the coarse phase's warm
    # starts; the untouched columns are cold-started once, at the switch.
    Pbuf = np.zeros((n, Nmc))
    Gbuf = None

    t0 = time.perf_counter()
    KaU = Kaapply(U)  # K_a U for the fixed probes: genotype-only, formed once
    KdU = Kdapply(U)  # K_d U   "        "
    WU = Wapply(U)    # W U     "        "
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
        s2a_, s2d_, s2gxg_, s2e_ = s_at
        mv = lambda B: _v_matvec(Z, Zd, wu, s2a_, s2d_, s2gxg_, s2e_, B)

        # --- x = V^{-1} y ---
        t0 = time.perf_counter()
        X_y, ok = _cg_batched(mv, yc, x0=xbuf, tol=cg_tol, maxiter=cg_maxiter)
        _add_time('solve_y', time.perf_counter() - t0)
        if not ok:
            return None, None, False
        x_ = X_y[:, 0]

        # data quadratics x'K_i x
        Kax_, Kdx_, Wx_ = Kaapply(x_), Kdapply(x_), Wapply(x_)

        # --- score traces tr(V^{-1} K_i): Hutchinson, exact probe solves ---
        Uc_ = U[:, :n_probe]
        t0 = time.perf_counter()
        Psol, ok = _cg_batched(mv, Uc_, x0=Pbuf[:, :n_probe], tol=cg_tol,
                               maxiter=cg_maxiter)
        _add_time('solve_probe_coarse' if coarse else 'solve_probe_fine',
                  time.perf_counter() - t0)
        if not ok:
            return None, None, False
        sc = np.array([
            0.5 * (x_ @ Kax_ - np.mean(np.sum(Psol * KaU[:, :n_probe], axis=0))),
            0.5 * (x_ @ Kdx_ - np.mean(np.sum(Psol * KdU[:, :n_probe], axis=0))),
            0.5 * (x_ @ Wx_ - np.mean(np.sum(Psol * WU[:, :n_probe], axis=0))),
            0.5 * (x_ @ x_ - np.mean(np.sum(Psol * Uc_, axis=0)))])

        # --- average information: A_ij = 0.5 (K_i x)' V^{-1}(K_j x) ---
        KX = np.column_stack([Kax_, Kdx_, Wx_, x_])
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


def MC_REML(Z, Zd, y, G, iters=30, Nmc=100, Nmc_coarse=15, cg_tol=1e-6,
            cg_maxiter=1000, seed=None, r=R_DEFAULT, verbose=False,
            w_route='storedA', A_dtype='float64', max_A_gb=16.0):
    """Wrapper: returns (s2a_hat, s2d_hat, s2gxg_hat, s2e_hat, AI, info).

    SEs from reml_se(AI, Nmc, info['at_bound']).  Nmc_coarse=None gives the
    single-phase fit.  w_route picks the W apply ('storedA' or 'bcast');
    A_dtype and max_A_gb set and bound the stored-A one (see mc_reml).  info
    is mc_reml's: converged, n_iters, n_reject, switch_it, at_bound.
    """
    s, AI, info = mc_reml(Z, Zd, y, G, iters=iters, Nmc=Nmc,
                          Nmc_coarse=Nmc_coarse, cg_tol=cg_tol,
                          cg_maxiter=cg_maxiter, seed=seed, verbose=verbose,
                          r=r, w_route=w_route, A_dtype=A_dtype,
                          max_A_gb=max_A_gb)
    s2a_hat, s2d_hat, s2gxg_hat, s2e_hat = s
    return s2a_hat, s2d_hat, s2gxg_hat, s2e_hat, AI, info
