# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.linalg import cholesky
from scipy.sparse.linalg import svds
import time

####################################################################
# FOUR variance components: additive + dominance + pooled within-gene
# epistasis + noise.
#
# Model (no fixed effects):
#     y = g_a + g_d + g_gxg + e ,
#     V = Var(y) = s2a K_a + s2d K_d + s2gxg W + s2e I ,
#     K_a = Z_a Z_a' / m                  (additive GRM, standardized dosages),
#     K_d = Z_d Z_d' / m                  (dominance GRM, standardized GCTA
#                                          dominance coding),
#     W   = (1/P) sum_g sum_{a<b in g} h_ab h_ab' ,  h_ab = Z_a .* Z_b ,
#           P = sum_g C(m_g, 2)           (pooled within-gene epistasis GRM,
#                                          UNSTANDARDIZED interactions).
#
# This is the two-component Lowrank_Wu_unstd pipeline with the ADDITIVE and
# DOMINANCE components added back, and nothing else changed: the epistasis
# kernel, its rank-r low-rank apply, the realized-variance factor c, the
# Hutchinson trace estimator and the damped AI-Newton optimizer are all carried
# over verbatim, so a run here differs from that sibling by the extra K_a and
# K_d ALONE.
#
# Defaults are s2a = 0.1, s2d = 0.1, s2gxg = 0.1, s2e = 0.7.
#
# THE EPISTASIS KERNEL IS BUILT FROM Z_a ONLY.  h_ab = Z_a .* Z_b pairs
# ADDITIVE columns; dominance enters through K_d and through nothing else.
# Dominance-by-dominance and additive-by-dominance interactions are NOT in this
# model -- they would be new kernels, not a new component on the existing W.
#
# COSTS.  K_a and K_d each enter as the thin mat-vec  K_i B = Z_i(Z_i'B)/m,
# O(n m c) -- strictly cheaper than the O(n m r c) epistasis apply they sit
# next to, so the extra components add no new order of cost to the estimator.
# Only the SIMULATION pays for them in n-by-n terms: two more dense Cholesky
# factors.  The estimator now holds TWO n-by-m designs instead of one.
#
# IDENTIFIABILITY.  K_a, K_d and W are far less collinear with each other than
# W is with I, but the components are not free: the AI matrix now has a
# 4-dimensional curvature to resolve, and small n / large m runs that were
# merely noisy in the two-component fit can become weakly identified here.
# K_d is the weakest of the three in practice -- dominance deviations are small
# and, at low MAF, their GRM is close to I -- so it is the component most
# likely to sit on the boundary.  The optimizer's Levenberg-Marquardt ridge +
# trust region + box clamp are unchanged and carry the same known boundary
# defect (see mc_reml).
####################################################################


# ------------------------------------------------------------------ designs
def _standardize_cols(M, stability_std=1e-12):
    """Column-standardize to mean 0, variance 1 (ddof=0).

    A near-constant column (std < stability_std) is left un-scaled to avoid a
    divide-by-zero; such columns contribute ~0 to the GRM anyway.
    """
    M = np.asarray(M, dtype=float)
    mu = M.mean(axis=0)
    sd = M.std(axis=0)
    sd = np.where(sd < stability_std, 1.0, sd)
    return (M - mu) / sd


def additive_design(real_data):
    """Additive design Z_a: column-standardized allele dosages.

    The SAME Z_a feeds two of the three genetic kernels: K_a = Z_a Z_a'/m is
    its GRM, and the epistasis interactions h_ab = Z_a .* Z_b are built from
    its columns.  Only K_d has a design of its own.
    """
    return _standardize_cols(real_data)


def dominance_design(real_data, maf_floor=1e-6):
    """Dominance design Z_d: column-standardized GCTA dominance coding.

    Genotypes are allele counts {0,1,2}; with p = column mean / 2 and q = 1 - p,
    the orthogonal dominance deviation maps

        0 -> -p/q ,   1 -> 1 ,   2 -> -q/p ,

    which has mean 0 under HWE, so the coding is orthogonal to the additive one
    in expectation -- that is what makes s2a and s2d separately estimable
    rather than two names for the same direction.  The columns are then
    standardized like the additive ones, so tr(K_d) = n exactly and c_d = 1
    (see compute_c_pooled).

    p is clipped to [maf_floor, 1 - maf_floor] so a monomorphic column cannot
    divide by zero; such a column is constant, hence ~0 after standardization,
    and contributes nothing to K_d.  Values are rounded to the nearest integer
    first: the coding is defined on hard calls, so imputed dosages are snapped
    rather than silently mis-coded.

    NOTE the ORTHOGONALITY IS ONLY IN EXPECTATION, and only under HWE.  In a
    finite sample Z_a'Z_d is not exactly 0, so K_a and K_d overlap a little;
    that overlap is real and is exactly what the AI matrix has to resolve.
    """
    X = np.rint(np.asarray(real_data, dtype=float))
    p = X.mean(axis=0) / 2.0
    p = np.clip(p, maf_floor, 1.0 - maf_floor)
    q = 1.0 - p
    D = (X == 0) * (-p / q) + (X == 1) * 1.0 + (X == 2) * (-q / p)
    return _standardize_cols(D)


def split_into_genes(Z, G):
    """Split the m SNP columns of Z into G contiguous gene blocks.

    Returns a list [Z_1, ..., Z_G] of column sub-matrices -- the SEVERAL Z, one
    per gene, that the pooled kernel consumes.  np.array_split handles m not
    divisible by G by making the first (m mod G) genes one SNP larger, so the
    genes may be unequal.

    The ADDITIVE kernel does not use this split at all: K_a is built from the
    whole Z, pooling every SNP with weight 1/m regardless of gene membership.
    """
    m = Z.shape[1]
    return [Z[:, cols] for cols in np.array_split(np.arange(m), G)]


# --------------------------------------- standardized GRMs (matrix-free)
# K_a and K_d have IDENTICAL structure -- Z_i Z_i' / m for a column-standardized
# design Z_i -- and differ only in which design is passed in.  One routine
# serves both; nothing here knows or cares whether it is holding dosages or
# dominance deviations.
def compute_KU(Zi, U):
    """Matrix-free product  K_i @ U  with  K_i = Z_i Z_i' / m.

    Zi is a column-standardized design (additive_design or dominance_design).
    U may be (n,) or (n, c); the result matches its shape.  O(n m c) time and
    O(n c) extra storage -- the n-by-n K_i is never formed, exactly as W is
    never formed.  This is the only route the ESTIMATOR uses; the simulation
    forms the GRMs densely because a Cholesky factor needs an explicit matrix.

    tr(K_i) = n EXACTLY for either design (every column has sample variance 1,
    so sum_t sum_a Z_ta^2 / m = n) -- the identity the unstandardized epistasis
    kernel gives up but the standardized ones keep.
    """
    Zi = np.asarray(Zi, dtype=float)
    U = np.asarray(U, dtype=float)
    m = Zi.shape[1]
    single = (U.ndim == 1)
    if single:
        U = U.reshape(-1, 1)
    out = (Zi @ (Zi.T @ U)) / m
    return out[:, 0] if single else out


def build_K(Zi):
    """Explicit GRM  K_i = Z_i Z_i' / m.  SIMULATION only (n-by-n)."""
    Zi = np.asarray(Zi, dtype=float)
    return (Zi @ Zi.T) / Zi.shape[1]


# ------------------------------------- epistasis GRM (dense, SIMULATION only)
def _within_gene_sum(Zg):
    """UN-normalized within-gene epistasis GRM  S = sum_{a<b} h_ab h_ab',
    with the UNSTANDARDIZED interaction column  h_ab = Z_a .* Z_b.

    There is NO centering and NO 1/sigma_ab scaling -- that is the whole
    difference from the standardized pipelines -- and no 1/p_g division either:
    the Pooled Model divides ONCE by the global pair total P (see
    build_W_pooled), giving every pair equal weight.  Returns (S, p_g).

    Uses the closed form  2 sum_{a<b} h_ab h_ab' = (K .* K) - D D'  with
    K = Z_g Z_g' and D = Z_g .* Z_g:

        S = 0.5 [ (K .* K) - D D' ] ,     O(n^2 m_g) time, O(n^2) storage.

    This is the SAME Hadamard square the low-rank operator truncates, so
    simulation and estimation differ ONLY by the truncation -- there is no
    other gap between them.  The identity is exact: it agrees with the literal
    C(m_g, 2)-term pair sum to machine precision (see verify_lowrank.py).
    """
    n, mg = Zg.shape
    pg = mg * (mg - 1) // 2
    K = Zg @ Zg.T
    D = Zg * Zg
    S = 0.5 * (K * K - D @ D.T)
    return S, pg


def build_W_pooled(Z_list):
    """Pooled WITHIN-gene pairwise-epistasis GRM, UNSTANDARDIZED interactions.

        W = (1/P) sum_{g=1}^G sum_{a<b in g} h_ab h_ab' ,   h_ab = Z_a .* Z_b ,
        P = sum_{g=1}^G C(m_g, 2) .

    A 1-SNP gene contributes no pair and is skipped.  n-by-n storage; O(n^2 m)
    time.  USED ONLY BY THE SIMULATION (the Cholesky factor needs an explicit
    matrix); estimation applies a rank-r TRUNCATION of the same W matrix-free
    via compute_WU_pooled.

    tr(W) is NOT n here.  With standardized SNPs, W_tt = (1/2P) sum_g
    [(sum_a Z_ta^2)^2 - sum_a Z_ta^4], so tr(W) = n only in expectation under
    independence and fluctuates with the genotype -- the standardized siblings'
    tr(W) = n identity does not survive dropping the 1/sigma_ab weight.  K_a,
    which IS standardized, keeps tr(K_a) = n exactly, so the two genetic
    components are not on the same trace scale.
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
    W /= P
    return W


# --------------------------------------- matrix-free low-rank pooled W apply
R_DEFAULT = 20          # truncation level r (the note's default; raise until
                        # the estimates stop moving -- LD-dependent, see header)


def setup_pooled(Z_list, r=R_DEFAULT):
    """SETUP phase: everything that depends on the genotype alone.

    Returns (F_list, P):

        F_list : per-gene state, one dict per gene (None for a 1-SNP gene)
        P      : sum_g C(m_g, 2)  -- the global pair total

    The leading eigenpairs of K_g = Z_g Z_g' come from a truncated SVD of Z_g;
    K_g is never formed.  ARPACK needs k < min(n, m_g) and returns singular
    values in ASCENDING order, so a gene small enough that r reaches full rank
    is factorized exactly instead -- which is also where the operator becomes
    EXACT.
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
    return F_list, P


def _gene_WU(gene, U, buf_elems):
    """ONE gene's un-normalized contribution,  2 S_g U, by the rank-r truncation.

        (K .* K) u ~ sum_{s=1}^r lam_s q_s .* (Z(Z'(q_s .* u))) ,

    """
    Zg, Dg, Q, lam = gene['Z'], gene['D'], gene['Q'], gene['lam']
    n, c = U.shape
    r = Q.shape[1]
    out = np.empty((n, c))
    Ql = Q * lam                                 # fold lam into the left factor

    cb = max(1, min(c, buf_elems // max(1, n * r)))
    for s in range(0, c, cb):
        e = min(s + cb, c)
        w = e - s
        Ub = U[:, s:e]

        Tb = (Q[:, :, None] * Ub[:, None, :]).reshape(n, r * w)   # q_s .* u
        Ob = (Zg @ (Zg.T @ Tb)).reshape(n, r, w)                  # K(q_s .* u)
        out[:, s:e] = np.einsum('ns,nsw->nw', Ql, Ob)             # lam q_s .* (.)

    out -= Dg @ (Dg.T @ U)                   # the a = b terms the pair sum drops
    return out


def compute_WU_pooled(Z_list, F_list, P, U, buf_elems=8_000_000):
    """Matrix-free product  W-hat @ U  for the pooled UNSTANDARDIZED epistasis
    GRM, W-hat being the rank-r truncation of W.

    U may be (n,) or (n, c); the result matches its shape.  Never forms W
    (n-by-n), any interaction block H_g (n-by-p_g), or anything m-by-m.
    Implements

        W-hat U = 1/(2P) sum_g 2 S-hat_g U ,

    """
    n = Z_list[0].shape[0]
    U = np.asarray(U, dtype=float)
    single = (U.ndim == 1)
    if single:
        U = U.reshape(n, 1)

    T1 = np.zeros_like(U)
    for gene in F_list:
        if gene is None:                     # 1-SNP gene: no pair, no work
            continue
        T1 += _gene_WU(gene, U, buf_elems)

    out = T1 / (2.0 * P)
    return out[:, 0] if single else out


# ------------------------------------------------------------- simulation
def simulate_Cholesky_4vc(real_data, G, s2a=0.1, s2d=0.1, s2gxg=0.1, s2e=0.7,
                          stability=1e-10):
    """Cholesky factors of the additive, dominance and pooled-epistasis
    covariances.

        La La'     = s2a   K_a ,   K_a = Z_a Z_a' / m ,
        Ld Ld'     = s2d   K_d ,   K_d = Z_d Z_d' / m ,
        Lgxg Lgxg' = s2gxg W .

    Returns (La, Ld, Lgxg, w_build_time), w_build_time being the wall-clock
    cost of the EPISTASIS kernel alone -- the O(n^2 m) term that dominates this
    job and the one the sibling pipelines report.  K_a and K_d cost O(n^2 m)
    too but with a single gemm each and no Hadamard square, so they are folded
    into the untimed remainder rather than given their own records.

    The three factors are built SEQUENTIALLY and each dense GRM is released as
    soon as its factor exists, so the peak is ~4 n-by-n arrays (two finished
    factors + one GRM + the factorization workspace) rather than 7.  At n =
    16000 that is ~8 GB instead of ~14 GB -- still the memory-critical job of
    the pipeline, and the reason it gets its own larger SLURM allocation.
    """
    Za = additive_design(real_data)
    Zd = dominance_design(real_data)
    n, m = Za.shape

    # --- epistasis first (the expensive one), then free W -------------------
    genes = split_into_genes(Za, G)          # several Z, one per gene

    t_start = time.perf_counter()
    W = build_W_pooled(genes)
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

    return La, Ld, Lgxg, w_build_time


def simulate_phenotype(La, Ld, Lgxg, n, s2a=0.1, s2d=0.1, s2gxg=0.1, s2e=0.7,
                       return_realized=False):
    """Draw one phenotype  y = g_a + g_d + g_gxg + e  from the four-component
    model,

        g_a   = La u1   ~ N(0, s2a K_a) ,   K_a = Z_a Z_a' / m ,
        g_d   = Ld u2   ~ N(0, s2d K_d) ,   K_d = Z_d Z_d' / m ,
        g_gxg = Lgxg u3 ~ N(0, s2gxg W) ,   W   = (1/P) H H' ,
        e     = sqrt(s2e) u4 ~ N(0, s2e I) .

    The four draws are INDEPENDENT, so the phenotype has exactly the law the
    estimator fits -- no rescaling to hit the target variances, matching the
    two-component pipeline (sampling-error removal is deliberately not done
    here: it would make the realized variances below degenerate).

    Since Lgxg Lgxg' = s2gxg W = (s2gxg / P) H H', the draw g_gxg has EXACTLY
    the law of  H gamma  with  gamma ~ N(0, (s2gxg / P) I_P) -- the effect-size
    model of realized_variance.pdf (V_gamma = s2gxg) -- without ever forming
    the n-by-P interaction design H.  The same argument makes g_a exactly
    Z_a beta with beta ~ N(0, (s2a/m) I_m), and g_d exactly Z_d delta with
    delta ~ N(0, (s2d/m) I_m).

    REALIZED VARIANCE.  s2gxg (V_gamma in the note) is the per-pair effect
    variance scaled by P; it is NOT the variance of the genetic value across
    individuals.  That quantity,

        V_ell := Var-hat(H gamma) = (1/n) g_gxg' P_c g_gxg ,   P_c = I - 11'/n ,

    is random (it changes with every gamma draw) and its expectation is
    E[V_ell] = c * s2gxg, c = (1/P) sum_{pairs} Var-hat(Z_a .* Z_b) -- see
    compute_c_pooled.  With return_realized=True the function also returns
    V_ell for THIS draw (ddof=0, matching the note's 1/n normalisation), so
    the simulation can record the target the estimator's c-corrected output
    c-hat * s2gxg-hat is meant to recover.

    The ADDITIVE and DOMINANCE components need no such correction: their
    columns ARE standardized, so their per-column variance factors are
    c_a = c_d = 1 exactly and E[Var-hat(Z_i eff)] = s2i with no rescaling (see
    compute_c_pooled's note).  Their realized variances are still RANDOM -- one
    effect draw is not its expectation -- so they are reported too, and so is
    the residual's:

        V_a := Var-hat(Z_a beta)  = var(a) ,   E[V_a] = s2a   (c_a = 1) ,
        V_d := Var-hat(Z_d delta) = var(d) ,   E[V_d] = s2d   (c_d = 1) ,
        V_e := Var-hat(e)         = var(e) ,   E[V_e] = s2e   (trivially) ,

    all ddof=0 like V_ell.  These are the draw-level targets s2a-hat, s2d-hat
    and s2e-hat should be compared against replicate by replicate; against the
    NOMINAL s2a, s2d and s2e they carry an extra O(1/sqrt n) scatter that is
    the draw's, not the estimator's.

    NOTE the four realized variances do NOT sum to var(y).  a, d, gxg and e are
    independent in expectation but their SAMPLE cross-products are not zero, so
    var(y) = V_a + V_d + V_ell + V_e + 2(cov-hat terms); at n = 16000 the cross
    terms are O(1/sqrt n) but they are not nothing.  Do not treat the four as a
    partition of the phenotypic variance.

    Returns y, or (y, V_ell, V_a, V_d, V_e) when return_realized is True.
    """
    u1 = np.random.randn(n)
    u2 = np.random.randn(n)
    u3 = np.random.randn(n)
    u4 = np.random.randn(n)

    a = La @ u1                           # additive effect   ~ N(0, s2a K_a)
    d = Ld @ u2                           # dominance effect  ~ N(0, s2d K_d)
    gxg = Lgxg @ u3                       # epistasis effect  ~ N(0, s2gxg W)
    e = np.sqrt(s2e) * u4                 # residual noise

    y = a + d + gxg + e
    if return_realized:
        # realized Var-hat of each drawn component, ddof=0
        return (y, float(gxg.var()), float(a.var()), float(d.var()),
                float(e.var()))
    return y


# ------------------------------------ realized-variance scale factor  c
def hwe_skewness(real_data):
    """Skewness of every column-standardized SNP under HWE, from allele
    frequency alone:

        s_i = (1 - 2 p_i) / sqrt(2 p_i (1 - p_i)) ,   p_i = mean(X_i) / 2 ,

    X being the 0/1/2 dosage matrix.  O(nm).

    A monomorphic SNP (p_i in {0, 1}) has no defined s_i and gets 0.  Such a
    column is identically 0 after standardization, so every pair it enters has
    Var-hat = 0 exactly, whereas the closed form in compute_c_pooled then
    counts that pair as 1 + 0 = 1: an O(#monomorphic / m) relative error in
    c.  Use MAF-filtered genotypes, as the rest of the pipeline assumes.
    """
    X = np.asarray(real_data, dtype=float)
    p = X.mean(axis=0) / 2.0
    denom = 2.0 * p * (1.0 - p)
    ok = denom > 0.0
    s = np.zeros(p.shape[0])
    s[ok] = (1.0 - 2.0 * p[ok]) / np.sqrt(denom[ok])
    return s


def compute_c_pooled(real_data, G, method='hwe'):
    """Scale factor c between the variance COMPONENT and the REALIZED variance
    of the EPISTASIS genetic value (realized_variance.pdf):

        E[ Var-hat(H gamma) ] = c * s2gxg ,
        c = (1/P) sum_g sum_{a<b in g} Var-hat(Z_a .* Z_b) ,  P = sum_g C(m_g, 2) ,

    for the pooled WITHIN-gene unstandardized kernel: only within-gene pairs
    make up H, exactly the pairs build_W_pooled / setup_pooled count in P.
    The REML fit of  V = s2a K_a + s2d K_d + s2gxg W + s2e I  with the raw
    (uncentered, unscaled) H estimates V_gamma = s2gxg; the estimate of the
    realized variance is then  V_ell-hat = c * s2gxg-hat.  Because H is not
    column-standardized, c != 1 in general -- it is exactly the
    per-pair-variance factor the standardized siblings divide out inside their
    kernel.

    THE ADDITIVE AND DOMINANCE COMPONENTS NEED NO ANALOGUE.  Their designs are
    Z_a and Z_d, whose columns are standardized to sample variance 1 (ddof=0)
    by construction, so c_a = (1/m) sum_i Var-hat(Z_a,i) = 1 EXACTLY, and
    likewise c_d = 1 -- not approximately, and with no HWE assumption (the HWE
    assumption in the dominance CODING is a separate matter: it decides what
    the deviation means, not how its column is scaled).  s2a-hat and s2d-hat
    are therefore already on the realized-variance scale and are carried
    through unmodified; there is no c_a or c_d to compute, to write out, or to
    get wrong.

    Takes the RAW 0/1/2 dosage matrix (the HWE route needs the allele
    frequencies, which the standardized Z no longer carries) and splits it
    into G contiguous genes the same way the kernel does.

    method='hwe'  (the note's estimator; DEFAULT, used by the pipeline)
        Under HWE the per-pair variance is  Var(Z_a Z_b) = 1 + r_ab s_a s_b,
        with s the HWE skewness (hwe_skewness) and r_ab the sample
        correlation.  R_g = Z_g' Z_g / n is symmetric with unit diagonal, so

            sum_{a<b} r_ab s_a s_b = 0.5 ( s' R_g s - ||s||^2 ) ,
            s' R_g s = ||Z_g s||^2 / n ,

        and  c-hat = 1 + (1 / 2P) sum_g ( ||Z_g s_g||^2 / n - ||s_g||^2 ) :
        one mat-vec per gene, O(nm) time, O(n) extra storage -- no R (m-by-m)
        and no H (n-by-P).  Its error against the exact c is the HWE
        approximation plus O(1/sqrt n) moment sampling, which the note's check
        finds ~30x smaller than the draw-to-draw scatter of Var-hat(H gamma)
        itself.

    method='exact'
        The mean sample variance of the P interaction columns, formed without
        H.  With D_g = Z_g .* Z_g and the per-gene Gram matrix Z_g' Z_g,

            sum_{a<b} mean(Z_a^2 Z_b^2) = (1 / 2n)  ( ||D_g 1||^2 - sum_{t,a} Z_ta^4 ) ,
            sum_{a<b} mean(Z_a Z_b)^2   = (1 / 2n^2)( ||Z_g' Z_g||_F^2 - sum_a (Z_a' Z_a)^2 ) ,

        O(n sum_g m_g^2) time, O(max m_g^2) storage.  This is the reference
        the closed form is measured against (the Cholesky job writes both);
        the note deliberately keeps it out of the estimation path.
    """
    if method not in ('hwe', 'exact'):
        raise ValueError(f"method must be 'hwe' or 'exact'; got {method!r}.")
    Z = additive_design(real_data)
    n, m = Z.shape
    genes = split_into_genes(Z, G)
    P = sum(Zg.shape[1] * (Zg.shape[1] - 1) // 2 for Zg in genes)
    if P == 0:
        raise ValueError("No gene block has >= 2 SNPs; increase m/G.")

    if method == 'hwe':
        s = hwe_skewness(real_data)
        s_genes = [s[cols] for cols in np.array_split(np.arange(m), G)]
        cross = 0.0
        for Zg, sg in zip(genes, s_genes):
            if Zg.shape[1] < 2:                  # no within-gene pair
                continue
            Zs = Zg @ sg                          # the ONE mat-vec per gene
            cross += (Zs @ Zs) / n - sg @ sg      # s' R_g s - ||s||^2
        return 1.0 + cross / (2.0 * P)

    total = 0.0
    for Zg in genes:
        if Zg.shape[1] < 2:
            continue
        Dg = Zg * Zg
        second = 0.5 * (np.sum(Dg.sum(axis=1) ** 2) - np.sum(Dg * Dg)) / n
        Gm = Zg.T @ Zg
        mean_sq = 0.5 * (np.sum(Gm * Gm) - np.sum(np.diag(Gm) ** 2)) / n ** 2
        total += second - mean_sq                 # sum_{a<b} Var-hat(Z_a Z_b)
    return total / P


# ------------------------------------------------ matrix-free linear algebra
def _v_matvec(Z, Zd, Z_list, F_list, P, s2a, s2d, s2gxg, s2e, B):
    """Apply V = s2a K_a + s2d K_d + s2gxg W + s2e I to B, every GRM
    matrix-free.

        V B = (s2a/m) Z_a(Z_a'B) + (s2d/m) Z_d(Z_d'B)
              + s2gxg (W-hat B) + s2e B .

    B may be (n,) or (n, c); the result matches its shape.  No n-by-n GRM
    exists: the additive and dominance terms are a pair of gemms each against
    their design, and the epistasis term is rebuilt from the cached low-rank
    factors on every call.
    """
    return (s2a * compute_KU(Z, B)
            + s2d * compute_KU(Zd, B)
            + s2gxg * compute_WU_pooled(Z_list, F_list, P, B)
            + s2e * B)


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


# ----------------------------------------------------------- MC AI-REML (k=4)
def mc_reml(Z, Zd, y, G, iters=30, Nmc=50, cg_tol=1e-6, cg_maxiter=1000,
            jitter=1e-8, tol=1e-8, lm=1e-3, step_frac=0.5, upper_mult=5.0,
            seed=None, verbose=False, r=R_DEFAULT):
    """Monte-Carlo AI-REML for V = s2a K_a + s2d K_d + s2gxg W + s2e I (pooled
    UNSTANDARDIZED W, applied by its rank-r truncation).

    Z is the column-standardized genotype (additive_design) and serves TWO of
    the three genetic kernels: K_a = Z Z'/m is applied straight from it, and it
    is split into G contiguous genes for the pooled epistasis setup (per-gene
    SVD factors, pair total P), built ONCE here outside the iteration.  Zd is
    the standardized dominance design (dominance_design), used by K_d and by
    nothing else -- in particular the interaction columns are additive-by-
    additive and never touch Zd.  Both designs are n-by-m and are the only
    large arrays the estimator holds.

    W u is applied by the note's O(n m r) low-rank operator and by nothing
    else: one copy of the Hadamard square (K .* K) is replaced by the rank-r
    eigendecomposition of K taken from the thin SVD of Z_g.  The operator is
    DETERMINISTIC -- W-hat is fixed by (genotype, r), identical across calls
    and replicates with no probe seed to manage -- and its error is a
    truncation BIAS set by the discarded eigenvalue tail: small under LD,
    large under linkage equilibrium, reduced only by raising r (see the module
    header).  K_a and K_d are EXACT: they are applied as written, with no
    truncation, so every discrepancy against the truth in s2a-hat or s2d-hat is
    sampling or coupling, never operator error.  verify_lowrank.py holds an
    independent exact implementation of the W apply to check this one against,
    including exactness at full rank.

    Score traces
    ------------
    tr(V^{-1} K_i) is estimated by HUTCHINSON and by nothing else: Nmc fixed
    Rademacher probes, K_a U, K_d U and W U each formed once up front, and Nmc
    CG solves for V^{-1}U per REML iteration (warm-started from the previous
    iteration).  The probes are fixed across iterations, so the objective is a
    deterministic function of s and the iteration converges to a fixed point.
    Each added component costs ONE extra probe product (formed once) and one
    extra right-hand side in the AI solve group -- the Nmc probe solves, which
    dominate, are unchanged.

    A KNOWN DEFECT OF THE SHARED OPTIMIZER, inherited deliberately.  When a
    replicate's likelihood peaks at a component = 0, that component sticks on
    the lower box clamp and the others then converge to the WRONG values: the
    AI-Newton step is solved jointly and never re-projected onto the free
    subspace, so the blocked direction keeps driving the free ones through the
    coupling terms.  With four components this is MORE likely to bite, not less
    -- s2a, s2d and s2gxg can each hit the boundary, and s2d is the likeliest
    of them, dominance being the weakest signal here -- but the optimizer block
    is byte-identical across the whole family and is left uncorrected so that a
    run here differs from a sibling run in the KERNELS alone; an active-set
    projection would fix it but would change exactly the replicates that make
    the runs comparable.  calc_stats.py counts the affected replicates.

    Returns
    -------
    s  : (4,) estimated (s2a, s2d, s2gxg, s2e).
    AI : (4, 4) final average-information matrix.
    """
    y = np.asarray(y, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float)
    Zd = np.asarray(Zd, dtype=float)
    n = y.shape[0]
    k = 4

    # --- setup: genotype-only, done once ---
    genes = split_into_genes(Z, G)
    F_list, P = setup_pooled(genes, r=r)
    Kaapply = lambda B: compute_KU(Z, B)
    Kdapply = lambda B: compute_KU(Zd, B)
    Wapply = lambda B: compute_WU_pooled(genes, F_list, P, B)

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

    KaU = Kaapply(U)  # K_a U for the fixed probes: genotype-only, formed once
    KdU = Kdapply(U)  # K_d U   "        "
    WU = Wapply(U)    # W U     "        "

    for it in range(iters):
        s2a, s2d, s2gxg, s2e = s
        matvec = lambda B: _v_matvec(Z, Zd, genes, F_list, P,
                                     s2a, s2d, s2gxg, s2e, B)

        # --- x = V^{-1} y ---
        xbuf = _cg_batched(matvec, yc, x0=xbuf, tol=cg_tol, maxiter=cg_maxiter)
        x = xbuf[:, 0]

        # data quadratics x'K_i x
        Kax = Kaapply(x)
        Kdx = Kdapply(x)
        Wx = Wapply(x)
        xKax = x @ Kax
        xKdx = x @ Kdx
        xWx = x @ Wx
        xIx = x @ x

        # --- score traces tr(V^{-1} K_i): Hutchinson, exact probe solves ---
        Pbuf = _cg_batched(matvec, U, x0=Pbuf, tol=cg_tol, maxiter=cg_maxiter)
        trV1Ka = np.mean(np.sum(Pbuf * KaU, axis=0))
        trV1Kd = np.mean(np.sum(Pbuf * KdU, axis=0))
        trV1W = np.mean(np.sum(Pbuf * WU, axis=0))
        trV1I = np.mean(np.sum(Pbuf * U, axis=0))

        score = np.array([0.5 * (xKax - trV1Ka),
                          0.5 * (xKdx - trV1Kd),
                          0.5 * (xWx - trV1W),
                          0.5 * (xIx - trV1I)])

        # --- average information: A_ij = 0.5 (K_i x)' V^{-1}(K_j x) ---
        KX = np.column_stack([Kax, Kdx, Wx, x])
        Gbuf = _cg_batched(matvec, KX, x0=Gbuf, tol=cg_tol, maxiter=cg_maxiter)
        AI = 0.5 * (KX.T @ Gbuf)
        AI = 0.5 * (AI + AI.T)

        # --- damped, bounded AI-Newton step ---------------------------------
        # W can still be nearly collinear with I, K_a with W, and K_d with I
        # (dominance deviations are close to independent noise at low MAF), so
        # the AI matrix can be near-singular and an undamped step explodes.
        # Levenberg-Marquardt ridge (scaled to AI) + trust region on the step +
        # a box clamp keep the path stable without perturbing well-identified
        # cases.  The constants are kept IDENTICAL to the two-component
        # pipeline's so a run is comparable to it line for line, even though
        # the extra components make the damping bind somewhat more often.
        dA = np.abs(np.diag(AI))
        ridge = lm * (dA.mean() + 1e-12)
        step = np.linalg.solve(AI + (ridge + jitter) * np.eye(k), score)

        mx = np.abs(step).max()
        max_step = step_frac * vary                  # trust region
        if mx > max_step:
            step *= max_step / mx
        s = np.clip(s + step, 1e-9, s_upper)         # box clamp

        if verbose:
            print(f"iter {it:2d}  s={s}  max|step|={np.abs(step).max():.3e}")
        if np.abs(step).max() < tol:
            break

    return s, AI


def MC_REML(Z, Zd, y, G, iters=30, Nmc=50, cg_tol=1e-6, cg_maxiter=1000,
            seed=None, r=R_DEFAULT, verbose=False):
    """Wrapper: returns (s2a_hat, s2d_hat, s2gxg_hat, s2e_hat, AI)."""
    s, AI = mc_reml(Z, Zd, y, G, iters=iters, Nmc=Nmc, cg_tol=cg_tol,
                    cg_maxiter=cg_maxiter, seed=seed, verbose=verbose, r=r)
    s2a_hat, s2d_hat, s2gxg_hat, s2e_hat = s
    return s2a_hat, s2d_hat, s2gxg_hat, s2e_hat, AI


# ----------------------------------------------------------- extra check
def variance_explained(Z, r):
    """Cumulative share of genotypic variance carried by the leading r PCs.

    rho(r) = sum_{s<=r} lam_s / sum_s lam_s,  lam_s the eigenvalues of Z Z'.

    The denominator is tr(Z'Z) = ||Z||_F^2 and needs no factorization, so only
    the leading r singular values are computed.  Returns a scalar.
    """
    full = min(Z.shape)
    r = min(r, full)
    if r >= full - 1:                       # ARPACK needs k < min(n, m)
        s = np.linalg.svd(Z, compute_uv=False)[:r]
    else:
        s = svds(Z, k=r, return_singular_vectors=False, v0=np.ones(full))
    return (s ** 2).sum() / np.einsum('ij,ij->', Z, Z)
