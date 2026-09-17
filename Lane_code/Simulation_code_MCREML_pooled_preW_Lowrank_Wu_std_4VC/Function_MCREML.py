# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.linalg import cholesky
from scipy.sparse.linalg import svds
import time

####################################################################
# FOUR variance components: additive + dominance + pooled within-gene
# epistasis + noise.  The epistasis kernel is C-NORMALIZED (this is the ONE
# thing that differs from the _unstd_4VC sibling; see NORMALIZATION below).
#
# Model (no fixed effects):
#     y = g_a + g_d + g_gxg + e ,
#     V = Var(y) = s2a K_a + s2d K_d + s2gxg W + s2e I ,
#     K_a = Z_a Z_a' / m                  (additive GRM, standardized dosages),
#     K_d = Z_d Z_d' / m                  (dominance GRM, standardized GCTA
#                                          dominance coding),
#     W   = W_raw / c ,
#     W_raw = (1/P) sum_g sum_{a<b in g} h_ab h_ab' ,  h_ab = Z_a .* Z_b ,
#           P = sum_g C(m_g, 2)           (pooled within-gene epistasis GRM,
#                                          UNSTANDARDIZED interactions).
#
# NORMALIZATION -- THE POINT OF THIS PIPELINE.
# The raw pooled kernel makes the REML component s2gxg and the REALIZED
# variance of the epistasis genetic value two different numbers:
#
#     E[ Var-hat(H gamma) ] = c * s2gxg ,
#     c = (1/P) sum_g sum_{a<b in g} Var-hat(Z_a .* Z_b) ,
#
# with c != 1 because the interaction columns h_ab are not standardized.  The
# _unstd_4VC sibling leaves that gap open and closes it AFTER the fit, by
# reporting V_l-hat = c-hat * s2gxg-hat.  Here the gap is closed BEFORE
# anything is drawn, by dividing the kernel itself once:
#
#     W = W_raw / c-hat   ==>   E[ Var-hat(H gamma) ] = (c / c-hat) s2gxg .
#
# So with s2a = s2d = s2gxg = 0.1 and s2e = 0.7 all four components have their
# realized variance centred on the value that was asked for, and no post-fit
# correction is applied to anything:  c_a = c_d = 1 already (both designs are
# column-standardized), c_e = 1 trivially, and c_gxg = c / c-hat, which is 1
# up to the accuracy of the estimator below.  s2gxg-hat is read directly as
# the realized epistasis variance.
#
# WHICH c-hat -- THE THREE ROUTES.  c is a GENOTYPE-ONLY constant and this
# module carries three ways to get it (realized_variance.pdf).  The one every
# kernel uses is fixed by the module constant C_METHOD and is the THIRD-MOMENT
# plug-in.  All three share the identity Var(Z_a Z_b) = 1 + r_ab s_a s_b for
# the per-pair variance; they differ in where the per-SNP skewness s comes
# from, and 'exact' skips the identity altogether:
#
#   'moment'  (C_METHOD -- THE ONE THE PIPELINE USES.)  s-hat_a = mean_t
#             (Z_ta^3), the EMPIRICAL third moment of the standardized column.
#             O(nm) time, O(n) storage, one mat-vec per gene -- no m-by-m R
#             and no n-by-P H -- and NO HWE assumption: the skewness is read
#             off the genotype rather than predicted from its allele
#             frequency, so departures from HWE move s-hat with the data
#             instead of biasing it.  Its only error is O(1/sqrt n) moment
#             sampling, the same order the note already accepts in c itself.
#
#   'hwe'     the note's closed form: the same formula with s_a = (1 - 2 p_a)
#             / sqrt(2 p_a (1 - p_a)) from the allele frequency, which is the
#             skewness of Z_a ONLY under Hardy-Weinberg.  Same O(nm) cost, and
#             it pays the HWE assumption on top of the sampling error.
#
#   'exact'   the literal mean per-pair sample variance, no identity and no
#             assumption, O(n sum_g m_g^2).  QUADRATIC in the gene size, i.e.
#             the cost the O(nm) routes exist to avoid, so it is not used as
#             the divisor.  It is the YARDSTICK the other two are scored
#             against and is reported next to every run.
#
# Keeping the divisor O(nm) is the point: c-hat then needs nothing a real
# analysis could not compute, and no number in this pipeline depends on a pass
# that is affordable only in simulation.  The residual scale gap c / c-hat is
# not hidden -- Simulate_Cholesky.py writes all three c's to
# result/c_<FILENAME>.txt, and the _matfree_ sibling's verify_lowrank.py measures the gap on the
# kernel itself.
#
# BOTH SIDES DIVIDE BY THE SAME c-hat.  The simulation Cholesky-factorizes
# s2gxg W_raw / c-hat and the estimator applies the rank-r truncation of the
# same W_raw / c-hat, both taking the divisor from pooled_c() -- ONE
# deterministic function of the same standardized genotype.  So the residual
# c / c-hat gap is COMMON to the two sides and cancels out of any comparison
# between them: simulation and estimation still differ by the TRUNCATION
# ALONE -- the invariant of this family -- and normalizing cannot mask
# truncation bias, a common positive scale factor changing neither the
# relative operator error nor which eigenvalues are discarded.
#
# WHAT NORMALIZING DOES NOT FIX.  tr(W) is still not n and W 1 is still not 0:
# dividing by a scalar cannot centre a kernel or restore a trace identity.  The
# realized-variance scale is the only thing c fixes, and it is the only thing
# claimed.
#
# THE EPISTASIS KERNEL IS BUILT FROM Z_a ONLY.  h_ab = Z_a .* Z_b pairs
# ADDITIVE columns; dominance enters through K_d and through nothing else.
# Dominance-by-dominance and additive-by-dominance interactions are NOT in this
# model -- they would be new kernels, not a new component on the existing W.
#
# COSTS.  K_a and K_d each enter as the thin mat-vec  K_i B = Z_i(Z_i'B)/m,
# O(n m c) -- strictly cheaper than the O(n^2 c) dense epistasis apply they sit
# next to, so the extra components add no new order of cost to the estimator.
# Only the SIMULATION pays for them in n-by-n terms: two more dense Cholesky
# factors.  The estimator holds TWO n-by-m designs and ONE n-by-n W.
#
# IDENTIFIABILITY.  K_a, K_d and W are far less collinear with each other than
# W is with I, but the components are not free: the AI matrix now has a
# 4-dimensional curvature to resolve, and small n / large m runs that were
# merely noisy in the two-component fit can become weakly identified here.
# K_d is the weakest of the three in practice -- dominance deviations are small
# and, at low MAF, their GRM is close to I -- so it is the component most
# likely to sit on the boundary.  The optimizer's Levenberg-Marquardt ridge +
# trust region + box clamp are unchanged and carry the same known boundary
# defect (see mc_reml).  NOTE the normalization RESCALES the epistasis
# direction: W = W_raw / c has a different norm from W_raw, so s2gxg lives on a
# different scale and the AI matrix is conditioned differently.  A run here is
# therefore NOT numerically step-for-step identical to the _unstd_4VC sibling's
# -- the two fits are reparametrizations of the same likelihood (s2gxg here
# equals c * s2gxg there at the optimum), not the same iteration.
#
# PRECOMPUTED W -- WHERE THIS DIRECTORY DIFFERS FROM ITS _matfree_ SIBLING.
# Everything above, and every line of the optimizer below, is the sibling's.
# The ONE change is how REML applies the epistasis kernel: the sibling rebuilds
# W u on every CG iteration from the O(n m r) low-rank LINEAR OPERATOR, while
# here the Cholesky job forms the kernel ONCE as a dense n-by-n matrix, caches
# it, and mc_reml applies it as a plain gemm  W @ B.  By default (r <= 0) the
# cached matrix is the EXACT normalized W the phenotype was drawn from, so
# simulation and estimation share the kernel outright and there is no
# truncation bias; r > 0 caches the dense matrix of the sibling's rank-r
# operator instead.  See the PRECOMPUTED estimation kernel section below.
####################################################################


# --------------------------------------------- operator-apply counting
# HOW MUCH LINEAR ALGEBRA ONE REPLICATE ACTUALLY DOES.  Apart from loading the
# precomputed W, this estimator's entire cost is APPLIES of the four linear
# operators K_a, K_d, W and I to a block of vectors, almost all of them
# inside the CG solves.  The wall-clock time already recorded per replicate
# mixes that count with BLAS threading, cache behaviour and machine load; the
# count below is the machine-independent half -- rerun the same replicate on
# any hardware and it comes out identical, because CG's iteration count is a
# deterministic function of (V, rhs, tol).
#
# WHAT IS COUNTED.  Two numbers per operator:
#
#   <op>_applies   how many times the operator was CALLED, each call on an
#                  (n, k) block;
#   <op>_columns   the total number of n-vectors pushed through it, i.e.
#                  sum of k over those calls.
#
# The columns are the cost-bearing number -- every one of the four operators
# is linear in k (a gemm pair against an n-by-m design for K_a and K_d, one
# n-by-n gemm for the precomputed W) -- while the applies say how well that work
# was batched.  A solve with Nmc = 50 right-hand sides is ONE apply and fifty
# columns, and costs fifty columns' worth of flops in one gemm rather than
# fifty.
#
# THE OPERATORS.
#   V   the composite s2a K_a + s2d K_d + s2gxg W-hat + s2e I, one _v_matvec
#       call.  This is "the" operator apply of the algorithm: CG sees only V.
#   K   compute_KU, i.e. K_a AND K_d together -- one routine serves both
#       designs and is not told which it holds.  Every V apply makes exactly
#       two K applies, so K_applies = 2 * V_applies + (the direct ones
#       mc_reml makes outside CG).
#   W   compute_WU_dense, the precomputed epistasis matrix.  Every V apply
#       makes exactly one, and it is the expensive one: O(n^2 k) against the
#       O(n m k) of a K apply.  Same key as the _matfree_ sibling's operator,
#       so the two pipelines' op-count files compare line for line.
# The identity operator costs a scalar multiply and is not counted.
#
# Also recorded: cg_solves (calls to _cg_batched), cg_iters (its inner
# iterations summed over those calls) and reml_iters (AI-REML iterations
# actually taken, which is <= iters because of the convergence break).
#
# The counters are GLOBAL and CUMULATIVE; mc_reml zeroes them on entry, so
# after one MC_REML call get_op_counts() describes exactly one replicate.
# Cost is two dict updates per apply, against a gemm -- unmeasurable.
_OP_COUNTS = {}


def reset_op_counts():
    """Zero every operator-apply counter.  mc_reml calls this on entry."""
    _OP_COUNTS.clear()


def get_op_counts():
    """Snapshot of the counters as a plain dict, safe to keep or write out.

    Read it AFTER a mc_reml / MC_REML call and it describes that call alone.
    Absent keys mean zero: nothing pre-registers a counter it never bumps.
    """
    return dict(_OP_COUNTS)


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
    _count_apply('K', U.shape[1])            # K_a and K_d share this routine
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
    C(m_g, 2)-term pair sum to machine precision (see the _matfree_ sibling's verify_lowrank.py).
    """
    n, mg = Zg.shape
    pg = mg * (mg - 1) // 2
    K = Zg @ Zg.T
    D = Zg * Zg
    S = 0.5 * (K * K - D @ D.T)
    return S, pg


def pooled_c_exact(Z_list):
    """The EXACT realized-variance factor c of the raw pooled kernel -- the
    YARDSTICK the O(nm) plug-ins are scored against.

        c = (1/P) sum_g sum_{a<b in g} Var-hat(Z_a .* Z_b) ,
        P = sum_g C(m_g, 2) ,

    Var-hat being the sample variance with ddof=0.  It is the constant that
    relates the REML variance COMPONENT to the REALIZED variance of the
    epistasis genetic value: for gamma ~ N(0, (s2gxg/P) I_P),

        E[ Var-hat(H gamma) ] = c * s2gxg ,

    so drawing from W_raw / c instead of W_raw makes that expectation exactly
    s2gxg.  Equivalently  c = tr(P_c W_raw) / n  with P_c = I - 11'/n, which is
    how the _matfree_ sibling's verify_lowrank.py re-derives it from the dense kernel.

    THIS IS NOT THE DIVISOR.  The kernels take theirs from pooled_c(), which
    defaults to the O(nm) third-moment plug-in (C_METHOD); see the module
    header.  What this function gives is the number that plug-in ESTIMATES, so
    c_exact / c-hat is the residual scale error of a run -- reported next to
    every one by Simulate_Cholesky.py and measured on the kernel itself by
    the _matfree_ sibling's verify_lowrank.py.

    Takes the list of COLUMN-STANDARDIZED gene blocks -- the same Z_g the
    kernel is built from, not the raw dosages.
    compute_c_pooled(method='exact') is this function with the
    standardize-and-split done for you; the two agree to the last bit because
    it delegates here.

    Formed WITHOUT the n-by-P interaction design H.  With D_g = Z_g .* Z_g and
    the per-gene Gram matrix Z_g' Z_g,

        sum_{a<b} mean(Z_a^2 Z_b^2) = (1 / 2n)  ( ||D_g 1||^2 - sum_{t,a} Z_ta^4 ) ,
        sum_{a<b} mean(Z_a Z_b)^2   = (1 / 2n^2)( ||Z_g' Z_g||_F^2 - sum_a (Z_a' Z_a)^2 ) ,

    i.e. O(n sum_g m_g^2) time and O(max_g m_g^2) storage: one gemm per gene.
    Cheap enough to run once per job here, but QUADRATIC in the gene size --
    the cost the plug-in routes exist to avoid, and the reason it is a
    diagnostic rather than the divisor.

    A 1-SNP gene contributes no pair and is skipped, exactly as in
    build_W_pooled and setup_pooled, so c averages over the same P pairs the
    kernel sums over.
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
        # c is an average of sample variances, so this can only happen if every
        # interaction column is constant (a degenerate / monomorphic genotype).
        # Dividing by it would silently produce inf or a sign flip.
        raise ValueError(f"Pooled realized-variance factor c = {c!r} is not "
                         f"positive; the genotype has no varying interaction "
                         f"column. Check the MAF filter.")
    return c


# ---------------------------------------- c: the O(nm) plug-in estimators
# The two plug-ins are the SAME closed form (realized_variance.pdf) fed a
# different per-SNP skewness s.  Under the identity
#
#     Var(Z_a Z_b) = 1 + r_ab s_a s_b ,
#
# summing over the within-gene pairs and using that R_g = Z_g' Z_g / n is
# symmetric with UNIT DIAGONAL (the columns are standardized, ddof=0),
#
#     c-hat = (1/P) sum_g sum_{a<b} (1 + r_ab s_a s_b)
#           = 1 + (1 / 2P) sum_g ( s_g' R_g s_g - ||s_g||^2 ) ,
#     s_g' R_g s_g = ||Z_g s_g||^2 / n ,
#
# i.e. ONE mat-vec per gene: O(n m_g) time and O(n) storage, with neither the
# m-by-m correlation matrix R nor the n-by-P interaction design H ever formed.
# Both are strictly cheaper than the O(n sum_g m_g^2) exact route and than
# anything else either job does.
#
# The routes differ ONLY in where s comes from: third_moment_skewness reads it
# off the standardized genotype (no assumption), hwe_skewness predicts it from
# the allele frequency (Hardy-Weinberg).  See the module header for why the
# former is the default.

C_METHOD = 'moment'     # the route pooled_c() -- and hence every kernel in
                        # this pipeline -- uses.  'moment' | 'hwe' | 'exact'.
                        # All three stay callable; only this one is ever on a
                        # code path that touches a kernel.


def third_moment_skewness(Z):
    """Per-SNP skewness read straight off the data, with NO HWE assumption:

        s-hat_i = (1/n) sum_t Z_ti^3 ,

    the third moment of the COLUMN-STANDARDIZED design.  Because those columns
    have mean 0 and sample variance 1 (ddof=0), the raw third moment IS the
    standardized skewness -- no centering or scaling is left to do.  O(nm)
    time, O(m) storage, and formed by an einsum so no n-by-m cube is
    materialized.

    THE PIPELINE'S ROUTE (C_METHOD = 'moment').  Against hwe_skewness it drops
    Hardy-Weinberg entirely: s is ESTIMATED from the genotype instead of
    PREDICTED from its allele frequency, so inbreeding, population structure
    and genotyping artefacts move s-hat with the data rather than biasing it.
    The price is O(1/sqrt n) sampling error in each s-hat_i -- the same order
    as the sample-variance error the note already accepts inside c itself, and
    it does not grow with m.

    Takes the standardized design (additive_design's output, or one of its
    gene blocks -- slicing columns changes neither their mean nor their
    variance), the SAME matrix the kernel is built from.  A monomorphic column
    is identically 0 after standardization and gets s-hat_i = 0; every pair it
    enters then contributes 1 + 0 = 1 to c-hat while its true contribution is
    Var-hat = 0, the one failure mode this route shares with the HWE one.  Use
    MAF-filtered genotypes, as the rest of the pipeline assumes.
    """
    Z = np.asarray(Z, dtype=float)
    return np.einsum('ij,ij,ij->j', Z, Z, Z) / Z.shape[0]


def hwe_skewness(real_data):
    """Skewness of every column-standardized SNP under HWE, from allele
    frequency alone:

        s_i = (1 - 2 p_i) / sqrt(2 p_i (1 - p_i)) ,   p_i = mean(X_i) / 2 ,

    X being the 0/1/2 dosage matrix.  O(nm).

    This is the third central moment of a Binomial(2, p_i) divided by its
    variance to the 3/2 -- exactly what third_moment_skewness measures
    empirically, but only WHEN HWE HOLDS.  The gap between the two is the HWE
    departure, and it is why 'moment' rather than 'hwe' is C_METHOD.

    Takes the RAW dosages: the allele frequency is the one thing the
    standardized design no longer carries.

    A monomorphic SNP (p_i in {0, 1}) has no defined s_i and gets 0.  Such a
    column is identically 0 after standardization, so every pair it enters has
    Var-hat = 0 exactly, whereas the closed form then counts that pair as
    1 + 0 = 1: an O(#monomorphic / m) relative error in c.  Use MAF-filtered
    genotypes, as the rest of the pipeline assumes.
    """
    X = np.asarray(real_data, dtype=float)
    p = X.mean(axis=0) / 2.0
    denom = 2.0 * p * (1.0 - p)
    ok = denom > 0.0
    s = np.zeros(p.shape[0])
    s[ok] = (1.0 - 2.0 * p[ok]) / np.sqrt(denom[ok])
    return s


def _split_like_genes(v, Z_list):
    """Cut a length-m per-SNP vector into the gene blocks of Z_list.

    Z_list comes from split_into_genes, whose blocks are CONTIGUOUS and in
    column order, so slicing v by the same block widths pairs s_i with the
    gene its SNP sits in.  Guards the total so a mismatched G can never
    silently misalign the skewness against the genotype.
    """
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

    One mat-vec per gene and nothing m-by-m or n-by-P.  A 1-SNP gene
    contributes no pair and is skipped, exactly as in build_W_pooled,
    setup_pooled and pooled_c_exact, so the average runs over the same P pairs
    the kernel sums over.
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
        # c-hat is an ESTIMATE, not an average of sample variances, so unlike
        # pooled_c_exact it can go non-positive: a strongly negative skewness
        # cross term drives the sum below -2P.  Dividing the kernel by it would
        # flip the sign of the epistasis component, so refuse instead.
        raise ValueError(f"Plug-in realized-variance factor c-hat = {c!r} "
                         f"(method={method!r}) is not positive; the skewness "
                         f"term overwhelms the leading 1.  Check the MAF "
                         f"filter, or fall back to method='exact'.")
    return c


def pooled_c_moment(Z_list):
    """c-hat by the THIRD-MOMENT plug-in -- the divisor this pipeline uses.

    s-hat from third_moment_skewness (empirical, no HWE), fed to the O(nm)
    closed form.  Takes the COLUMN-STANDARDIZED gene blocks, the same argument
    pooled_c_exact takes and the same matrices the kernel is built from, so
    the divisor is a deterministic function of exactly what is being
    normalized -- which is what lets simulation and estimation arrive at the
    same number without passing it between them.
    """
    return _pooled_c_plugin(Z_list,
                            [third_moment_skewness(Zg) for Zg in Z_list],
                            'moment')


def pooled_c_hwe(Z_list, real_data):
    """c-hat by the HWE closed form -- the note's estimator, DIAGNOSTIC here.

    s from hwe_skewness (allele frequency, Hardy-Weinberg), fed to the same
    O(nm) closed form as pooled_c_moment; the two differ in nothing else, so
    their gap is the HWE departure alone.  Needs the RAW dosages alongside the
    standardized blocks because the allele frequencies are not recoverable
    from the latter.
    """
    s = hwe_skewness(real_data)
    return _pooled_c_plugin(Z_list, _split_like_genes(s, Z_list), 'hwe')


def pooled_c(Z_list, method=C_METHOD, real_data=None):
    """THE DIVISOR.  Every kernel in this pipeline gets its c from here.

    Z_list is the list of COLUMN-STANDARDIZED gene blocks; method is one of

        'moment'  third-moment plug-in, O(nm)   -- the default (C_METHOD)
        'hwe'     HWE closed form,       O(nm)  -- needs real_data
        'exact'   mean per-pair sample variance, O(n sum_g m_g^2)

    all documented in the module header.  build_W_pooled (simulation) and
    setup_pooled (estimation) both call this with the default and with gene
    blocks derived from the same genotype by the same route, so they cannot
    disagree about the divisor -- there is nothing to pass between the two
    jobs and nothing to cache.

    method is a keyword with a module-level default ON PURPOSE: the pipeline
    exposes no switch for it (no CLI flag, no environment variable), so a run
    cannot half-change routes.  Editing C_METHOD changes both sides at once,
    or neither.
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
    """Pooled WITHIN-gene pairwise-epistasis GRM, UNSTANDARDIZED interactions,
    C-NORMALIZED.

        W_raw = (1/P) sum_{g=1}^G sum_{a<b in g} h_ab h_ab' ,  h_ab = Z_a .* Z_b ,
        P     = sum_{g=1}^G C(m_g, 2) ,
        W     = W_raw / c-hat ,   c-hat = pooled_c(Z_list) .

    The division by c-hat is the ONE change from the _unstd_4VC sibling.  It
    makes E[Var-hat(H gamma)] = (c / c-hat) s2gxg for gamma ~ N(0, (s2gxg/P)
    I_P), so the phenotype this kernel generates realizes the variance that was
    asked for up to the plug-in's own error (see the module header).  c-hat is
    a positive genotype-only scalar, so W stays symmetric PSD and its
    eigenvectors are untouched -- only the scale moves.

    A 1-SNP gene contributes no pair and is skipped.  n-by-n storage; O(n^2 m)
    time, plus the O(nm) of c-hat, which is negligible beside it.  Built by the
    SIMULATION (the Cholesky factor needs an explicit matrix), which also
    caches it as the PRECOMPUTED estimation kernel when r <= 0; with r > 0 the
    cached kernel is instead build_W_lowrank_dense's rank-r truncation of this
    same normalized W, dividing by the SAME c-hat.

    With return_c=True returns (W, c) -- the divisor is worth recording next to
    the run (Simulate_Cholesky.py writes it out), since it is the number that
    puts s2gxg and the realized variance on one scale.

    tr(W) is NOT n here, and normalizing does not make it so: c-hat estimates
    the mean CENTERED per-pair variance, which is not the trace identity.  Nor is
    W 1 = 0 -- there is still no centering.  K_a and K_d, which ARE built from
    standardized designs, keep tr(K) = n exactly, so the genetic components
    remain on different TRACE scales even though they are now on a common
    REALIZED-VARIANCE scale.  The latter is what the estimates are read on.
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


# ----------------------------------- PRECOMPUTED estimation kernel (dense)
# THE ONE THING THIS PIPELINE CHANGES against the _matfree_ sibling.  There,
# REML applies the epistasis kernel through the O(n m r) low-rank LINEAR
# OPERATOR compute_WU_pooled, rebuilt from per-gene SVD factors inside every CG
# iteration.  Here the kernel is formed ONCE as a dense n-by-n matrix by the
# Cholesky job, cached to disk, loaded by every replicate, and applied as a
# plain  W @ B  (O(n^2 k) per apply).  Everything else -- the c-normalization,
# the K_a / K_d applies, Hutchinson traces, the two-phase schedule, the trust
# region, the feasibility guard, the operator counters -- is the sibling's code.
#
# WHICH MATRIX IS PRECOMPUTED is set by r:
#
#   r <= 0  (R_DEFAULT)  the EXACT c-normalized kernel W = W_raw / c-hat -- the
#                        very matrix the Cholesky job factorizes to draw g_gxg,
#                        saved before it is freed.  Simulation and estimation
#                        then use literally the same W and there is NO
#                        truncation bias at all; W is PSD by construction.
#   r >  0               the dense matrix of the sibling's rank-r operator,
#                        W-hat = (1/(2 P c)) sum_g [ (Q_g Lam_g Q_g') .* K_g
#                        - D_g D_g' ].  Same matrix the sibling's operator
#                        applies (to machine precision, see verify_preW.py), so
#                        a run here reproduces the sibling's truncation bias
#                        with a dense apply in place of the operator.  NOT PSD
#                        in general, exactly as in the sibling.
#
# THE CACHE IS THIS PIPELINE'S OWN.  It lives under <DIR>/W/ and its name
# carries (mode, n, m, G, r), never the shared stored_genotype/W_*.npy that the
# dense gxg family writes if-not-exists and keys by genotype only -- a kernel
# from another pipeline must never be picked up silently.  It is OVERWRITTEN
# (atomically) by every Cholesky job rather than written if-not-exists, so the
# W a run fits is always the W its own Cholesky job built.  A JSON sidecar
# records (mode, n, m, G, r, c-hat, P, psd) and the estimator refuses a kernel
# whose sidecar disagrees with its arguments, or whose c-hat differs from the
# one it recomputes from the genotype.
R_DEFAULT = 0           # 0 = exact precomputed W; r > 0 = dense rank-r W-hat
W_CACHE_KERNEL = "pooled_unstd_cnorm_4VC"   # written into the sidecar


def w_cache_tag(r):
    """'rexact' for the exact kernel (r <= 0), else 'r<r>'."""
    return "rexact" if r <= 0 else f"r{int(r)}"


def w_cache_paths(root, mode, n, m, G, r):
    """(npy, json) paths of the cached estimation kernel for one genotype,
    gene split and truncation level."""
    base = f"{root}/W_{mode}_n{n}_m{m}_G{G}_{w_cache_tag(r)}"
    return base + ".npy", base + ".json"


def setup_pooled(Z_list, r):
    """Per-gene rank-r factors of K_g = Z_g Z_g' -- the sibling's SETUP phase.

    Returns (F_list, P, c):

        F_list : per-gene state, one dict per gene (None for a 1-SNP gene)
        P      : sum_g C(m_g, 2)  -- the global pair total
        c      : the realized-variance factor c-hat (pooled_c, i.e. the
                 O(nm) third-moment plug-in), the SAME divisor build_W_pooled
                 uses, so the truncated kernel is a truncation of the SAME
                 normalized W the simulation drew from.

    The leading eigenpairs of K_g come from a truncated SVD of Z_g; K_g is
    never formed here.  ARPACK needs k < min(n, m_g) and returns singular
    values in ASCENDING order, so a gene small enough that r reaches full rank
    is factorized exactly instead -- which is also where W-hat becomes EXACT.
    Used here only to BUILD the dense W-hat (r > 0); the REML loop never sees
    these factors.
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
    return F_list, P, pooled_c(Z_list)


def build_W_lowrank_dense(Z_list, r, return_c=False):
    """DENSE rank-r truncated kernel  W-hat,  the matrix of the _matfree_
    sibling's compute_WU_pooled:

        (K .* K) u ~ sum_{s<=r} lam_s q_s .* (K (q_s .* u))
                   = ((Q Lam Q') .* K) u ,

    so, per gene,  2 S-hat_g = (Q_g Lam_g Q_g') .* K_g - D_g D_g'  and

        W-hat = (1 / (2 P c)) sum_g 2 S-hat_g ,   c = pooled_c(Z_list) .

    O(n^2 (m_g + r)) per gene and O(n^2) storage -- the same order as the exact
    build, and paid ONCE by the Cholesky job.  Symmetrized at the end: the two
    Hadamard factors are symmetric, but (Q Lam) Q' is a general gemm and is not
    bitwise symmetric, and CG and the spectral-range routine both assume a
    symmetric matrix.  The change is at round-off level.

    At r >= min(n, m_g) for every gene W-hat equals build_W_pooled's W to
    machine precision.  Below that it is symmetric but NOT PSD in general.
    """
    F_list, P, c = setup_pooled(Z_list, r=r)
    n = Z_list[0].shape[0]
    W = np.zeros((n, n))
    for gene in F_list:
        if gene is None:                         # 1-SNP gene: no pair
            continue
        Zg, Dg, Q, lam = gene['Z'], gene['D'], gene['Q'], gene['lam']
        S = Zg @ Zg.T                            # K_g
        S *= (Q * lam) @ Q.T                     # (Q Lam Q') .* K_g
        S -= Dg @ Dg.T                           # the a = b terms
        W += S
        del S
    W /= (2.0 * P * c)
    W += W.T
    W *= 0.5
    return (W, c) if return_c else W


def save_W_cache(W, meta, npy_path, json_path):
    """Write the estimation kernel and its sidecar ATOMICALLY.

    Each file goes to a unique temporary name in the same directory and is
    os.replace'd into place, so a reader never sees a half-written array and
    two Cholesky jobs on the same genotype (different variance targets, same
    kernel) cannot interleave their bytes.  Overwrites any previous kernel.
    """
    import json
    import os
    import tempfile
    d = os.path.dirname(npy_path)
    os.makedirs(d, exist_ok=True)
    for path, write in ((npy_path, lambda f: np.save(f, W)),
                        (json_path, lambda f: f.write(
                            json.dumps(meta, indent=2).encode()))):
        fd, tmp = tempfile.mkstemp(dir=d, prefix=".tmp_W_")
        try:
            with os.fdopen(fd, 'wb') as f:
                write(f)
            os.replace(tmp, path)
        except BaseException:
            if os.path.exists(tmp):
                os.remove(tmp)
            raise


def load_W_cache(npy_path, json_path, mode, n, m, G, r, c_expected):
    """Load a cached estimation kernel, refusing one that is not THIS run's.

    Checks the sidecar against (kernel, mode, n, m, G, r), the array shape
    against (n, n), and the recorded c-hat against c_expected -- the divisor
    the caller recomputed from the genotype it is about to fit with (O(nm), one
    mat-vec per gene).  The last check is what catches a kernel left over from
    a genotype file that has since been regenerated under the same name.

    Returns (W, meta).
    """
    import json
    import os
    for p in (npy_path, json_path):
        if not os.path.exists(p):
            raise FileNotFoundError(
                f"precomputed kernel file not found: {p}\n"
                f"Run the Cholesky step with the same (mode, n, m, G, r) first.")
    with open(json_path) as f:
        meta = json.load(f)
    want = {'kernel': W_CACHE_KERNEL, 'mode': mode, 'n': int(n), 'm': int(m),
            'G': int(G), 'r': int(max(r, 0))}
    bad = {k: (meta.get(k), v) for k, v in want.items() if meta.get(k) != v}
    if bad:
        raise ValueError(f"{json_path} does not describe this run "
                         f"(field: (cached, wanted)) {bad}")
    c_cached = float(meta['c_norm'])
    if abs(c_cached - c_expected) > 1e-10 * abs(c_expected):
        raise ValueError(f"cached kernel was divided by c-hat={c_cached!r} but "
                         f"this genotype gives {c_expected!r}: the kernel was "
                         f"built from a different genotype.  Rerun the "
                         f"Cholesky step.")
    W = np.load(npy_path)
    if W.shape != (n, n):
        raise ValueError(f"{npy_path} has shape {W.shape}, expected {(n, n)}.")
    return W, meta


def compute_WU_dense(W, U):
    """Dense product  W @ U  with the PRECOMPUTED estimation kernel.

    U may be (n,) or (n, k); the result matches its shape.  O(n^2 k) -- one
    gemm -- against the sibling operator's O(n m r k).  Counted as one 'W'
    apply on k columns, exactly as the operator was, so the op-count files of
    the two pipelines line up key for key.
    """
    U = np.asarray(U, dtype=float)
    _count_apply('W', 1 if U.ndim == 1 else U.shape[1])
    return W @ U

# ------------------------------------------------------------- simulation
def simulate_Cholesky_4vc(real_data, G, s2a=0.1, s2d=0.1, s2gxg=0.1, s2e=0.7,
                          stability=1e-10, r=R_DEFAULT, save_W_est=None):
    """Cholesky factors of the additive, dominance and pooled-epistasis
    covariances.

        La La'     = s2a   K_a ,   K_a = Z_a Z_a' / m ,
        Ld Ld'     = s2d   K_d ,   K_d = Z_d Z_d' / m ,
        Lgxg Lgxg' = s2gxg W ,     W   = W_raw / c .

    The epistasis factor is built from the C-NORMALIZED W, so the g_gxg it
    draws has E[Var-hat(g_gxg)] = s2gxg exactly -- the same target the additive
    and dominance factors already hit, their designs being standardized.  That
    is the whole difference from the _unstd_4VC sibling; La and Ld are
    untouched.

    PRECOMPUTING THE ESTIMATION KERNEL.  save_W_est, if given, is called ONCE
    as  save_W_est(W_est, c, psd, build_time)  with the matrix REML will fit:

        r <= 0  W itself -- the very array Lgxg was factorized from, handed
                over before it is freed, so no second build and no chance of
                the two sides drifting apart.  psd = True.
        r >  0  build_W_lowrank_dense(genes, r), built after W is freed.
                psd = False (the truncation is not PSD in general).

    The callback is expected to write the array to disk; this function keeps
    no reference to it afterwards.

    Returns (La, Ld, Lgxg, w_build_time, c, w_est_build_time), w_build_time
    being the wall-clock cost of the EPISTASIS kernel alone -- the O(n^2 m)
    term that dominates this job and the one the sibling pipelines report --
    and w_est_build_time that of the cached estimation kernel (equal to
    w_build_time when r <= 0, the matrix being the same one).  K_a and K_d cost
    O(n^2 m) too but with a single gemm each and no Hadamard square, so they
    are folded into the untimed remainder rather than given their own records.

    The three factors are built SEQUENTIALLY and each dense GRM is released as
    soon as its factor exists, so the peak is ~4 n-by-n arrays (two finished
    factors + one GRM + the factorization workspace) rather than 7.  At n =
    16000 that is ~8 GB instead of ~14 GB -- still the memory-critical job of
    the pipeline, and the reason it gets its own larger SLURM allocation.
    Building W-hat for r > 0 happens while only Lgxg is held, so it stays
    inside the same peak.
    """
    Za = additive_design(real_data)
    Zd = dominance_design(real_data)
    n, m = Za.shape

    # --- epistasis first (the expensive one), then free W -------------------
    genes = split_into_genes(Za, G)          # several Z, one per gene

    # The c-normalization happens INSIDE build_W_pooled and is inside the
    # timed block: it is part of building this kernel, and pretending it is
    # free would misreport the simulation's cost.  It is one gemm per gene
    # against the O(n^2 m) Hadamard square, so the timing barely moves.
    t_start = time.perf_counter()
    W, c = build_W_pooled(genes, return_c=True)
    w_build_time = time.perf_counter() - t_start

    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)

    # --- the PRECOMPUTED estimation kernel, then free W ---------------------
    if r <= 0:
        w_est_build_time = w_build_time      # the SAME matrix, built once
        if save_W_est is not None:
            save_W_est(W, c, True, w_est_build_time)
        del W
    else:
        del W
        t_start = time.perf_counter()
        W_hat, c_hat = build_W_lowrank_dense(genes, r, return_c=True)
        w_est_build_time = time.perf_counter() - t_start
        if c_hat != c:
            raise RuntimeError(f"build_W_lowrank_dense divided by {c_hat!r} but "
                               f"build_W_pooled by {c!r}: the two kernels are "
                               f"not on the same scale.")
        if save_W_est is not None:
            save_W_est(W_hat, c, False, w_est_build_time)
        del W_hat

    # --- additive --------------------------------------------------------
    Ka = build_K(Za)
    La = cholesky(s2a * Ka + stability * np.eye(n), lower=True)
    del Ka

    # --- dominance -------------------------------------------------------
    Kd = build_K(Zd)
    Ld = cholesky(s2d * Kd + stability * np.eye(n), lower=True)
    del Kd

    return La, Ld, Lgxg, w_build_time, c, w_est_build_time


def _force_var(v, target, name):
    """Rescale v so its ddof=0 sample variance is exactly target.

    A target of 0 means the component is absent: its Cholesky factor is the
    zero matrix, the draw is identically zero, and sqrt(0/0) is not a scale
    factor -- return the zero vector, which already has the required variance.
    A nonzero target with a degenerate draw is a real fault (a rank-0 factor
    where one was expected) and raises rather than dividing by zero.
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
    """Draw one phenotype  y = g_a + g_d + g_gxg + e  from the four-component
    model,

        g_a   = La u1   ~ N(0, s2a K_a) ,   K_a = Z_a Z_a' / m ,
        g_d   = Ld u2   ~ N(0, s2d K_d) ,   K_d = Z_d Z_d' / m ,
        g_gxg = Lgxg u3 ~ N(0, s2gxg W) ,   W   = (1/(P c)) H H' ,
        e     = sqrt(s2e) u4 ~ N(0, s2e I) .

    The four draws are INDEPENDENT.  With force_realized=False that is the
    whole story and the phenotype has exactly the law the estimator fits.

    force_realized=True (THE DEFAULT) RESCALES EACH COMPONENT so its realized
    variance equals its target EXACTLY:

        a *= sqrt(s2a / Var-hat(a)) ,   ... and likewise d, gxg, e,

    Var-hat at ddof=0, the same normalisation as the realized variances
    described below.  Read what this does and does not do:

      - It removes the O(1/sqrt n) scatter of each realized variance around
        its target.  Every replicate now has Var-hat(a) = s2a, Var-hat(d) =
        s2d, Var-hat(gxg) = s2gxg and Var-hat(e) = s2e to machine precision,
        so an estimate's deviation from the nominal target is the ESTIMATOR's
        error alone and no longer the estimator's error plus the draw's.  That
        is the point of the switch.
      - It also absorbs the residual c / c-hat factor on the epistasis
        component.  E[Var-hat(g_gxg)] was s2gxg only up to the O(nm) plug-in's
        accuracy; forcing the realized variance makes it s2gxg exactly, so the
        plug-in error no longer reaches the phenotype at all.
      - IT CHANGES THE GENERATING LAW.  Each component is divided by the
        square root of its own sample variance, which is a random quantity
        correlated with the draw, so a is no longer N(0, s2a K_a) and y is no
        longer N(0, V).  The scaled vector lies on a sphere in the metric
        Var-hat(.) = s2a rather than being Gaussian around it.  REML is
        therefore fitting a model the simulation does not obey.  At these n
        the effect is small -- the scale factor is 1 + O(1/sqrt n) -- but it is
        a MISSPECIFICATION, not a variance reduction, and the estimator's
        sampling distribution under it is not the one REML theory describes.
      - IT MAKES THE REALIZED-VARIANCE COLUMNS DEGENERATE.  Columns 6-9 of the
        per-replicate result row are now constants, so the paired
        (estimate - realized) diagnostic in calc_stats.py collapses onto the
        unpaired (estimate - target) one and its standard deviation stops
        being a separate number.  Nothing breaks; the two summaries simply
        become the same summary.
      - A component whose target is ZERO cannot be rescaled (its draw is
        identically zero) and is left at zero, which is the right answer.

    force_realized=False restores the plain draw, and is what the sibling
    pipelines still do -- a run here is not step-for-step comparable to theirs
    with the default on.

    Since Lgxg Lgxg' = s2gxg W = (s2gxg / (P c)) H H', the draw g_gxg has
    EXACTLY the law of  H gamma  with  gamma ~ N(0, (s2gxg / (P c)) I_P) --
    equivalently  (H / sqrt(c)) gamma_0  with gamma_0 ~ N(0, (s2gxg/P) I_P),
    the effect-size model of realized_variance.pdf applied to the NORMALIZED
    interaction design -- without ever forming the n-by-P design H.  The same
    argument makes g_a exactly Z_a beta with beta ~ N(0, (s2a/m) I_m), and g_d
    exactly Z_d delta with delta ~ N(0, (s2d/m) I_m).

    REALIZED VARIANCE -- ALL FOUR COMPONENTS NOW AGREE WITH THEIR TARGET.
    The variance of a genetic value across individuals,

        V_ell := Var-hat(g_gxg) = (1/n) g_gxg' P_c g_gxg ,   P_c = I - 11'/n ,

    is random (it changes with every effect draw) and, for the RAW kernel, has
    expectation c * s2gxg with c != 1 (see pooled_c_exact).  Because W is
    divided by the estimate c-hat of that same c here,

        E[V_ell] = (c / c-hat) s2gxg ,

    which is s2gxg to the accuracy of the O(nm) plug-in (C_METHOD; the ratio
    is written out per run, see Simulate_Cholesky.py).  So the epistasis
    component joins the other three, which needed no correction to begin with:

        V_a := Var-hat(Z_a beta)  = var(a) ,   E[V_a] = s2a   (c_a = 1) ,
        V_d := Var-hat(Z_d delta) = var(d) ,   E[V_d] = s2d   (c_d = 1) ,
        V_e := Var-hat(e)         = var(e) ,   E[V_e] = s2e   (trivially) ,

    all at ddof=0, matching the note's 1/n normalisation.  With s2a = s2d =
    s2gxg = 0.1 and s2e = 0.7 the four realized variances are centred on
    0.1, 0.1, 0.1, 0.7 -- the first three exactly, the epistasis one up to
    c / c-hat.  That is the point of this pipeline, and the reason no post-fit
    c correction is applied anywhere downstream.

    Under force_realized=False they are still RANDOM: one draw is not its
    expectation, and each scatters around its target by O(1/sqrt n).  With
    return_realized=True all four are returned for THIS draw, so every estimate
    can be compared against what its own replicate actually realized rather
    than only against the nominal target.  Under force_realized=True the four
    returned values are the four targets, by construction.

    NOTE the four realized variances do NOT sum to var(y), and FORCING THEM
    DOES NOT FIX THAT.  a, d, gxg and e are independent in expectation but
    their SAMPLE cross-products are not zero, so

        var(y) = V_a + V_d + V_ell + V_e + 2 (cov-hat terms) ,

    and rescaling each component individually leaves those cross terms exactly
    where they were -- it pins the four diagonal terms and touches none of the
    six off-diagonal ones.  So var(y) is still not s2a + s2d + s2gxg + s2e, and
    its residual scatter is now entirely the cross terms' (O(1/sqrt n)).  Do
    not treat the four as a partition of the phenotypic variance under either
    setting.

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

    if force_realized:
        # Pin each component's realized variance to its target.  ddof=0, to
        # match Var-hat(.) = (1/n) g' P_c g as defined above and as the
        # realized columns report it.  See the docstring: this is a change of
        # generating law, not a variance reduction.
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
def compute_c_pooled(real_data, G, method=C_METHOD):
    """Realized-variance factor c of the RAW pooled kernel, from the raw 0/1/2
    dosage matrix -- the from-genotype entry point to the three routes.

        E[ Var-hat(H gamma) ] = c * s2gxg ,
        c = (1/P) sum_g sum_{a<b in g} Var-hat(Z_a .* Z_b) ,  P = sum_g C(m_g, 2) ,

    for the pooled WITHIN-gene unstandardized kernel: only within-gene pairs
    make up H, exactly the pairs build_W_pooled / setup_pooled count in P.

    This is pooled_c() with the standardize-and-split done for you -- it builds
    the additive design and cuts it into the same G contiguous genes the kernel
    uses, then dispatches on method:

        method='moment'  (DEFAULT, C_METHOD -- the number the kernels divide
                          by)  third-moment plug-in, O(nm).  s-hat_a =
                          mean_t(Z_ta^3) read off the genotype, no HWE.
        method='hwe'      the note's closed form, O(nm).  Same formula, s from
                          the allele frequency under Hardy-Weinberg.  Its gap
                          against 'moment' is the HWE departure.
        method='exact'    the literal mean per-pair sample variance,
                          O(n sum_g m_g^2).  The yardstick: its gap against
                          'moment' is the plug-in's total error, and therefore
                          the residual scale factor c / c-hat of a run.

    It takes the RAW dosages because the HWE route needs the allele
    frequencies, which the standardized Z no longer carries.  The Cholesky job
    calls all three and writes them to result/c_<FILENAME>.txt, which is the
    only place the comparison is made; the kernels never come through here.

    WHAT IT IS FOR HERE.  The _unstd_4VC sibling calls this AFTER the fit, to
    turn s2gxg-hat into a realized-variance estimate c-hat * s2gxg-hat.  This
    pipeline divides the KERNEL by c-hat instead, so the fit already lives on
    the realized-variance scale and NO post-fit correction is applied.

    THE DIVISOR IS NOT TAKEN FROM HERE.  The kernels call pooled_c on the
    standardized gene blocks directly, so simulation and estimation cannot
    disagree about which c they used and neither can drift from this function.
    method defaults to C_METHOD accordingly: an accidental call with the
    default returns the same number the kernels divide by.

    THE ADDITIVE AND DOMINANCE COMPONENTS NEED NO ANALOGUE.  Their designs are
    Z_a and Z_d, whose columns are standardized to sample variance 1 (ddof=0)
    by construction, so c_a = (1/m) sum_i Var-hat(Z_a,i) = 1 EXACTLY, and
    likewise c_d = 1 -- not approximately, and with no HWE assumption (the HWE
    assumption in the dominance CODING is a separate matter: it decides what
    the deviation means, not how its column is scaled).  After the epistasis
    normalization all three genetic components sit on one realized-variance
    scale, which is exactly what makes the four estimates directly comparable
    to the four targets.
    """
    Z = additive_design(real_data)
    genes = split_into_genes(Z, G)
    return pooled_c(genes, method=method, real_data=real_data)


# ------------------------------------------------------------ linear algebra
def _v_matvec(Z, Zd, W, s2a, s2d, s2gxg, s2e, B):
    """Apply V = s2a K_a + s2d K_d + s2gxg W + s2e I to B, K_a and K_d
    matrix-free, W the PRECOMPUTED dense C-NORMALIZED epistasis kernel.

        V B = (s2a/m) Z_a(Z_a'B) + (s2d/m) Z_d(Z_d'B)
              + s2gxg (W B) + s2e B ,     W = W_raw / c  (or its rank-r W-hat).

    B may be (n,) or (n, k); the result matches its shape.  The additive and
    dominance terms are a pair of gemms each against their design, exactly as
    in the _matfree_ sibling; the epistasis term is one dense gemm against the
    cached W, which already carries the 1/c.  No n-by-n GRM other than W exists.

    THE UNIT OF WORK OF THE WHOLE ESTIMATOR.  CG sees V and nothing else, so
    counting these calls (and the columns in them) counts the algorithm; each
    one expands into two K applies and one W apply, which are counted
    separately by the routines below.  See the counter block at the top.
    """
    _count_apply('V', 1 if np.ndim(B) == 1 else np.shape(B)[1])
    return (s2a * compute_KU(Z, B)
            + s2d * compute_KU(Zd, B)
            + s2gxg * compute_WU_dense(W, B)
            + s2e * B)


def _cg_batched(matvec, Bmat, x0=None, tol=1e-6, maxiter=1000):
    """Conjugate gradient for the SPD system V X = Bmat.

    matvec : callable  X -> V X   (accepts / returns (n, c) arrays).
    Bmat   : (n, c) right-hand sides -- all c columns solved together with
             per-column CG scalars, so one V-pass advances every column.
    x0     : (n, c) warm start (e.g. the previous REML iteration's solution).

    Costs 1 + (iterations taken) applies of matvec, every one on all c columns
    at once -- the residual setup below is the +1.  Both the solve and its
    iterations are counted (cg_solves, cg_iters); the applies themselves are
    counted inside _v_matvec, so the two views must agree.
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

    NOT a PSD routine, deliberately.  The rank-r truncated epistasis operator
    W-hat is symmetric but NOT positive semidefinite: dropping the tail of each
    gene's eigendecomposition leaves a Hadamard-square reconstruction that can
    go negative, badly so under linkage equilibrium at small r (measured at
    r = 20, m = 1000, G = 10: lam_min = -5.73 against lam_max = +0.51, so the
    operator is DOMINATED by its negative eigenvalue).  A routine that assumed
    PSD here -- e.g. a plain power iteration clamped at 0 -- would report 0 for
    that end and the feasibility guard built on it would not guard at all.

    Two shifted power iterations on ONE column:

      1. power-iterate on A, giving the dominant-MAGNITUDE eigenvalue mu1 as a
         SIGNED Rayleigh quotient (so mu1 = -5.73 above, not +5.73);
      2. power-iterate on A - mu1 I, whose dominant-magnitude eigenvalue is the
         one FURTHEST from mu1, i.e. the opposite extreme; add mu1 back.

    Returns the sorted pair.  Deterministic (fixed seed, fixed start) and
    computed once per fit, against the Nmc-column solves that dominate every
    REML iteration: unmeasurable.
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
    """Standard errors of the AI-REML estimate, with the Monte-Carlo inflation.

        SE_p = sqrt( [AI^{-1}]_pp ) * sqrt(1 + 1/Nmc)

    The first factor is the usual observed-information standard error: AI is
    the average information, so AI^{-1} is the asymptotic covariance of s-hat.

    THE SECOND FACTOR IS WHAT MONTE-CARLO AI-REML ADDS (BOLT-REML
    supplementary note 2.3).  Exact REML matches the observed data's
    quadratics to their EXPECTED values.  MC REML cannot form those
    expectations, so it matches them instead to an average over Nmc simulated
    reference datasets -- that is exactly what a Hutchinson trace estimate is.
    The reference average carries sampling error of the same kind as the
    observed data's and 1/Nmc of its size, so the estimator behaves like one
    fitting one real dataset against Nmc simulated ones and its variance is
    inflated by (1 + 1/Nmc): 1% at Nmc = 100, hence 0.5% on the SE.

    The correction uses the FINE Nmc alone.  The coarse phase only sets the
    fine phase's starting point; it contributes nothing to the final AI, which
    is computed at the full Nmc by construction (the fit cannot stop in the
    coarse phase -- passing tol_ll_coarse switches it, it does not end it).

    Returns a (k,) array.  This is the ONLY place SEs should come from.
    """
    return np.sqrt(np.diag(np.linalg.inv(AI))) * np.sqrt(1.0 + 1.0 / Nmc)


# ----------------------------------------------------------- MC AI-REML (k=4)
def mc_reml(Z, Zd, W, y, iters=30, Nmc=100, Nmc_coarse=15, cg_tol=1e-6,
            cg_maxiter=1000, jitter=1e-8, tol_ll=1e-4, tol_ll_coarse=1e-2,
            lm=1e-3, eta1=1e-4, eta2=0.99, alpha1=0.25, alpha2=3.5,
            upper_mult=1.5, s_init=(0.05, 0.05, 0.05, 0.85),
            seed=None, verbose=False, w_psd=False):
    """Monte-Carlo AI-REML for V = s2a K_a + s2d K_d + s2gxg W + s2e I (pooled
    UNSTANDARDIZED, C-NORMALIZED W = W_raw / c, PRECOMPUTED as a dense matrix).

    WHAT s2gxg MEANS HERE.  Because the kernel carries the 1/c-hat, s2gxg is
    the REALIZED epistasis variance -- E[Var-hat(g_gxg)] = (c / c-hat) s2gxg --
    and is compared directly with the target and with the replicate's realized
    value.  There is no post-fit c correction to apply, and applying one would
    double-count.  The divisor is baked into the cached W by the Cholesky job
    (pooled_c on the same gene blocks, same C_METHOD), and Simulate_MCREML.py
    refuses a cache whose c-hat differs from the one it recomputes from the
    genotype -- so the plug-in's error is COMMON to both sides and cancels.

    Z is the column-standardized genotype (additive_design), used here for
    K_a = Z Z'/m only; Zd is the standardized dominance design
    (dominance_design), used by K_d only.  Both are applied MATRIX-FREE exactly
    as in the _matfree_ sibling.  The gene split lives entirely inside the
    precomputed W, so this routine takes no G.

    W is the PRECOMPUTED (n, n) epistasis kernel -- THE ONE DIFFERENCE from the
    sibling, which applies W u through the O(n m r) low-rank operator on every
    CG iteration.  Here W u is a single dense gemm, O(n^2) per column, and
    nothing is rebuilt from the genotype.  Which matrix W is was decided when
    it was cached (see the PRECOMPUTED estimation kernel section): the EXACT
    normalized kernel the phenotype was drawn from (r <= 0; no truncation bias,
    PSD), or the sibling's rank-r truncation formed densely (r > 0; the same
    deterministic truncation bias the sibling has, NOT PSD in general).

    w_psd says whether W is PSD BY CONSTRUCTION.  True (the exact kernel) pins
    lam_min(W) = 0 in the feasibility bound, exactly as K_a and K_d are pinned;
    False estimates both ends of W's spectrum, which is what the truncated
    kernel needs.  Passing False for a PSD W is safe, just less tight.

    K_a and K_d are EXACT: they are applied as written, with no truncation, so
    every discrepancy against the truth in s2a-hat or s2d-hat is sampling or
    coupling, never operator error.  verify_preW.py checks this routine against
    the _matfree_ sibling's mc_reml on the same data.

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

    TWO-PHASE MONTE-CARLO SCHEDULE (BOLT-REML supplementary note 3.5).  The
    probe solves are the whole cost of an iteration, and the early iterations
    -- which merely walk the iterate into the neighbourhood of the optimum --
    do not need an accurate trace to do it.  So the fit runs in two phases:

      COARSE  Nmc_coarse probes, stopped at dLL_pred < tol_ll_coarse (1e-2,
              i.e. within ~0.1 SE of the NOISY objective's optimum);
      FINE    all Nmc probes, started from the coarse iterate and stopped at
              the usual dLL_pred < tol_ll.

    THE PROBES ARE DRAWN ONCE AT THE FULL Nmc AND THE COARSE PHASE USES THE
    FIRST Nmc_coarse COLUMNS.  This is not the same as drawing Nmc_coarse
    probes and redrawing Nmc of them at the switch, and the difference is the
    point:

      - K_a U, K_d U and W U are formed ONCE, at full width.  A redraw would
        pay that setup twice, and the W product is the expensive one.
      - The switch EXTENDS the objective rather than replacing it: the fine
        trace is the coarse trace's Nmc_coarse terms plus the rest, with the
        first Nmc_coarse contributions bit-for-bit unchanged.  A redraw moves
        the fixed point discontinuously, and the fine phase would be walking
        toward a different optimum than the one it just stopped at.
      - Because the objective moves continuously, THE TRUST-REGION STATE
        CARRIES OVER: Delta is NOT reset at the switch, and the whole
        iteration stays a deterministic function of s given (seed, Nmc).

    The switch does NOT take the step that triggered it.  That dLL_pred was
    measured on the coarse objective -- a noisier surface, whose optimum is
    not the fine one's -- so the iteration recomputes score and AI at the full
    Nmc and re-tests before moving.  It costs one extra pair of solve groups
    and buys a first fine step chosen by the objective actually being
    optimized.

    Set Nmc_coarse=None (or >= Nmc) for the single-phase fit.

    Adaptive trust region
    ---------------------
    The step is an AI-Newton step confined to a trust region whose radius
    ADAPTS to how well the quadratic model has been predicting (BOLT-REML
    supplementary note 3.2 and 3.4).  This replaces the fixed step_frac * var(y)
    cap and the oscillation damping that used to sit beside it; both are gone.

    Per iteration:

        p          Newton step, constrained to ||diag(AI) * p|| <= Delta
        dLL_pred   p'g - 0.5 p'AI p          (the model's predicted gain)
        dLL_approx 0.5 p'(g(s) + g(s+p))      (the TRAPEZOID rule)
        rho        dLL_approx / dLL_pred
        accept if rho > eta1; shrink Delta if rho < eta1, grow it if rho > eta2

    WHY THE TRAPEZOID.  A trust region needs the ACTUAL gain to compare
    against the predicted one, and the actual gain is a difference of REML log
    likelihoods -- which this estimator cannot compute: log|V| is exactly the
    quantity the whole matrix-free construction avoids.  But the likelihood is
    the line integral of its own gradient along the step, and the gradient is
    the score, which IS computable.  Evaluating that integral by the trapezoid
    rule gives dLL_approx from the two endpoint scores and nothing else.  It is
    exact for a quadratic likelihood and the error is the third derivative's,
    which is the same order as the quadratic model's own error -- so it is
    accurate exactly where rho needs to be trusted.

    WHY THE SCALED NORM.  The radius binds ||diag(AI) * p||, not ||p||.  The
    four components sit on wildly different curvature scales here -- the
    epistasis direction is the flattest by a wide margin, and s2e the
    stiffest -- so a radius on the raw norm would be a loose constraint on the
    stiff directions and a crushing one on the flat direction that actually
    needs to move.  Scaling by diag(AI) puts the constraint in units of
    predicted likelihood change, where the components are comparable.

    COST.  One extra score evaluation per iteration -- one extra batch of
    n_probe CG solves.  But on an ACCEPTED step that score is the next
    iteration's current score, so it is reused rather than recomputed: in
    steady state the fit still pays one score per iteration and only REJECTED
    steps cost extra.  (The AI at the trial point is computed with it, since
    an accepted step needs it anyway; on a rejection those 4 extra right-hand
    sides are wasted, against n_probe probe solves -- 4% at Nmc = 100.)

    THE CONSTRAINED STEP IS AN APPROXIMATION.  BOLT solves the trust-region
    subproblem properly (three NLopt solvers, best of).  Here the problem is
    4-dimensional and the code instead takes the unconstrained Newton step and,
    if it violates the radius, SCALES IT BACK along the same direction.  That
    is not the constrained optimum: the true solution of
    min -p'g + 0.5 p'AI p  s.t. ||Dp|| <= Delta  is
    p(mu) = (AI + mu D'D)^{-1} g, which ROTATES toward the gradient as the
    radius tightens, and rotating is most of what a trust region buys on an
    ill-conditioned AI.  Scaling back keeps the Newton DIRECTION and only
    shortens it, so a binding radius here gives a worse step than BOLT's would.
    It is still a descent direction and still bounded, so the accept/reject
    machinery remains valid; it just needs more iterations when the radius
    binds.  With Delta = Inf initially and the radius shrinking only after a
    rejection, it does not bind at all on a clean replicate.

    Staying inside the positive-definite cone
    -----------------------------------------
    s_lower lets the three GENETIC components go negative (see the bounds
    block below), and a sufficiently negative one takes

        V = s2a K_a + s2d K_d + s2gxg W + s2e I

    out of the positive-definite cone.  Everything downstream then breaks at
    once and silently: CG is not a valid solver for an indefinite system, so
    x = V^{-1} y and the probe solves are meaningless, AI = 0.5 KX' V^{-1} KX
    loses its PSD guarantee and acquires a LARGE negative eigenvalue (observed:
    -9100 against positive ones of 2000-4000), and the Newton step computed
    from it is not an ascent direction -- dLL_pred comes out NEGATIVE.  The
    ridge cannot repair this: it is scaled to mean|diag(AI)|, so it is O(1)
    against an O(1000) negative eigenvalue.

    This bit in practice on NULL runs (s2gxg = 0) at n = m = 1000, where the
    likelihood peaks at or below a boundary and the iterate is pushed negative.

    IT IS NOT ONLY NEGATIVE COMPONENTS THAT DO THIS -- when W is the r > 0
    truncation.  (The exact precomputed kernel is PSD, and with w_psd=True it
    is guarded exactly like K_a and K_d.)  K_a and K_d are Gram
    matrices and exactly PSD, but W-hat is the rank-r TRUNCATION of the pooled
    kernel and is NOT: under linkage equilibrium at small r its most negative
    eigenvalue dominates (measured at r = 20, m = 1000, G = 10: lam_min =
    -5.73 against lam_max = +0.51).  So a POSITIVE s2gxg takes V indefinite
    just as readily as a negative s2a does, and a guard written on the
    assumption that every kernel is PSD does not guard the epistasis direction
    at all.  Observed: an iterate with all four components POSITIVE at
    lam_min(V) = -0.208.

    The fit therefore enforces FEASIBILITY directly, on the two-sided bound

        lam_min(V) >= s2e + sum_i [ s_i lam_min(K_i) if s_i > 0
                                    else s_i lam_max(K_i) ] ,

    with each kernel's spectral range found once by shifted power iteration
    (genotype-only, one column, unmeasurable against the Nmc-column solves).
    A trial point is accepted into the iteration only if that bound is
    comfortably positive.  A step that would leave the cone is BACKTRACKED
    along its own direction until it does not; if no backtrack works, it is
    rejected and the trust region shrinks, which is precisely the machinery
    already in place for a step the model got wrong.  The iterate can still
    approach the boundary -- the negative room that makes a NULL run two-sided
    is intact -- it simply cannot cross it.

    With feasibility enforced, AI is PSD by construction and AI + ridge is
    strictly PD, so dLL_pred > 0 is an INVARIANT rather than a hope.  The
    assertion on it is kept as a cheap check that the invariant holds.

    ONE INTERACTION TO KNOW ABOUT.  dLL_pred is tested AFTER the radius
    constraint, following BOLT's pseudocode -- so a heavily shrunk Delta makes
    the step small and can drive dLL_pred under tol_ll while the iterate is
    still far from the optimum, stopping the fit early.  On a clean replicate
    Delta stays Inf and the tested step is the raw Newton step, which is where
    the dLL_pred interpretation (within sqrt(tol_ll) SE of the optimum) holds.
    After repeated rejections it does not.

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

    Operator-apply count
    --------------------
    The counters described at the top of the module are ZEROED on entry, so
    get_op_counts() after this call reports exactly what THIS fit did: how
    many times V, K and W were applied, to how many columns each, plus the CG
    and REML iteration counts.  Simulate_MCREML.py writes them out per
    replicate and the combine job averages them, next to the wall-clock time.
    They are deterministic given (genotype, y, r, Nmc, seed, tolerances) --
    unlike the time, which is not.

    Returns
    -------
    s  : (4,) estimated (s2a, s2d, s2gxg, s2e).
    AI : (4, 4) final average-information matrix, at the FULL Nmc.  Standard
         errors come from reml_se(AI, Nmc), which applies the Monte-Carlo
         inflation sqrt(1 + 1/Nmc) -- a raw sqrt(diag(inv(AI))) understates
         them.
    """
    reset_op_counts()                  # counts describe THIS replicate alone
    y = np.asarray(y, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float)
    Zd = np.asarray(Zd, dtype=float)
    n = y.shape[0]
    k = 4

    # --- setup: the kernel is already built; nothing genotype-side to do ---
    W = np.asarray(W, dtype=float)
    if W.shape != (n, n):
        raise ValueError(f"W has shape {W.shape}, expected {(n, n)} to match y.")
    Kaapply = lambda B: compute_KU(Z, B)
    Kdapply = lambda B: compute_KU(Zd, B)
    Wapply = lambda B: compute_WU_dense(W, B)

    # SPECTRAL RANGE of the three kernels, for the feasibility bound on
    # lam_min(V).  Genotype-only, computed ONCE outside the iteration.
    #
    # K_a = Z Z'/m and K_d = Z_d Z_d'/m are Gram matrices and therefore PSD
    # EXACTLY, so their lower end is 0 by construction and is written as 0
    # rather than estimated -- a power iteration would return a tiny +-1e-15
    # and the sign of that is noise.  The EXACT precomputed W is a sum of Gram
    # matrices divided by a positive c-hat, PSD for the same reason, and with
    # w_psd=True its lower end is pinned at 0 too.  The rank-r W-hat is a
    # DIFFERENT matter: symmetric but NOT PSD, and under linkage equilibrium at
    # small r its most negative eigenvalue dominates the spectrum entirely
    # (r = 20, m = 1000, G = 10: -5.73 against +0.51).  Both ends of it are
    # therefore estimated, and the lower one is expected to be negative.
    _, lmax_a = _spectral_range(Kaapply, n)
    _, lmax_d = _spectral_range(Kdapply, n)
    lmin_w, lmax_w = _spectral_range(Wapply, n)
    if w_psd:
        lmin_w = 0.0
    lmin_k = np.array([0.0, 0.0, lmin_w])       # K_a, K_d PSD exactly
    lmax_k = np.array([lmax_a, lmax_d, lmax_w])

    def lam_min_V_bound(s_vec):
        """Guaranteed lower bound on lam_min(V) at s_vec.  Positive => V is PD.

        lam_min(sum_i s_i K_i + s2e I) >= s2e + sum_i min over the spectrum of
        s_i K_i, and that minimum is s_i lam_min(K_i) when s_i > 0 and
        s_i lam_max(K_i) when s_i < 0.  Using lam_max for BOTH signs -- i.e.
        assuming every kernel PSD, so that only a negative component can lower
        lam_min -- is wrong here and silently disables the guard on the
        epistasis direction, which is the one that needs it: a POSITIVE s2gxg
        times a negative lam_min(W-hat) takes V indefinite just as readily.
        """
        g = np.asarray(s_vec[:3], dtype=float)
        worst = np.where(g > 0.0, g * lmin_k, g * lmax_k)
        return float(s_vec[3] + worst.sum())

    vary = y.var()
    # UPPER bound, 1.5 * var(y) (was 5.0).  The four components decompose the
    # phenotypic variance, so in expectation they SUM to var(y) and no single
    # one should approach it, let alone exceed it: an s2e above var(y) is not a
    # large estimate, it is a diverged one.  5 * var(y) was loose enough that a
    # run could settle at such a point and still look like a converged fit.
    s_upper = upper_mult * vary

    # --- LOWER bound, PER COMPONENT -- no floor on the three genetic ones ----
    # s2a, s2d and s2gxg are free to go NEGATIVE, but only down to 0.2*s_upper
    # in magnitude (= 0.3*var(y) at upper_mult=1.5), NOT the full -s_upper the
    # symmetric box used to allow; only s2e keeps the 1e-9 floor.  The negative
    # room exists to make a NULL run two-sided, and that needs a band around
    # zero of the order of the estimate's own standard error -- roughly 0.1 in
    # these runs -- not a range wider than the total phenotypic variance.  The
    # extra room bought nothing statistically and cost positive-definiteness.
    # A floor at ~0 makes
    # the estimator one-sided exactly where a NULL run puts the truth: roughly
    # half the sampling distribution is folded onto the floor, so the null mean
    # is biased up by construction and its spread is not the estimator's.  A
    # negative estimate is read as "0, plus the noise around it".
    #
    # s2e KEEPS ITS FLOOR because V = s2a K_a + s2d K_d + s2gxg W + s2e I needs
    # something to hold it inside the positive-definite cone: the three kernels
    # are PSD and singular, so with all four free nothing does.  That floor is
    # NECESSARY, not SUFFICIENT -- a sufficiently negative genetic component
    # still takes V indefinite, and lam_max(K_a) is the largest of the three,
    # so s2a is the one that does it soonest.
    s_lower = np.array([-0.2 * s_upper, -0.2 * s_upper, -0.2 * s_upper, 1e-9])

    rng = np.random.default_rng(seed)
    # Drawn ONCE at the FULL width.  The coarse phase reads U[:, :Nmc_coarse];
    # see the two-phase note in the docstring for why this is not the same as
    # drawing Nmc_coarse probes and redrawing at the switch.
    U = rng.choice([-1.0, 1.0], size=(n, Nmc))      # Rademacher probes in {+-1}
    if Nmc_coarse is None or Nmc_coarse >= Nmc:
        n_probe, coarse = Nmc, False                # single-phase fit
    else:
        n_probe, coarse = Nmc_coarse, True
    switch_it = None                                # iteration the switch fired

    # --- starting point, as FRACTIONS OF var(y) -----------------------------
    # s_init is read as a split of the total phenotypic variance, so the start
    # is vary * s_init and the estimator stays scale-equivariant: rescaling y
    # by a rescales every component by a^2 and the iteration is the same one.
    # The default (0.05, 0.05, 0.05, 0.85) starts the three genetic components
    # SMALL and puts the rest in the residual, rather than the equal split
    # vary/k = 0.25 each this used to use.  Two consequences worth knowing:
    #
    #   - It starts much closer to the truth for the runs this pipeline does
    #     (S2A = S2D = S2GXG = 0.1, S2E = 0.7), and closer still under a null,
    #     so fewer AI-Newton steps are spent walking the residual down.
    #   - It is the ONE place this variant's optimizer now differs from its
    #     siblings', which still start at vary/k.  A run here is therefore not
    #     step-for-step comparable to theirs even on identical data.  The fixed
    #     point is unchanged -- this moves where the iteration starts, not what
    #     it converges to -- but a replicate that ends on a bound, or one that
    #     stops early on the dLL_pred test, can land differently.
    #
    # Both bounds still apply to the start: every entry must lie in
    # [1e-9, upper_mult * vary], which any fraction in (0, upper_mult] does.
    s = vary * np.asarray(s_init, dtype=float)
    if s.shape != (k,):
        raise ValueError(f"s_init must have {k} entries, got {s.shape}")
    AI = np.eye(k)

    yc = y.reshape(n, 1)
    xbuf = None       # warm-start buffers for the CG solve groups
    # Pbuf is held at the FULL Nmc so the coarse phase warm-starts columns
    # 0..Nmc_coarse-1 and the fine phase inherits them; the columns the coarse
    # phase never touched stay zero and are cold-started once, at the switch.
    Pbuf = np.zeros((n, Nmc))
    Gbuf = None

    KaU = Kaapply(U)  # K_a U for the fixed probes: genotype-only, formed once
    KdU = Kdapply(U)  # K_d U   "        "
    WU = Wapply(U)    # W U     "        "

    def eval_grad_ai(s_at):
        """Score and average information at s_at, at the CURRENT n_probe.

        Lifted out of the loop because the trust region needs the score at a
        TRIAL point as well as at the current one, and the two must be the same
        function of s -- same probes, same tolerances -- or rho compares two
        different objectives.  Reads n_probe and the warm-start buffers from
        the enclosing scope and updates the buffers in place; a trial point's
        solutions are a perfectly good warm start for whatever comes next,
        accepted or not.
        """
        nonlocal xbuf, Gbuf
        s2a_, s2d_, s2gxg_, s2e_ = s_at
        mv = lambda B: _v_matvec(Z, Zd, W, s2a_, s2d_, s2gxg_, s2e_, B)

        # --- x = V^{-1} y ---
        xbuf = _cg_batched(mv, yc, x0=xbuf, tol=cg_tol, maxiter=cg_maxiter)
        x_ = xbuf[:, 0]

        # data quadratics x'K_i x
        Kax_, Kdx_, Wx_ = Kaapply(x_), Kdapply(x_), Wapply(x_)

        # --- score traces tr(V^{-1} K_i): Hutchinson, exact probe solves ---
        Uc_ = U[:, :n_probe]
        Psol = _cg_batched(mv, Uc_, x0=Pbuf[:, :n_probe], tol=cg_tol,
                           maxiter=cg_maxiter)
        Pbuf[:, :n_probe] = Psol
        sc = np.array([
            0.5 * (x_ @ Kax_ - np.mean(np.sum(Psol * KaU[:, :n_probe], axis=0))),
            0.5 * (x_ @ Kdx_ - np.mean(np.sum(Psol * KdU[:, :n_probe], axis=0))),
            0.5 * (x_ @ Wx_ - np.mean(np.sum(Psol * WU[:, :n_probe], axis=0))),
            0.5 * (x_ @ x_ - np.mean(np.sum(Psol * Uc_, axis=0)))])

        # --- average information: A_ij = 0.5 (K_i x)' V^{-1}(K_j x) ---
        KX = np.column_stack([Kax_, Kdx_, Wx_, x_])
        Gbuf = _cg_batched(mv, KX, x0=Gbuf, tol=cg_tol, maxiter=cg_maxiter)
        ai = 0.5 * (KX.T @ Gbuf)
        return sc, 0.5 * (ai + ai.T)

    # Delta = Inf: the first step is the unconstrained Newton step, and the
    # radius only ever comes into existence once a step has been rejected.  A
    # replicate that never overshoots never sees a finite radius and follows
    # exactly the path an untruncated AI-Newton iteration would.
    Delta = np.inf
    score, AI = eval_grad_ai(s)          # the current gradient, carried forward
    n_reject = 0

    for it in range(iters):
        _count_event('reml_iters')

        # --- AI-Newton step inside the adaptive trust region ----------------
        # W can still be nearly collinear with I, K_a with W, and K_d with I
        # (dominance deviations are close to independent noise at low MAF), so
        # AI can be near-singular and an undamped step explodes.  The
        # Levenberg-Marquardt ridge (scaled to AI) handles the singularity; the
        # trust region below handles the overshoot, and the box clamp is the
        # last backstop.
        dA = np.abs(np.diag(AI))
        ridge = lm * (dA.mean() + 1e-12)
        # The ridge must make AI + ridge I strictly POSITIVE DEFINITE, which the
        # LM term alone does not guarantee: it is scaled to mean|diag(AI)| and
        # cannot lift an eigenvalue that is negative by a comparable amount.
        # Under the feasibility guard below AI is PSD and this adds nothing;
        # it is kept so a marginal case degrades into a small step rather than
        # a wrong direction.  k = 4, so the eigendecomposition is free.
        lam_lo = float(np.linalg.eigvalsh(AI)[0])
        if lam_lo < 0.0:
            ridge += abs(lam_lo) * (1.0 + 1e-6)
        step = np.linalg.solve(AI + (ridge + jitter) * np.eye(k), score)

        # Radius constraint, on the SCALED norm ||diag(AI) * p|| (see the
        # docstring): scale the Newton step back along its own direction if it
        # violates it.  An APPROXIMATION to the constrained optimum, which
        # would rotate the direction toward the gradient rather than only
        # shorten it -- flagged in the docstring, and inactive while Delta is
        # Inf, which is the whole of a clean replicate.
        dnorm = np.linalg.norm(np.diag(AI) * step)
        if dnorm > Delta and dnorm > 0.0:
            step *= Delta / dnorm
            dnorm = Delta

        # --- CONVERGENCE: BOLT-REML predicted log-likelihood gain ------------
        # Loh et al. 2015 Nat Genet, supplementary note 3.3-3.4.  The local
        # quadratic model around s is
        #     l(s + d) ~= l(s) + score'd - 0.5 d' AI d,
        # so dLL_pred is the model's own remaining height above the iterate.
        #
        # Interpretation: AI is the observed information, so dLL_pred is
        # 0.5 * sum_p ((s_opt - s)_p / SE_p)^2 in the model's metric, and
        # dLL_pred < tol_ll puts the iterate within ~sqrt(tol_ll) standard
        # errors of the optimum -- 1e-2 SE at the default 1e-4.  That is a
        # statement about the ESTIMATE's accuracy, unlike a max|step| test,
        # which is a statement about step size in raw variance units and so
        # depends on the scale of y.  Two inner products, no inverse.
        #
        # The test fires BEFORE the step is taken, so the final iterate is the
        # one the test passed on and the last iteration does not move.  It is
        # tested on the RADIUS-CONSTRAINED step, following BOLT's pseudocode --
        # see the docstring for when that matters.
        dLL_pred = float(score @ step - 0.5 * step @ (AI @ step))
        if verbose:
            print(f"iter {it:2d}  dLL_pred={dLL_pred:.6e}"
                  f"  phase={'coarse' if coarse else 'fine'}(S={n_probe})"
                  f"  Delta={Delta:.4g}", flush=True)
        # d = M^{-1} score with M = AI + ridge STRICTLY positive definite (the
        # ridge is lifted above -lam_min(AI) above), so
        # dLL_pred = score'M^{-1}score - 0.5 d'AI d > 0, and scaling d back by
        # t in (0, 1] keeps it positive.  This is an INVARIANT, not a hope:
        # the feasibility guard keeps V in the PD cone, which keeps AI PSD.
        # It is checked rather than assumed because the consequence of it
        # failing silently is an estimate that looks like a number.
        assert dLL_pred > -1e-8 * max(1.0, abs(dLL_pred)), (
            f"dLL_pred={dLL_pred:.6e} < 0 at iter {it}: AI + ridge is not "
            f"positive definite despite the feasibility guard "
            f"(eigs(AI)={np.linalg.eigvalsh(AI)}, ridge={ridge:.6g}, "
            f"lam_min(V) bound={lam_min_V_bound(s):.6g})")

        # In the COARSE phase the active tolerance is the loose one, and
        # passing it switches the schedule rather than ending the fit.  The
        # step that triggered the switch is DISCARDED: its dLL_pred was
        # measured on the Nmc_coarse-probe objective, and the fine objective's
        # score and AI at this same iterate are different numbers.  Delta
        # carries over untouched -- the radius is a property of how well the
        # QUADRATIC MODEL has been predicting, which the probe count does not
        # change.
        if coarse:
            if dLL_pred < tol_ll_coarse:
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
        # The box clamp is applied to the TRIAL point, so the move actually
        # evaluated is s_try - s and not necessarily step.  rho and the radius
        # update are computed on that actual move: a clamped step that the
        # model did not predict is exactly the case the trust region exists to
        # notice.
        # FEASIBILITY: never hand CG a V that is not positive definite.  The
        # box clamp alone does not prevent it -- s_lower deliberately allows
        # negative genetic components -- so the step is backtracked along its
        # own direction until the bound on lam_min(V) is comfortably positive.
        # The margin is relative to s2e, the only component guaranteed
        # positive, so it scales with the data like everything else here.
        s_try = np.clip(s + step, s_lower, s_upper)
        shrink = 0
        while (lam_min_V_bound(s_try) <= 1e-4 * max(s_try[3], 1e-12)
               and shrink < 40):
            step = 0.5 * step
            s_try = np.clip(s + step, s_lower, s_upper)
            shrink += 1
        if lam_min_V_bound(s_try) <= 1e-4 * max(s_try[3], 1e-12):
            # Even a vanishing step is infeasible, which means the CURRENT
            # iterate is on the boundary.  Reject, shrink the radius hard and
            # let the next iteration try a different direction from the same
            # (still feasible) point.
            n_reject += 1
            Delta = alpha1 * np.linalg.norm(np.diag(AI) * step)
            if verbose:
                print(f"iter {it:2d}  INFEASIBLE even at 2^-{shrink} of the "
                      f"step (lam_min(V) bound <= 0); rejected, Delta->{Delta:.4g}",
                      flush=True)
            continue
        if shrink and verbose:
            print(f"iter {it:2d}  step backtracked 2^-{shrink} to keep V "
                  f"positive definite (lam_min(V) bound="
                  f"{lam_min_V_bound(s_try):.4g})", flush=True)

        taken = s_try - s
        score_try, AI_try = eval_grad_ai(s_try)

        # Model breakdown guard: if the gradient GREW by more than 2x, the
        # quadratic model is not describing this region at all and the
        # trapezoid gain is meaningless too -- reject outright rather than
        # divide two unreliable numbers.
        gnorm, gnorm_try = np.linalg.norm(score), np.linalg.norm(score_try)
        if gnorm_try > 2.0 * gnorm:
            rho = -1.0
        else:
            # TRAPEZOID rule for the actual gain: l(s+p) - l(s) is the line
            # integral of the score along the step, and the score is the one
            # thing here that IS computable (the likelihood itself needs
            # log|V|, which this estimator never forms).
            dLL_approx = float(taken @ (0.5 * (score + score_try)))
            rho = dLL_approx / dLL_pred if dLL_pred != 0.0 else -1.0

        taken_dnorm = np.linalg.norm(np.diag(AI) * taken)
        if rho > eta1:
            # Accept.  The trial score/AI become the current ones -- this is
            # what keeps the steady-state cost at ONE score evaluation per
            # iteration despite the extra gradient.
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

    if verbose:
        print(f"switch at iter {switch_it}; rejected {n_reject} step(s); "
              f"final Delta={Delta:.4g}; SE={reml_se(AI, Nmc)}", flush=True)
    return s, AI


def MC_REML(Z, Zd, W, y, iters=30, Nmc=100, Nmc_coarse=15, cg_tol=1e-6,
            cg_maxiter=1000, seed=None, verbose=False, w_psd=False):
    """Wrapper: returns (s2a_hat, s2d_hat, s2gxg_hat, s2e_hat, AI).

    W is the precomputed dense epistasis kernel (load_W_cache); w_psd comes
    from its sidecar.  SEs from reml_se(AI, Nmc).  Nmc_coarse=None gives the
    single-phase fit.
    """
    s, AI = mc_reml(Z, Zd, W, y, iters=iters, Nmc=Nmc, Nmc_coarse=Nmc_coarse,
                    cg_tol=cg_tol, cg_maxiter=cg_maxiter, seed=seed,
                    verbose=verbose, w_psd=w_psd)
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
