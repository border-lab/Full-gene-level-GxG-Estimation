# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.linalg import cholesky, cho_solve, LinAlgError
from scipy.linalg.lapack import dpotri
import time

####################################################################
# FOUR variance components: additive + dominance + pooled within-gene
# epistasis + noise, fitted by EXACT AI-REML.  The epistasis kernel is
# C-NORMALIZED and PRECOMPUTED as a dense n-by-n matrix by the Cholesky job.
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
# EXACT AI-REML -- WHERE THIS DIRECTORY DIFFERS FROM
# Simulation_code_MCREML_pooled_preW_Lowrank_Wu_std_4VC (the "MC sibling").
# The simulation side -- designs, kernels, c-normalization, Cholesky factors,
# the force_realized phenotype, the W cache -- is the sibling's code, and so is
# the optimizer's OUTER loop: starting point, box, feasibility guard, LM ridge,
# adaptive trust region with the trapezoid gain, and the dLL_pred stopping
# rule.  What changed is how ONE score / AI evaluation is computed:
#
#   MC (sibling)  V applied matrix-free inside CG; tr(V^{-1} K_i) a HUTCHINSON
#                 estimate over Nmc Rademacher probe solves (coarse S = 15
#                 phase, then fine S = Nmc); SEs inflated by sqrt(1 + 1/Nmc).
#   EXACT (here)  V formed densely and Cholesky-factorized, V = L L'.
#                 x = V^{-1} y by two triangular solves, V^{-1} from the factor
#                 (LAPACK potri), and every trace EXACT:
#                     tr(V^{-1} K_i) = <V^{-1}, K_i>_F  (both symmetric).
#                 No probes, no CG tolerance, no seed, no two-phase schedule:
#                 the fit is a deterministic function of (genotype, y).
#
# Consequences worth knowing:
#
#   - W is always the EXACT kernel the phenotype was drawn from.  The
#     sibling's r > 0 option (a dense rank-r W-hat) is not carried over: exact
#     traces on a truncated kernel would put a deterministic truncation bias
#     back in through the other door.  W is PSD by construction.
#   - COST.  O(n^3) per score evaluation (potrf + potri, ~n^3 flops, plus a
#     few O(n^2) ddots and one n-by-4 gemm for AI), and O(n^2 m) once per fit
#     to form K_a and K_d densely.  MEMORY ~7 n-by-n arrays (K_a, K_d, W, V/L,
#     V^{-1} and transients), 56 n^2 bytes: 0.45 GB at n = 8000, 14 GB at
#     n = 16000.
#   - The log-likelihood is now COMPUTABLE (log|V| = 2 sum log diag L).  The
#     verbose trace reports the exact gain next to the trapezoid one, but the
#     trust region still accepts / rejects on the trapezoid rho, so the
#     iteration logic is the sibling's line for line.
#   - A V that leaves the PD cone makes the Cholesky FAIL LOUDLY instead of
#     silently corrupting CG; ai_reml turns that into a rejected step.
#   - With no Monte-Carlo reference sample the SE is the plain
#     sqrt(diag(AI^{-1})) -- no (1 + 1/Nmc) inflation.
#
# NORMALIZATION.  The raw pooled kernel makes the REML component s2gxg and the
# REALIZED variance of the epistasis genetic value two different numbers:
#
#     E[ Var-hat(H gamma) ] = c * s2gxg ,
#     c = (1/P) sum_g sum_{a<b in g} Var-hat(Z_a .* Z_b) ,
#
# with c != 1 because the interaction columns h_ab are not standardized.  The
# gap is closed BEFORE anything is drawn, by dividing the kernel itself once:
#
#     W = W_raw / c-hat   ==>   E[ Var-hat(H gamma) ] = (c / c-hat) s2gxg .
#
# So with s2a = s2d = s2gxg = 0.1 and s2e = 0.7 all four components have their
# realized variance centred on the value that was asked for, and no post-fit
# correction is applied to anything:  c_a = c_d = 1 already (both designs are
# column-standardized), c_e = 1 trivially, and c_gxg = c / c-hat.
#
# WHICH c-hat -- THE THREE ROUTES.  c is a GENOTYPE-ONLY constant and this
# module carries three ways to get it (realized_variance.pdf).  All three share
# the identity Var(Z_a Z_b) = 1 + r_ab s_a s_b for the per-pair variance; they
# differ in where the per-SNP skewness s comes from, and 'exact' skips the
# identity altogether:
#
#   'moment'  (C_METHOD -- THE ONE THE PIPELINE USES.)  s-hat_a = mean_t
#             (Z_ta^3), the EMPIRICAL third moment of the standardized column.
#             O(nm) time, no HWE assumption; its only error is O(1/sqrt n)
#             moment sampling.
#   'hwe'     the note's closed form, s_a = (1 - 2 p_a) / sqrt(2 p_a (1 - p_a))
#             from the allele frequency -- the skewness of Z_a ONLY under HWE.
#   'exact'   the literal mean per-pair sample variance, O(n sum_g m_g^2).
#             The YARDSTICK the other two are scored against; reported next
#             to every run, never used as the divisor.
#
# BOTH SIDES DIVIDE BY THE SAME c-hat.  The Cholesky job factorizes s2gxg W
# and caches that same W; the estimator recomputes c-hat from the genotype
# (O(nm)) and refuses a cached kernel whose recorded divisor differs.  The
# residual c / c-hat gap is written to result/c_<FILENAME>.txt, and
# force_realized absorbs it at the draw.
#
# WHAT NORMALIZING DOES NOT FIX.  tr(W) is still not n and W 1 is still not 0:
# dividing by a scalar cannot centre a kernel or restore a trace identity.
#
# THE EPISTASIS KERNEL IS BUILT FROM Z_a ONLY.  h_ab = Z_a .* Z_b pairs
# ADDITIVE columns; dominance enters through K_d and through nothing else.
#
# IDENTIFIABILITY.  K_a, K_d and W are far less collinear with each other than
# W is with I, but the AI matrix has a 4-dimensional curvature to resolve, and
# small n / large m runs can be weakly identified.  K_d is the weakest of the
# three in practice.  Exact traces remove the Monte-Carlo noise from the
# objective, NOT this: a flat likelihood stays flat.  The optimizer's LM ridge
# + trust region + box clamp are unchanged and carry the same known boundary
# defect (see ai_reml).
####################################################################


# ------------------------------------------------------------ cost counting
# HOW MUCH LINEAR ALGEBRA ONE REPLICATE ACTUALLY DOES, machine-independently.
# The MC sibling counted operator APPLIES (V / K / W columns, CG iterations),
# the cost unit of a matrix-free fit.  Exact AI-REML does no applies worth
# counting: its whole per-iteration cost is one dense Cholesky of V (plus the
# inverse from it) per score / AI evaluation.  So the counters record that:
#
#   reml_iters        AI-REML iterations taken (<= iters: convergence break).
#   V_factorizations  score / AI evaluations, each one potrf + potri, O(n^3).
#                     = 1 (the start) + one per trial point actually evaluated.
#   rejected_steps    trial steps turned down -- by the trust region, by the
#                     feasibility guard, or by a failed Cholesky.  Those that
#                     reached a factorization bought it for nothing.
#
# The counters are GLOBAL and CUMULATIVE; ai_reml zeroes them on entry, so
# after one AI_REML call get_op_counts() describes exactly one replicate.  The
# per-rep file format is the sibling's, so combine_code.sh's key-agnostic awk
# reduction is unchanged.
_OP_COUNTS = {}


def reset_op_counts():
    """Zero every counter.  ai_reml calls this on entry."""
    _OP_COUNTS.clear()


def get_op_counts():
    """Snapshot of the counters as a plain dict, safe to keep or write out.

    Read it AFTER an ai_reml / AI_REML call and it describes that call alone.
    Absent keys mean zero: nothing pre-registers a counter it never bumps.
    """
    return dict(_OP_COUNTS)


def _count_event(name, k=1):
    """Record k occurrences of a counted event."""
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


# --------------------------------------------- standardized GRMs (dense)
# K_a and K_d have IDENTICAL structure -- Z_i Z_i' / m for a column-standardized
# design Z_i -- and differ only in which design is passed in.  The MC sibling
# applied them matrix-free inside CG; exact REML needs V as an explicit matrix
# to factorize, so here the simulation AND the estimator form them densely,
# once each, from the same design -- K_a and K_d are identical on both sides.
def build_K(Zi):
    """Explicit GRM  K_i = Z_i Z_i' / m  (n-by-n), O(n^2 m).

    tr(K_i) = n EXACTLY for either design (every column has sample variance 1,
    so sum_t sum_a Z_ta^2 / m = n) -- the identity the unstandardized epistasis
    kernel gives up but the standardized ones keep.
    """
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

    The identity is exact: it agrees with the literal C(m_g, 2)-term pair sum
    to machine precision (see the _matfree_ sibling's verify_lowrank.py).
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
    build_W_pooled, so c averages over the same P pairs the
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
    contributes no pair and is skipped, exactly as in build_W_pooled
    and pooled_c_exact, so the average runs over the same P pairs
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

    all documented in the module header.  build_W_pooled (simulation) calls
    this with the default; the estimator does not rebuild W but RECOMPUTES
    c-hat here from the same genotype by the same route, and load_W_cache
    refuses a cached kernel whose recorded divisor differs.

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
    time, plus the O(nm) of c-hat, which is negligible beside it.  Built ONCE
    by the SIMULATION (the Cholesky factor needs an explicit matrix), which
    also caches it as the PRECOMPUTED kernel ai_reml fits.

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
# The Cholesky job forms the EXACT c-normalized W once -- it needs it anyway
# for Lgxg -- writes it to this pipeline's own W/ directory before freeing it,
# and every AI-REML replicate loads it.  Simulation and estimation therefore
# use literally the same matrix: no rebuild, no truncation, no drift.
#
# THE CACHE IS THIS PIPELINE'S OWN.  It lives under <DIR>/W/ and its name
# carries (mode, n, m, G), never the shared stored_genotype/W_*.npy that the
# dense gxg family writes if-not-exists and keys by genotype only -- a kernel
# from another pipeline must never be picked up silently.  It is OVERWRITTEN
# (atomically) by every Cholesky job rather than written if-not-exists, so the
# W a run fits is always the W its own Cholesky job built.  A JSON sidecar
# records (kernel, mode, n, m, G, c-hat, psd) and the estimator refuses a
# kernel whose sidecar disagrees with its arguments, or whose c-hat differs
# from the one it recomputes from the genotype.
W_CACHE_KERNEL = "pooled_unstd_cnorm_4VC"   # written into the sidecar


def w_cache_paths(root, mode, n, m, G):
    """(npy, json) paths of the cached estimation kernel for one genotype and
    gene split."""
    base = f"{root}/W_{mode}_n{n}_m{m}_G{G}"
    return base + ".npy", base + ".json"


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


def load_W_cache(npy_path, json_path, mode, n, m, G, c_expected):
    """Load a cached estimation kernel, refusing one that is not THIS run's.

    Checks the sidecar against (kernel, mode, n, m, G), the array shape against
    (n, n), and the recorded c-hat against c_expected -- the divisor the caller
    recomputed from the genotype it is about to fit with (O(nm), one mat-vec
    per gene).  The last check is what catches a kernel left over from a
    genotype file that has since been regenerated under the same name.

    Returns (W, meta).
    """
    import json
    import os
    for p in (npy_path, json_path):
        if not os.path.exists(p):
            raise FileNotFoundError(
                f"precomputed kernel file not found: {p}\n"
                f"Run the Cholesky step with the same (mode, n, m, G) first.")
    with open(json_path) as f:
        meta = json.load(f)
    want = {'kernel': W_CACHE_KERNEL, 'mode': mode, 'n': int(n), 'm': int(m),
            'G': int(G)}
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


# ------------------------------------------------------------- simulation
def simulate_Cholesky_4vc(real_data, G, s2a=0.1, s2d=0.1, s2gxg=0.1, s2e=0.7,
                          stability=1e-10, save_W=None):
    """Cholesky factors of the additive, dominance and pooled-epistasis
    covariances.

        La La'     = s2a   K_a ,   K_a = Z_a Z_a' / m ,
        Ld Ld'     = s2d   K_d ,   K_d = Z_d Z_d' / m ,
        Lgxg Lgxg' = s2gxg W ,     W   = W_raw / c .

    The epistasis factor is built from the C-NORMALIZED W, so the g_gxg it
    draws has E[Var-hat(g_gxg)] = (c / c-hat) s2gxg -- the same target the
    additive and dominance factors already hit, their designs being
    standardized.

    PRECOMPUTING THE ESTIMATION KERNEL.  save_W, if given, is called ONCE as
    save_W(W, c, build_time) with the very array Lgxg was factorized from,
    handed over before it is freed -- so there is no second build and no chance
    of the simulated and the fitted kernel drifting apart.  The callback is
    expected to write it to disk; this function keeps no reference to it.

    Returns (La, Ld, Lgxg, w_build_time, c), w_build_time being the wall-clock
    cost of the EPISTASIS kernel alone -- the O(n^2 m) Hadamard square that
    dominates this job.  K_a and K_d cost O(n^2 m) too but with a single gemm
    each, so they are folded into the untimed remainder.

    The three factors are built SEQUENTIALLY and each dense GRM is released as
    soon as its factor exists, so the peak is ~4 n-by-n arrays (two finished
    factors + one GRM + the factorization workspace) rather than 7.
    """
    Za = additive_design(real_data)
    Zd = dominance_design(real_data)
    n, m = Za.shape

    # --- epistasis first (the expensive one), then free W -------------------
    genes = split_into_genes(Za, G)          # several Z, one per gene

    # The c-normalization happens INSIDE build_W_pooled and is inside the
    # timed block: it is part of building this kernel.
    t_start = time.perf_counter()
    W, c = build_W_pooled(genes, return_c=True)
    w_build_time = time.perf_counter() - t_start

    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)

    # --- hand the SAME W to the cache, then free it -------------------------
    if save_W is not None:
        save_W(W, c, w_build_time)
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
    make up H, exactly the pairs build_W_pooled counts in P.

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
def _spectral_range(apply, n, iters=80, tol=1e-7, seed=0):
    """(lam_min, lam_max) of a SYMMETRIC operator given only its apply.

    The MC sibling's routine, unchanged, so the feasibility bound in ai_reml is
    computed exactly as it is there.  Two shifted power iterations on ONE
    column:

      1. power-iterate on A, giving the dominant-MAGNITUDE eigenvalue mu1 as a
         SIGNED Rayleigh quotient;
      2. power-iterate on A - mu1 I, whose dominant-magnitude eigenvalue is the
         one FURTHEST from mu1, i.e. the opposite extreme; add mu1 back.

    Returns the sorted pair.  Deterministic (fixed seed, fixed start).  All
    three kernels here are PSD, so only the upper end is used.  It costs a few
    hundred n-by-n gemvs once per fit -- negligible against one Cholesky of V.

    A power iteration approaches lam_max FROM BELOW, so a bound built on it can
    be marginally optimistic.  The sibling had no way to notice; here a V that
    is not PD makes the Cholesky fail loudly, and ai_reml rejects that step.
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


def _factor_V(Ka, Kd, W, s_at):
    """Cholesky factor, inverse and log-determinant of

        V = s2a K_a + s2d K_d + s2gxg W + s2e I ,   s_at = (s2a, s2d, s2gxg, s2e).

    THE UNIT OF WORK OF THE WHOLE ESTIMATOR, O(n^3): potrf (n^3/3 flops) and
    potri (2 n^3/3).  Counted as one V_factorization.

    Raises LinAlgError if V is not positive definite -- the loud failure the
    sibling's CG could not give.

    Returns (L, Vinv, logdet): L lower triangular, Vinv the FULL symmetric
    inverse, logdet = log|V| = 2 sum_t log L_tt.
    """
    s2a, s2d, s2gxg, s2e = (float(v) for v in s_at)
    n = W.shape[0]
    V = s2a * Ka
    V += s2d * Kd
    V += s2gxg * W
    V.flat[::n + 1] += s2e
    _count_event('V_factorizations')
    L = cholesky(V, lower=True, overwrite_a=True, check_finite=False)
    logdet = 2.0 * float(np.sum(np.log(np.diag(L))))
    Vinv, info = dpotri(L, lower=1)
    if info != 0:
        raise LinAlgError(f"potri failed (info={info}): V is singular.")
    # potri writes the inverse into the LOWER triangle only; rebuild the full
    # symmetric matrix from that half (np.tril copies, so no aliasing).
    Vinv = np.tril(Vinv)
    Vinv += np.tril(Vinv, -1).T
    return L, Vinv, logdet


def _score_ai(Ka, Kd, W, y, s_at):
    """EXACT score, average information and log-likelihood at s_at.

        x        = V^{-1} y ,
        score_i  = 0.5 ( x' K_i x - tr(V^{-1} K_i) ) ,     K_4 = I ,
        AI_ij    = 0.5 (K_i x)' V^{-1} (K_j x) ,
        logL     = -0.5 ( log|V| + y' V^{-1} y )          (constant dropped).

    The sibling's formulas -- there is no fixed effect, so the REML projection
    P is V^{-1} itself and REML coincides with ML -- with the one Monte-Carlo
    quantity replaced: tr(V^{-1} K_i) is the exact Frobenius inner product
    <V^{-1}, K_i> (both symmetric), ONE ddot over n^2 entries, where the
    sibling had a Hutchinson mean over Nmc probe solves.

    Module-level rather than a closure inside ai_reml so verify_exact.py can
    check it against finite differences of logL.

    Returns (score (4,), AI (4, 4), logL).
    """
    L, Vinv, logdet = _factor_V(Ka, Kd, W, s_at)
    x = cho_solve((L, True), y, check_finite=False)
    del L
    Kax, Kdx, Wx = Ka @ x, Kd @ x, W @ x

    # --- score: exact traces ---
    tr = np.array([np.vdot(Vinv, Ka), np.vdot(Vinv, Kd), np.vdot(Vinv, W),
                   np.trace(Vinv)])
    quad = np.array([x @ Kax, x @ Kdx, x @ Wx, x @ x])
    sc = 0.5 * (quad - tr)

    # --- average information: A_ij = 0.5 (K_i x)' V^{-1} (K_j x) ---
    KX = np.column_stack([Kax, Kdx, Wx, x])
    ai = 0.5 * (KX.T @ (Vinv @ KX))

    ll = -0.5 * (logdet + float(y @ x))
    return sc, 0.5 * (ai + ai.T), ll


def reml_se(AI):
    """Standard errors of the AI-REML estimate:  SE_p = sqrt( [AI^{-1}]_pp ).

    AI is the average information at the optimum, so AI^{-1} is the
    asymptotic covariance of s-hat.  Unlike the MC sibling's reml_se there is
    NO sqrt(1 + 1/Nmc) factor: the traces are exact, so the fit matches the
    data's quadratics to their true expectations rather than to an average
    over simulated reference datasets.

    Returns a (k,) array.  This is the ONLY place SEs should come from.
    """
    return np.sqrt(np.diag(np.linalg.inv(AI)))


# ------------------------------------------------------------ exact AI-REML
def ai_reml(Z, Zd, W, y, iters=30, jitter=1e-8, tol_ll=1e-4,
            lm=1e-3, eta1=1e-4, eta2=0.99, alpha1=0.25, alpha2=3.5,
            upper_mult=1.5, s_init=(0.05, 0.05, 0.05, 0.85),
            verbose=False, w_psd=True):
    """EXACT AI-REML for V = s2a K_a + s2d K_d + s2gxg W + s2e I (pooled
    UNSTANDARDIZED, C-NORMALIZED W = W_raw / c, PRECOMPUTED as a dense matrix).

    WHAT s2gxg MEANS HERE.  Because the kernel carries the 1/c-hat, s2gxg is
    the REALIZED epistasis variance -- E[Var-hat(g_gxg)] = (c / c-hat) s2gxg --
    and is compared directly with the target.  There is no post-fit c
    correction to apply, and applying one would double-count.

    Z is the column-standardized genotype (additive_design), Zd the
    standardized dominance design (dominance_design); K_a and K_d are formed
    DENSELY from them once, on entry.  W is the (n, n) precomputed kernel
    (load_W_cache); the gene split lives entirely inside it, so this routine
    takes no G.  w_psd pins lam_min(W) = 0 in the feasibility bound (True for
    the exact kernel, which is a sum of Gram matrices over a positive c-hat);
    False estimates both ends, which is safe, just less tight.

    What is exact and what is not
    -----------------------------
    Every score / AI evaluation (_score_ai) Cholesky-factorizes V and computes
    x = V^{-1} y, all four traces tr(V^{-1} K_i) and the AI matrix EXACTLY, up
    to floating point.  The OPTIMIZER around it is the MC sibling's mc_reml,
    minus the pieces that only existed because of the Monte Carlo:

      REMOVED  Hutchinson probes, the Nmc / Nmc_coarse two-phase schedule and
               its switch, CG (tolerance, warm starts), the probe seed, the
               sqrt(1 + 1/Nmc) SE inflation.
      KEPT     s_init as fractions of var(y); the per-component box; the
               spectral feasibility guard with backtracking; the LM ridge
               lifted above -lam_min(AI); the scaled-norm trust region with
               the TRAPEZOID gain and its accept / reject / radius rules; the
               gradient-growth guard; the dLL_pred < tol_ll stopping rule.
      ADDED    a failed Cholesky at a trial point is a rejected step (see
               below); the exact log-likelihood is carried along and printed
               in the verbose trace, but NOT used to accept or reject.

    So on the same data the two estimators walk the same kind of path; they
    differ by the Hutchinson noise in the sibling's score, and nothing else.

    Adaptive trust region
    ---------------------
    Per iteration:

        p          Newton step, constrained to ||diag(AI) * p|| <= Delta
        dLL_pred   p'g - 0.5 p'AI p          (the model's predicted gain)
        dLL_approx 0.5 p'(g(s) + g(s+p))      (the TRAPEZOID rule)
        rho        dLL_approx / dLL_pred
        accept if rho > eta1; shrink Delta if rho < eta1, grow it if rho > eta2

    The trapezoid is exact for a quadratic likelihood and its error is the
    third derivative's, the same order as the quadratic model's own.  The
    sibling used it because log|V| was unaffordable; here log|V| is free, but
    rho stays on the trapezoid so that the iteration logic is identical and a
    difference between the two pipelines is the trace estimator and nothing
    else.  (verbose prints dLL_exact next to it -- they agree closely.)

    The radius binds the SCALED norm ||diag(AI) * p||: the four components sit
    on very different curvature scales (epistasis flattest, s2e stiffest), and
    the scaling puts the constraint in units of predicted likelihood change.
    An over-long step is scaled back along the Newton direction, an
    APPROXIMATION to the constrained optimum (which would rotate toward the
    gradient).  With Delta = Inf initially it does not bind on a clean
    replicate.

    Staying inside the positive-definite cone
    -----------------------------------------
    s_lower lets the three GENETIC components go negative (a NULL run needs a
    two-sided estimator), and a sufficiently negative one takes V out of the
    PD cone.  The fit enforces feasibility on the bound

        lam_min(V) >= s2e + sum_i [ s_i lam_min(K_i) if s_i > 0
                                    else s_i lam_max(K_i) ] ,

    each kernel's spectral range found once by shifted power iteration.  A
    trial step that would leave the cone is BACKTRACKED along its own direction
    until the bound is comfortably positive; if no backtrack works it is
    rejected and the trust region shrinks.

    NEW HERE: lam_max from power iteration is approached from below, so the
    bound can be marginally optimistic.  If the Cholesky of V at a trial point
    then fails, the step is rejected exactly like an infeasible one -- where
    the sibling's CG would have returned a silently meaningless solve.

    With feasibility enforced, AI is PSD and AI + ridge strictly PD, so
    dLL_pred > 0 is an INVARIANT; the assertion on it is a cheap check.

    ONE INTERACTION TO KNOW ABOUT.  dLL_pred is tested AFTER the radius
    constraint, following BOLT's pseudocode -- so a heavily shrunk Delta makes
    the step small and can drive dLL_pred under tol_ll while the iterate is
    still far from the optimum, stopping the fit early.  On a clean replicate
    Delta stays Inf and the tested step is the raw Newton step.

    A KNOWN DEFECT OF THE SHARED OPTIMIZER, inherited deliberately.  When a
    replicate's likelihood peaks on a box bound, that component sticks there
    and the others converge to the WRONG values: the AI-Newton step is solved
    jointly and never re-projected onto the free subspace.  Left uncorrected
    so that a run here differs from the MC sibling's in the TRACES alone.
    calc_stats.py counts the affected replicates.

    Cost counting
    -------------
    The counters at the top of the module are ZEROED on entry, so
    get_op_counts() afterwards reports THIS fit: reml_iters, V_factorizations
    (each O(n^3)) and rejected_steps.  They are deterministic given
    (genotype, y, tolerances) -- unlike the wall-clock time.

    Returns
    -------
    s  : (4,) estimated (s2a, s2d, s2gxg, s2e).
    AI : (4, 4) final average-information matrix.  SEs from reml_se(AI).
    """
    reset_op_counts()                  # counts describe THIS replicate alone
    y = np.asarray(y, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float)
    Zd = np.asarray(Zd, dtype=float)
    n = y.shape[0]
    k = 4

    # --- setup: dense K_a, K_d once; W arrives precomputed -------------------
    W = np.asarray(W, dtype=float)
    if W.shape != (n, n):
        raise ValueError(f"W has shape {W.shape}, expected {(n, n)} to match y.")
    if Z.shape[0] != n or Zd.shape[0] != n:
        raise ValueError(f"designs have {Z.shape[0]} / {Zd.shape[0]} rows, "
                         f"expected {n} to match y.")
    Ka = build_K(Z)
    Kd = build_K(Zd)

    # SPECTRAL RANGE of the three kernels, for the feasibility bound on
    # lam_min(V).  Genotype-only, computed ONCE outside the iteration.  K_a and
    # K_d are Gram matrices, PSD exactly, so their lower end is written as 0
    # rather than estimated (a power iteration would return +-1e-15 noise);
    # the exact W is PSD for the same reason and w_psd=True pins it likewise.
    _, lmax_a = _spectral_range(lambda B: Ka @ B, n)
    _, lmax_d = _spectral_range(lambda B: Kd @ B, n)
    lmin_w, lmax_w = _spectral_range(lambda B: W @ B, n)
    if w_psd:
        lmin_w = 0.0
    lmin_k = np.array([0.0, 0.0, lmin_w])       # K_a, K_d PSD exactly
    lmax_k = np.array([lmax_a, lmax_d, lmax_w])

    def lam_min_V_bound(s_vec):
        """Guaranteed lower bound on lam_min(V) at s_vec (up to the power
        iteration's accuracy).  Positive => V is PD."""
        g = np.asarray(s_vec[:3], dtype=float)
        worst = np.where(g > 0.0, g * lmin_k, g * lmax_k)
        return float(s_vec[3] + worst.sum())

    vary = y.var()
    # UPPER bound, 1.5 * var(y).  The four components decompose the phenotypic
    # variance, so in expectation they SUM to var(y) and no single one should
    # approach it: an s2e above var(y) is not a large estimate, it is a
    # diverged one.
    s_upper = upper_mult * vary

    # --- LOWER bound, PER COMPONENT -- no floor on the three genetic ones ----
    # s2a, s2d and s2gxg may go NEGATIVE, down to 0.2*s_upper in magnitude
    # (= 0.3*var(y) at upper_mult=1.5): enough of a band around zero for a NULL
    # run to be two-sided, not a range wider than the phenotypic variance.  A
    # floor at ~0 would fold half the null sampling distribution onto it and
    # bias the null mean up by construction.
    #
    # s2e KEEPS ITS FLOOR because V needs something to hold it inside the PD
    # cone: the three kernels are PSD and singular.  That floor is NECESSARY,
    # not SUFFICIENT -- the feasibility guard does the rest.
    s_lower = np.array([-0.2 * s_upper, -0.2 * s_upper, -0.2 * s_upper, 1e-9])

    # --- starting point, as FRACTIONS OF var(y) -----------------------------
    # vary * s_init keeps the estimator scale-equivariant: rescaling y by a
    # rescales every component by a^2 and the iteration is the same one.  The
    # default starts the genetic components SMALL and the rest in the residual.
    s = vary * np.asarray(s_init, dtype=float)
    if s.shape != (k,):
        raise ValueError(f"s_init must have {k} entries, got {s.shape}")

    # Delta = Inf: the first step is the unconstrained Newton step, and the
    # radius only ever comes into existence once a step has been rejected.
    Delta = np.inf
    score, AI, ll = _score_ai(Ka, Kd, W, y, s)   # carried forward
    n_reject = 0
    if verbose:
        print(f"start s={s}  logL={ll:.6f}", flush=True)

    for it in range(iters):
        _count_event('reml_iters')

        # --- AI-Newton step inside the adaptive trust region ----------------
        # W can be nearly collinear with I, K_a with W, and K_d with I, so AI
        # can be near-singular and an undamped step explodes.  The LM ridge
        # (scaled to AI) handles the singularity; the trust region handles the
        # overshoot, and the box clamp is the last backstop.  The ridge is
        # lifted above -lam_min(AI) so AI + ridge I is strictly PD even in a
        # marginal case.  k = 4, so the eigendecomposition is free.
        dA = np.abs(np.diag(AI))
        ridge = lm * (dA.mean() + 1e-12)
        lam_lo = float(np.linalg.eigvalsh(AI)[0])
        if lam_lo < 0.0:
            ridge += abs(lam_lo) * (1.0 + 1e-6)
        step = np.linalg.solve(AI + (ridge + jitter) * np.eye(k), score)

        # Radius constraint on the SCALED norm; inactive while Delta is Inf.
        dnorm = np.linalg.norm(np.diag(AI) * step)
        if dnorm > Delta and dnorm > 0.0:
            step *= Delta / dnorm
            dnorm = Delta

        # --- CONVERGENCE: BOLT-REML predicted log-likelihood gain ------------
        # dLL_pred = score'd - 0.5 d'AI d is the local quadratic model's own
        # remaining height above the iterate; dLL_pred < tol_ll puts the
        # iterate within ~sqrt(tol_ll) standard errors of the optimum (1e-2 SE
        # at 1e-4).  Tested BEFORE the step is taken, on the radius-constrained
        # step, following BOLT's pseudocode.
        dLL_pred = float(score @ step - 0.5 * step @ (AI @ step))
        if verbose:
            print(f"iter {it:2d}  dLL_pred={dLL_pred:.6e}  Delta={Delta:.4g}",
                  flush=True)
        assert dLL_pred > -1e-8 * max(1.0, abs(dLL_pred)), (
            f"dLL_pred={dLL_pred:.6e} < 0 at iter {it}: AI + ridge is not "
            f"positive definite despite the feasibility guard "
            f"(eigs(AI)={np.linalg.eigvalsh(AI)}, ridge={ridge:.6g}, "
            f"lam_min(V) bound={lam_min_V_bound(s):.6g})")
        if dLL_pred < tol_ll:
            break

        # --- trial point, feasibility backtracking --------------------------
        # The box clamp is applied to the TRIAL point, so the move actually
        # evaluated is s_try - s and not necessarily step; rho and the radius
        # update use that actual move.  The margin is relative to s2e, the only
        # component guaranteed positive, so it scales with the data.
        s_try = np.clip(s + step, s_lower, s_upper)
        shrink = 0
        while (lam_min_V_bound(s_try) <= 1e-4 * max(s_try[3], 1e-12)
               and shrink < 40):
            step = 0.5 * step
            s_try = np.clip(s + step, s_lower, s_upper)
            shrink += 1
        if lam_min_V_bound(s_try) <= 1e-4 * max(s_try[3], 1e-12):
            # Even a vanishing step is infeasible: the CURRENT iterate is on
            # the boundary.  Reject, shrink the radius hard and let the next
            # iteration try a different direction from the same point.
            n_reject += 1
            _count_event('rejected_steps')
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
        taken_dnorm = np.linalg.norm(np.diag(AI) * taken)
        try:
            score_try, AI_try, ll_try = _score_ai(Ka, Kd, W, y, s_try)
        except LinAlgError:
            # The spectral bound passed but V is not PD after all -- the power
            # iteration's lam_max was a hair low.  Same treatment as an
            # infeasible step: reject and shrink onto it.
            n_reject += 1
            _count_event('rejected_steps')
            Delta = alpha1 * taken_dnorm
            if verbose:
                print(f"iter {it:2d}  Cholesky of V FAILED at s_try={s_try}; "
                      f"rejected, Delta->{Delta:.4g}", flush=True)
            continue

        # Model breakdown guard: if the gradient GREW by more than 2x, the
        # quadratic model is not describing this region at all -- reject
        # outright rather than divide two unreliable numbers.
        gnorm, gnorm_try = np.linalg.norm(score), np.linalg.norm(score_try)
        if gnorm_try > 2.0 * gnorm:
            rho = -1.0
        else:
            # TRAPEZOID rule for the actual gain: the line integral of the
            # score along the step (see the docstring for why not ll_try - ll).
            dLL_approx = float(taken @ (0.5 * (score + score_try)))
            rho = dLL_approx / dLL_pred if dLL_pred != 0.0 else -1.0

        dLL_exact = ll_try - ll
        if rho > eta1:
            # Accept.  The trial score/AI/logL become the current ones, so the
            # steady-state cost is ONE factorization per iteration.
            s, score, AI, ll = s_try, score_try, AI_try, ll_try
        else:
            n_reject += 1
            _count_event('rejected_steps')
        if rho < eta1:
            Delta = alpha1 * taken_dnorm          # shrink onto the bad step
        elif rho >= eta2:
            Delta = max(Delta, alpha2 * taken_dnorm)   # model is good: expand

        if verbose:
            print(f"iter {it:2d}  s={s}  rho={rho:+.4f}"
                  f"  dLL_exact={dLL_exact:+.4e}"
                  f"  {'ACCEPT' if rho > eta1 else 'REJECT'}"
                  f"  |D p|={taken_dnorm:.4g}  Delta->{Delta:.4g}", flush=True)

    if verbose:
        print(f"rejected {n_reject} step(s); final Delta={Delta:.4g}; "
              f"logL={ll:.6f}; SE={reml_se(AI)}", flush=True)
    return s, AI


def AI_REML(Z, Zd, W, y, iters=30, verbose=False, w_psd=True):
    """Wrapper: returns (s2a_hat, s2d_hat, s2gxg_hat, s2e_hat, AI).

    W is the precomputed dense epistasis kernel (load_W_cache); w_psd comes
    from its sidecar.  SEs from reml_se(AI).
    """
    s, AI = ai_reml(Z, Zd, W, y, iters=iters, verbose=verbose, w_psd=w_psd)
    s2a_hat, s2d_hat, s2gxg_hat, s2e_hat = s
    return s2a_hat, s2d_hat, s2gxg_hat, s2e_hat, AI
