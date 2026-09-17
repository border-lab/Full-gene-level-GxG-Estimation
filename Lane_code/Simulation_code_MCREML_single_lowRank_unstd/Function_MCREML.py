# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.linalg import cholesky
from scipy.sparse.linalg import svds
import time


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
    """Additive design Z_a: column-standardized allele dosages."""
    return _standardize_cols(real_data)


def split_into_genes(Z, G):
    """Split the m SNP columns of Z into G contiguous gene blocks.

    Returns a list [Z_1, ..., Z_G] of column sub-matrices -- the SEVERAL Z, one
    per gene, that the pooled kernel consumes.  np.array_split handles m not
    divisible by G by making the first (m mod G) genes one SNP larger, so the
    genes may be unequal.
    """
    m = Z.shape[1]
    return [Z[:, cols] for cols in np.array_split(np.arange(m), G)]


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
    time with the default 'hadamard' route.  USED ONLY BY THE SIMULATION (the
    Cholesky factor needs an explicit matrix); estimation applies a rank-r
    TRUNCATION of the same W matrix-free via compute_WU_pooled.

    tr(W) is NOT n here.  With standardized SNPs, W_tt = (1/2P) sum_g
    [(sum_a Z_ta^2)^2 - sum_a Z_ta^4], so tr(W) = n only in expectation under
    independence and fluctuates with the genotype -- the standardized siblings'
    tr(W) = n identity does not survive dropping the 1/sigma_ab weight.
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
def simulate_Cholesky_gxg(real_data, G, s2gxg=0.5, s2e=0.5, stability=1e-10):
    """Cholesky factor of the pooled epistasis covariance.

        Lgxg Lgxg' = s2gxg W .

    """
    Za = additive_design(real_data)
    n, m = Za.shape

    genes = split_into_genes(Za, G)          # several Z, one per gene

    t_start = time.perf_counter()
    W = build_W_pooled(genes)
    w_build_time = time.perf_counter() - t_start

    Lgxg = cholesky(s2gxg * W + stability * np.eye(n), lower=True)
    return Lgxg, w_build_time


def simulate_phenotype(Lgxg, n, s2gxg=0.5, s2e=0.5, return_realized=False):
    """Draw one phenotype  y = g + e  from the pooled epistasis model,

        g = Lgxg u1 ~ N(0, s2gxg W) ,   W = (1/P) H H' ,
        e = sqrt(s2e) u2 ~ N(0, s2e I) .

    Since Lgxg Lgxg' = s2gxg W = (s2gxg / P) H H', the draw g has EXACTLY the
    law of  H gamma  with  gamma ~ N(0, (s2gxg / P) I_P)  -- the effect-size
    model of realized_variance.pdf (V_gamma = s2gxg) -- without ever forming
    the n-by-P interaction design H.

    REALIZED VARIANCE.  s2gxg (V_gamma in the note) is the per-pair effect
    variance scaled by P; it is NOT the variance of the genetic value across
    individuals.  That quantity,

        V_ell := Var-hat(H gamma) = (1/n) g' P_c g ,   P_c = I - 11'/n ,

    is random (it changes with every gamma draw) and its expectation is
    E[V_ell] = c * s2gxg, c = (1/P) sum_{pairs} Var-hat(Z_a .* Z_b) -- see
    compute_c_pooled.  With return_realized=True the function also returns
    V_ell for THIS draw (ddof=0, matching the note's 1/n normalisation), so
    the simulation can record the target the estimator's c-corrected output
    c-hat * s2gxg-hat is meant to recover.

    Returns y, or (y, V_ell) when return_realized is True.
    """
    u1 = np.random.randn(n)
    u2 = np.random.randn(n)

    gxg = Lgxg @ u1                       # epistasis effect ~ N(0, s2gxg W)
    e = np.sqrt(s2e) * u2                 # residual noise

    y = gxg + e
    if return_realized:
        return y, float(gxg.var())        # realized Var-hat(H gamma), ddof=0
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
    of the genetic value (realized_variance.pdf):

        E[ Var-hat(H gamma) ] = c * s2gxg ,
        c = (1/P) sum_g sum_{a<b in g} Var-hat(Z_a .* Z_b) ,  P = sum_g C(m_g, 2) ,

    for the pooled WITHIN-gene unstandardized kernel: only within-gene pairs
    make up H, exactly the pairs build_W_pooled / setup_pooled count in P.
    The REML fit of  V = s2gxg W + s2e I  with the raw (uncentered, unscaled)
    H estimates V_gamma = s2gxg; the estimate of the realized variance is
    then  V_ell-hat = c * s2gxg-hat.  Because H is not column-standardized,
    c != 1 in general -- it is exactly the per-pair-variance factor the
    standardized siblings divide out inside their kernel.

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
        (Untitled.ipynb) finds ~30x smaller than the draw-to-draw scatter of
        Var-hat(H gamma) itself.

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
def _v_matvec(Z_list, F_list, P, s2gxg, s2e, B):
    """Apply V = s2gxg W + s2e I to B, with W applied matrix-free (rank-r).

        V B = s2gxg (W B) + s2e B .

    B may be (n,) or (n, c); the result matches its shape.  The epistasis term
    is rebuilt from the cached low-rank factors on every call -- no n-by-n W
    exists.
    """
    return s2gxg * compute_WU_pooled(Z_list, F_list, P, B) + s2e * B


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


# ----------------------------------------------------------- MC AI-REML (k=2)
def mc_reml(Z, y, G, iters=30, Nmc=50, cg_tol=1e-6, cg_maxiter=1000,
            jitter=1e-8, tol=1e-8, lm=1e-3, step_frac=0.5, upper_mult=5.0,
            seed=None, verbose=False, r=R_DEFAULT):
    """Monte-Carlo AI-REML for V = s2gxg W + s2e I (pooled UNSTANDARDIZED W,
    applied by its rank-r truncation).

    Z is the column-standardized genotype; it is split into G contiguous genes
    and the pooled setup (per-gene SVD factors, pair total P) is built ONCE
    here, outside the iteration.

    W u is applied by the note's O(n m r) low-rank operator and by nothing
    else: one copy of the Hadamard square (K .* K) is replaced by the rank-r
    eigendecomposition of K taken from the thin SVD of Z_g.  The operator is
    DETERMINISTIC -- W-hat is fixed by (genotype, r), identical across calls
    and replicates with no probe seed to manage -- and its error is a
    truncation BIAS set by the discarded eigenvalue tail: small under LD,
    large under linkage equilibrium, reduced only by raising r (see the module
    header).  verify_lowrank.py holds an independent exact implementation to
    check this one against, including exactness at full rank.

    Score traces
    ------------
    tr(V^{-1} K_i) is estimated by HUTCHINSON and by nothing else: Nmc fixed
    Rademacher probes, W U formed once up front, and Nmc CG solves for V^{-1}U
    per REML iteration (warm-started from the previous iteration).  The probes
    are fixed across iterations, so the objective is a deterministic function
    of s and the iteration converges to a fixed point.

    The stochastic-Lanczos-quadrature alternative that used to sit here is
    gone: it evaluated the SAME Hutchinson target from one cached quadrature
    and needed a conditioning guard plus a fallback to these very solves, so
    the only estimator left is the exact-solve one.  (Its Ritz values also
    doubled as the PSD early-warning on W-hat; without them, an indefinite
    truncation at small r shows up as CG stalling instead -- raise r.)

    A KNOWN DEFECT OF THE SHARED OPTIMIZER, inherited deliberately.  When a
    replicate's likelihood peaks at s2gxg = 0, s2gxg sticks on the lower box
    clamp and s2e then converges to the WRONG value: the AI-Newton step is
    solved jointly and never re-projected onto the free subspace, so the blocked
    gxg direction keeps driving s2e through the coupling term.  It settles where
    (AI^{-1} score)_e = 0 rather than where score_e = 0, i.e. above var(y).
    This is NOT a property of the low-rank operator -- the optimizer block is
    byte-identical across the whole family -- and it does not touch s2gxg,
    which is 0 either way.  It is left uncorrected so that a run here differs
    from a sibling run in the W apply ALONE; an active-set projection would fix
    it but would change exactly the replicates that make the two comparable.

    Returns
    -------
    s  : (2,) estimated (s2gxg, s2e).
    AI : (2, 2) final average-information matrix.
    """
    y = np.asarray(y, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float)
    n = y.shape[0]
    k = 2

    # --- setup: genotype-only, done once ---
    genes = split_into_genes(Z, G)
    F_list, P = setup_pooled(genes, r=r)
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

    WU = Wapply(U)    # W U for the fixed probes: genotype-only, formed once

    for it in range(iters):
        s2gxg, s2e = s
        matvec = lambda B: _v_matvec(genes, F_list, P, s2gxg, s2e, B)

        # --- x = V^{-1} y ---
        xbuf = _cg_batched(matvec, yc, x0=xbuf, tol=cg_tol, maxiter=cg_maxiter)
        x = xbuf[:, 0]

        # data quadratics x'K_i x
        Wx = Wapply(x)
        xWx = x @ Wx
        xIx = x @ x

        # --- score traces tr(V^{-1} K_i): Hutchinson, exact probe solves ---
        Pbuf = _cg_batched(matvec, U, x0=Pbuf, tol=cg_tol, maxiter=cg_maxiter)
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
            print(f"iter {it:2d}  s={s}  max|step|={np.abs(step).max():.3e}")
        if np.abs(step).max() < tol:
            break

    return s, AI


def MC_REML(Z, y, G, iters=30, Nmc=50, cg_tol=1e-6, cg_maxiter=1000, seed=None,
            r=R_DEFAULT, verbose=False):
    """Wrapper: returns (s2gxg_hat, s2e_hat, AI)."""
    s, AI = mc_reml(Z, y, G, iters=iters, Nmc=Nmc, cg_tol=cg_tol,
                    cg_maxiter=cg_maxiter, seed=seed, verbose=verbose, r=r)
    s2gxg_hat, s2e_hat = s
    return s2gxg_hat, s2e_hat, AI


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


