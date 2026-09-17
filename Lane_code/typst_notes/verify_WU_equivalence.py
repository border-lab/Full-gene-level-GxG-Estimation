# -*- coding: utf-8 -*-
"""
Show that the IMPLICIT (matrix-free) W-apply equals the EXPLICIT dense-W apply.

Both compute the action  W u  of the pairwise-epistasis GRM

    W = (1/p) sum_{a<b} h_ab h_ab',                       p = m(m-1)/2
    h_ab = standardize(Z_a . Z_b)                         (mean 0, var 1)

where Z is the column-standardized genotype and "." is the element-wise
(Hadamard) product of two SNP columns.

  * EXPLICIT: build the n-by-n W by summing h_ab h_ab' over all pairs, then W @ U.
              O(n^2 p) to build, O(n^2 c) to apply, O(n^2) storage.

  * IMPLICIT: apply W u WITHOUT forming W (or the n-by-p matrix of h_ab columns),
              using only m-by-m weight matrices contracted from Z.
              O(n m^2 c) to apply, O(m^2) storage — the large-m / small-n regime.

The two functions below are written from the model, not copied from the
pipelines; the script checks compute_WU_implicit == W_explicit @ U to
floating-point precision.

Run:  python verify_WU_equivalence.py

----------------------------------------------------------------------------
Derivation of the implicit apply
----------------------------------------------------------------------------
Write the standardized interaction column for pair (a,b) as

    h_ab = (d_ab - mu_ab) / sigma_ab ,     d_ab = Z_a . Z_b  (Hadamard),

with, over the n individuals (ddof=0),

    mu_ab   = E[Z_a Z_b]                    = (Z^T Z)_{ab} / n              =: E_ab
    sigma_ab^2 = Var(Z_a Z_b)
              = E[(Z_a Z_b)^2] - E[Z_a Z_b]^2
              = ((Z.Z)^T (Z.Z))_{ab} / n  -  E_ab^2                        =: Vp_ab

Then  W u = (1/p) sum_{a<b} h_ab (h_ab^T u).  Substituting h_ab and expanding
the product (d_ab - mu_ab)(d_ab - mu_ab)^T u gives four terms.  Define the
m-by-m weight matrices (diagonals zeroed, since a < b are the off-diagonal
pairs and each unordered pair is counted once via the 1/2 factor)

    S = 1 / Vp ,     R = E / Vp ,     T = E.^2 / Vp .

For an individual i and column u, the four terms are, per pair (a,b):

    d_ab (d_ab^T u) S_ab   -> Z_ia Z_ib * S_ab * sum_j Z_ja Z_jb u_j
    - mu_ab d_ab^T u  R'   ... etc.

Collapsing the sums over pairs into matrix contractions M = Z^T diag(u) Z
yields the storage-free form used below (identical algebra to the pipelines'
compute_WU, re-derived here).
"""

import numpy as np


# --------------------------------------------------------------------------
# Shared design: column-standardize the genotype (mean 0, variance 1, ddof=0).
# --------------------------------------------------------------------------
def standardize_cols(M, eps=1e-12):
    M = np.asarray(M, dtype=float)
    mu = M.mean(axis=0)
    sd = M.std(axis=0)
    sd = np.where(sd < eps, 1.0, sd)
    return (M - mu) / sd


# --------------------------------------------------------------------------
# EXPLICIT operator: build W densely from every standardized pair, apply W @ U.
# This is the reference — a literal transcription of the model definition.
# --------------------------------------------------------------------------
def build_W_explicit(Z, pair_batch=4096):
    """W = (1/p) sum_{a<b} h_ab h_ab', formed as a dense n-by-n matrix."""
    n, m = Z.shape
    p = m * (m - 1) // 2
    ia, ib = np.triu_indices(m, k=1)          # all unordered pairs a < b

    W = np.zeros((n, n))
    for s in range(0, p, pair_batch):
        e = min(s + pair_batch, p)
        H = Z[:, ia[s:e]] * Z[:, ib[s:e]]     # raw Hadamard products, (n, batch)
        mu = H.mean(axis=0)
        sd = H.std(axis=0, ddof=0)
        ok = sd > 1e-10                       # drop degenerate (constant) pairs
        H[:, ok] = (H[:, ok] - mu[ok]) / sd[ok]
        H[:, ~ok] = 0.0
        W += H @ H.T
    return W / p


def apply_W_explicit(Z, U):
    """Explicit action: form dense W, return W @ U."""
    return build_W_explicit(Z) @ U


# --------------------------------------------------------------------------
# IMPLICIT operator: apply W u with no n-by-n (or n-by-p) object, via the
# m-by-m weight matrices S, R, T contracted from Z.  See the header derivation.
# --------------------------------------------------------------------------
def weight_matrices(Z):
    """m-by-m (S, R, T) encoding the per-pair standardization of h_ab.

        E  = E[Z_a Z_b]        = (Z^T Z) / n
        Vp = Var(Z_a Z_b)      = ((Z.Z)^T (Z.Z)) / n - E^2
        S = 1/Vp ,  R = E/Vp ,  T = E^2/Vp ,   diagonals zeroed (a != b).
    """
    n = Z.shape[0]
    E = (Z.T @ Z) / n
    D = Z * Z                                  # element-wise square, (n, m)
    Vp = (D.T @ D) / n - E * E
    S = 1.0 / Vp
    R = E / Vp
    T = (E * E) / Vp
    np.fill_diagonal(S, 0.0)
    np.fill_diagonal(R, 0.0)
    np.fill_diagonal(T, 0.0)
    return S, R, T


def apply_W_implicit(Z, U, S, R, T):
    """Matrix-free action  W @ U  built from (S, R, T); never forms W or H.

    For each column u, the m-by-m contraction  M = Z^T diag(u) Z = Z^T (u . Z)
    carries the pair sums; the four expansion terms of
    sum_{a<b} h_ab (h_ab^T u) collapse to per-individual reductions over Z.
    Cost O(n m^2) per column — the M formation dominates.
    """
    n, m = Z.shape
    c = U.shape[1]
    p = m * (m - 1) // 2

    # Terms that do NOT depend on u through M (only through the column sum of u).
    ZR = Z @ R                                 # (n, m)
    quad_R = 0.5 * np.sum(Z * ZR, axis=1, keepdims=True)   # 0.5 * z_i^T R z_i, (n,1)
    sum_T = 0.5 * np.sum(T)                                 # scalar
    colsum_U = np.sum(U, axis=0, keepdims=True)            # (1, c) = sum_i u_i

    term_mu = quad_R @ colsum_U                # (n, c)
    term_TT = (sum_T * np.ones((n, 1))) @ colsum_U         # (n, c)

    term_dd = np.zeros((n, c))                 # d_ab d_ab^T u part (weight S)
    term_dm = np.zeros((n, c))                 # cross term        (weight R)
    for k in range(c):
        u = U[:, k]
        M = Z.T @ (u[:, None] * Z)             # (m, m) = Z^T diag(u) Z
        term_dd[:, k] = 0.5 * np.sum(Z * (Z @ (S * M)), axis=1)
        term_dm[:, k] = 0.5 * np.sum(R * M)    # scalar broadcast over individuals

    return (term_dd - term_mu - term_dm + term_TT) / p


# --------------------------------------------------------------------------
# Genotype generators
# --------------------------------------------------------------------------
def gaussian_geno(n, m, rng):
    """Continuous N(0,1) 'dosages' — a clean, non-degenerate stress test."""
    return rng.standard_normal((n, m))


def snp_geno(n, m, rng, maf_lo=0.05, maf_hi=0.5):
    """0/1/2 minor-allele counts under HWE with random MAFs (real-SNP-like)."""
    maf = rng.uniform(maf_lo, maf_hi, size=m)
    G = rng.binomial(2, maf[None, :], size=(n, m)).astype(float)
    for j in range(m):                         # redraw any monomorphic column
        while G[:, j].std() < 1e-8:
            G[:, j] = rng.binomial(2, maf[j], size=n)
    return G


# --------------------------------------------------------------------------
# Comparison
# --------------------------------------------------------------------------
def compare_once(n, m, c, geno_fn, rng, tol=1e-8):
    Z = standardize_cols(geno_fn(n, m, rng))
    U = rng.choice([-1.0, 1.0], size=(n, c))   # Rademacher block

    W = build_W_explicit(Z)
    WU_exp = W @ U

    S, R, T = weight_matrices(Z)
    WU_imp = apply_W_implicit(Z, U, S, R, T)

    abs_err = np.max(np.abs(WU_imp - WU_exp))
    denom = np.max(np.abs(WU_exp))
    rel_err = abs_err / denom if denom > 0 else abs_err

    # feeding the identity through the implicit apply must rebuild W itself
    W_imp = apply_W_implicit(Z, np.eye(n), S, R, T)
    op_err = np.max(np.abs(W_imp - W))

    return dict(n=n, m=m, c=c, geno=geno_fn.__name__,
                abs_err=abs_err, rel_err=rel_err, op_err=op_err,
                passed=rel_err < tol)


# --------------------------------------------------------------------------
# Identities behind the optimized apply (Wu_complexity.typ, Strategies 1 & 2)
#
# NAMING: this file's `S` (= 1/Vp) is the note's `V`.  The note's
# `S = V .* M` is the expression `S * M` here.  Below, `Sn` denotes the note's S.
# --------------------------------------------------------------------------
def check_dg_identity(Z, Sn):
    """dg(Z Sn Z^T) computed three ways — all must agree.

    Level 1  np.diag(Z Sn Z^T)                  O(n^2 m) time, O(n^2) space  [never do this]
    Level 2  (Z .* (Z Sn)) 1                    O(n m^2) time, O(n m) space  [current code]
    Level 3  2 sum_{a<b} Sn_ab Z_ia Z_ib        O(n m^2 / 2), no n-by-m buffer [Strategy 2b]

    Level 3 needs Sn symmetric with zero diagonal; levels 1 and 2 hold for any
    square Sn.
    """
    n, m = Z.shape
    lvl1 = np.diag(Z @ Sn @ Z.T)
    lvl2 = np.sum(Z * (Z @ Sn), axis=1)
    ia, ib = np.triu_indices(m, k=1)
    lvl3 = 2.0 * (Z[:, ia] * Z[:, ib]) @ Sn[ia, ib]
    return {
        "dg: level1(diag of ZSZ') vs level2(row-sums)": _err(lvl1, lvl2),
        "dg: level2(row-sums)     vs level3(a<b pairs)": _err(lvl2, lvl3),
    }


def check_strategy1(Z, R, T, u):
    """t_3 needs no M:  1'(R .* M)1 = <R,M> = u' dg(Z R Z^T) = u' v_R."""
    M = Z.T @ (u[:, None] * Z)
    v_R = np.sum(Z * (Z @ R), axis=1)            # dg(Z R Z^T), precomputable
    return {
        "Str1: v_R = dg(Z R Z^T)":            _err(v_R, np.diag(Z @ R @ Z.T)),
        "Str1: 1'(R.*M)1 = <R,M>":            _err(np.sum(R * M), np.trace(R @ M)),
        "Str1: <R,M>     = u' v_R":           _err(np.sum(R * M), u @ v_R),
    }


def check_strategy2a(Z, u):
    """M = Z' diag(u) Z = A'A - B'B, the sign-split syrk form."""
    M_ref = Z.T @ (u[:, None] * Z)
    pos, neg = u > 0, u < 0
    A = Z[pos] * np.sqrt(u[pos])[:, None]
    B = Z[neg] * np.sqrt(-u[neg])[:, None]
    out = {"Str2a: M = A'A - B'B (sign split)": _err(M_ref, A.T @ A - B.T @ B)}

    try:                                          # the actual BLAS route
        from scipy.linalg.blas import dsyrk
        tri = dsyrk(1.0, A, trans=1) - dsyrk(1.0, B, trans=1)   # upper triangle
        M_syrk = tri + tri.T - np.diag(np.diag(tri))
        out["Str2a: M via two dsyrk calls"] = _err(M_ref, M_syrk)
    except ImportError:                           # pragma: no cover
        pass
    return out


def apply_W_final(Z, U, V, v_R, s_T):
    """The boxed formula of Wu_complexity.typ, in block form.

        W U = 1/(2p) [ T_1 - v_R s_U' + 1 (s_T s_U - U' v_R)' ],
        T_1[:,k] = dg(Z (V .* M_k) Z^T),   M_k = Z' diag(u_k) Z.

    Takes the precomputed v_R = dg(Z R Z^T) and s_T = 1'T1; R and T themselves
    are no longer needed.
    """
    n, m = Z.shape
    p = m * (m - 1) // 2
    s_U = U.sum(axis=0)                                     # (c,)

    T_1 = np.empty_like(U, dtype=float)
    for k in range(U.shape[1]):
        M = Z.T @ (U[:, k, None] * Z)
        T_1[:, k] = np.sum(Z * (Z @ (V * M)), axis=1)       # dg(Z (V.*M) Z^T)

    rank1 = np.outer(v_R, s_U)                              # v_R s_U'
    scal = s_T * s_U - U.T @ v_R                            # (c,)
    return (T_1 - rank1 + np.outer(np.ones(n), scal)) / (2 * p)


def _err(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    denom = max(np.max(np.abs(a)), 1e-300)
    return np.max(np.abs(a - b)) / denom


def verify_note_identities(n=180, m=40, c=6, seed=7, tol=1e-10):
    """Check every algebraic step the note relies on, plus the final formula."""
    rng = np.random.default_rng(seed)
    Z = standardize_cols(snp_geno(n, m, rng))
    U = rng.choice([-1.0, 1.0], size=(n, c))
    u = U[:, 0]

    V, R, T = weight_matrices(Z)          # this file's S is the note's V
    M = Z.T @ (u[:, None] * Z)
    Sn = V * M                            # the note's S = V .* M

    checks = {}
    checks.update(check_dg_identity(Z, Sn))
    checks.update(check_strategy1(Z, R, T, u))
    checks.update(check_strategy2a(Z, u))

    # the note's boxed formula vs. the four-term implicit apply vs. dense W
    v_R = np.sum(Z * (Z @ R), axis=1)
    s_T = np.sum(T)
    WU_final = apply_W_final(Z, U, V, v_R, s_T)
    checks["Final: boxed formula = 4-term implicit"] = _err(
        WU_final, apply_W_implicit(Z, U, V, R, T))
    checks["Final: boxed formula = dense W @ U"] = _err(
        WU_final, build_W_explicit(Z) @ U)

    print(f"\nIdentities from Wu_complexity.typ   (n={n}, m={m}, c={c})")
    print("-" * 84)
    ok = True
    for name, err in checks.items():
        ok &= err < tol
        print(f"  {name:<48} rel_err = {err:>10.3e}  "
              f"{'OK' if err < tol else 'FAIL'}")
    print("-" * 84)
    if not ok:
        raise SystemExit("SOME IDENTITIES FAILED")
    print(f"ALL IDENTITIES HOLD  (rel_err < {tol:g})")
    return ok


def main():
    rng = np.random.default_rng(0)
    tol = 1e-8
    cases = [
        (200, 30, 1),
        (200, 30, 8),
        (300, 50, 5),
        (500, 80, 16),
        (400, 120, 4),
        (150, 200, 10),     # m > n: the matrix-free target regime
    ]
    geno_fns = [gaussian_geno, snp_geno]

    print(f"{'geno':>14} {'n':>5} {'m':>5} {'c':>4} "
          f"{'abs_err':>12} {'rel_err':>12} {'op_err(W)':>12}  status")
    print("-" * 84)

    all_pass = True
    for geno_fn in geno_fns:
        for (n, m, c) in cases:
            r = compare_once(n, m, c, geno_fn, rng, tol=tol)
            all_pass &= r["passed"]
            print(f"{r['geno']:>14} {n:>5} {m:>5} {c:>4} "
                  f"{r['abs_err']:>12.3e} {r['rel_err']:>12.3e} "
                  f"{r['op_err']:>12.3e}  {'OK' if r['passed'] else 'FAIL'}")

    print("-" * 84)
    if all_pass:
        print(f"ALL CASES PASS  (max relative error < tol = {tol:g})")
        print("=> implicit  apply_W_implicit(Z, U, S, R, T)  "
              "==  explicit  build_W_explicit(Z) @ U")
    else:
        print("SOME CASES FAILED — the two W-apply operators disagree.")
        raise SystemExit(1)

    verify_note_identities()


if __name__ == "__main__":
    main()
