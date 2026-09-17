# -*- coding: utf-8 -*-
"""Check the EXACT AI-REML estimator at small n, m.

    python verify_exact.py [--n 400] [--m 60] [--G 3] [--nmc 200]
                           [--ref_dir ../Simulation_code_MCREML_pooled_preW_Lowrank_Wu_std_4VC]

  1. build_W_pooled's W is PSD; simulate_Cholesky_4vc hands save_W the very W
     it factorized; the save/load cache round-trips bit for bit, and
     load_W_cache refuses a wrong G and a wrong c-hat.

  2. _factor_V returns the true inverse and log-determinant of V (against
     numpy's dense inv / slogdet), and an exactly symmetric V^{-1}.

  3. At an arbitrary feasible point, _score_ai's score equals an INDEPENDENT
     dense-numpy gradient of logL = -0.5 (log|V| + y'V^{-1}y) and a central
     finite difference of it, and its AI equals 0.5 (K_i x)'V^{-1}(K_j x)
     formed with np.linalg.solve.

  4. ai_reml reaches the likelihood maximum: its estimate and logL match
     scipy L-BFGS-B run on the same independent logL.

  5. (if the MC sibling directory exists) the sibling's Monte-Carlo mc_reml on
     the same y and W lands within Monte-Carlo error of the exact fit.
"""
import argparse
import importlib.util
import os
import tempfile

import numpy as np
from scipy.optimize import minimize

import Function_AIREML as F

HERE = os.path.dirname(os.path.abspath(__file__))


def _rel(a, b):
    return np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300)


def genotype(n, m, seed=1):
    rng = np.random.default_rng(seed)
    p = rng.uniform(0.1, 0.5, size=m)
    return rng.binomial(2, p, size=(n, m)).astype(float)


def ref_loglik(Ks, y, s, grad=True):
    """logL and its gradient by plain dense numpy -- no Cholesky, no potri."""
    n = y.shape[0]
    K4 = list(Ks) + [np.eye(n)]
    V = sum(si * Ki for si, Ki in zip(s, K4))
    sign, logdet = np.linalg.slogdet(V)
    if sign <= 0:
        return -np.inf, None
    x = np.linalg.solve(V, y)
    ll = -0.5 * (logdet + y @ x)
    if not grad:
        return ll, None
    g = np.array([0.5 * (x @ K @ x - np.trace(np.linalg.solve(V, K)))
                  for K in K4])
    return ll, g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=400)
    ap.add_argument('--m', type=int, default=60)
    ap.add_argument('--G', type=int, default=3)
    ap.add_argument('--nmc', type=int, default=200)
    ap.add_argument('--ref_dir', default=os.path.join(
        HERE, "..", "Simulation_code_MCREML_pooled_preW_Lowrank_Wu_std_4VC"))
    a = ap.parse_args()
    n, m, G = a.n, a.m, a.G
    s0 = np.array([0.2, 0.1, 0.3, 0.4])

    SNP = genotype(n, m)
    Z = F.additive_design(SNP)
    Zd = F.dominance_design(SNP)
    genes = F.split_into_genes(Z, G)

    # ---- 1 ---------------------------------------------------------------
    W, c = F.build_W_pooled(genes, return_c=True)
    lmin = np.linalg.eigvalsh(W)[0]
    assert lmin > -1e-10 * np.abs(W).max(), lmin
    got = {}
    La, Ld, Lgxg, _, c_sim = F.simulate_Cholesky_4vc(
        SNP, G, *s0, save_W=lambda We, cc, t: got.update(W=We.copy(), c=cc))
    assert c_sim == c and got['c'] == c and np.array_equal(got['W'], W)
    tmp = tempfile.mkdtemp(prefix="vexact_")
    npy, js = F.w_cache_paths(tmp, "Test", n, m, G)
    meta = {'kernel': F.W_CACHE_KERNEL, 'mode': "Test", 'n': n, 'm': m, 'G': G,
            'exact': True, 'psd': True, 'c_norm': c, 'c_method': F.C_METHOD,
            'build_time': 0.0}
    F.save_W_cache(got['W'], meta, npy, js)
    W_loaded, _ = F.load_W_cache(npy, js, "Test", n, m, G, c)
    assert np.array_equal(W_loaded, W)
    for bad in (dict(G=G + 1), dict(c_expected=c * 1.001)):
        kw = dict(mode="Test", n=n, m=m, G=G, c_expected=c)
        kw.update(bad)
        try:
            F.load_W_cache(npy, js, **kw)
        except (ValueError, FileNotFoundError):
            continue
        raise AssertionError(f"load_W_cache accepted a mismatched cache {bad}")
    print(f"[1] W PSD (lam_min {lmin:.1e}); Cholesky hands over the factorized "
          f"W; cache round-trips bit-exact; wrong G / c-hat refused")

    # ---- 2 ---------------------------------------------------------------
    Ka, Kd = F.build_K(Z), F.build_K(Zd)
    Ks = (Ka, Kd, W)
    V = s0[0] * Ka + s0[1] * Kd + s0[2] * W + s0[3] * np.eye(n)
    _, Vinv, logdet = F._factor_V(Ka, Kd, W, s0)
    e_inv = _rel(Vinv, np.linalg.inv(V))
    e_ld = abs(logdet - np.linalg.slogdet(V)[1])
    print(f"[2] _factor_V: rel err inverse {e_inv:.1e}, |log|V| err| {e_ld:.1e}, "
          f"symmetric {np.array_equal(Vinv, Vinv.T)}")
    assert e_inv < 1e-10 and e_ld < 1e-8 and np.array_equal(Vinv, Vinv.T)

    # ---- 3 ---------------------------------------------------------------
    np.random.seed(7)
    y = F.simulate_phenotype(La, Ld, Lgxg, n, *s0)
    s_pt = s0 * np.array([1.3, 0.7, 1.2, 0.9])
    sc, ai, ll = F._score_ai(Ka, Kd, W, y, s_pt)
    ll_ref, g_ref = ref_loglik(Ks, y, s_pt)
    h = 1e-5
    fd = np.array([(ref_loglik(Ks, y, s_pt + h * e, grad=False)[0]
                    - ref_loglik(Ks, y, s_pt - h * e, grad=False)[0]) / (2 * h)
                   for e in np.eye(4)])
    V = sum(si * Ki for si, Ki in zip(s_pt, list(Ks) + [np.eye(n)]))
    x = np.linalg.solve(V, y)
    KX = np.column_stack([Ka @ x, Kd @ x, W @ x, x])
    ai_ref = 0.5 * KX.T @ np.linalg.solve(V, KX)
    e_g, e_fd, e_ai = _rel(sc, g_ref), _rel(sc, fd), _rel(ai, ai_ref)
    e_ll = abs(ll - ll_ref) / abs(ll_ref)
    print(f"[3] score vs dense gradient {e_g:.1e}, vs finite difference "
          f"{e_fd:.1e}; AI {e_ai:.1e}; logL {e_ll:.1e}")
    assert e_g < 1e-9 and e_fd < 1e-5 and e_ai < 1e-9 and e_ll < 1e-10

    # ---- 4 ---------------------------------------------------------------
    s_hat, AI = F.ai_reml(Z, Zd, W, y)
    cnt = F.get_op_counts()
    se = F.reml_se(AI)
    vary = y.var()
    su = 1.5 * vary
    bounds = [(-0.2 * su, su)] * 3 + [(1e-9, su)]

    def negll(s):
        ll_, g_ = ref_loglik(Ks, y, s)
        if not np.isfinite(ll_):
            return 1e10, np.zeros(4)
        return -ll_, -g_

    res = minimize(negll, vary * np.array([0.05, 0.05, 0.05, 0.85]), jac=True,
                   method='L-BFGS-B', bounds=bounds,
                   options=dict(ftol=1e-15, gtol=1e-9, maxiter=5000))
    ll_hat = ref_loglik(Ks, y, s_hat, grad=False)[0]
    dz = np.abs(s_hat - res.x) / se
    print(f"[4] ai_reml  {np.round(s_hat, 5)}  logL {ll_hat:.6f}  SE {np.round(se, 4)}\n"
          f"    L-BFGS-B {np.round(res.x, 5)}  logL {-res.fun:.6f}\n"
          f"    max|diff|/SE {dz.max():.2e}; counts {cnt}")
    assert ll_hat >= -res.fun - 1e-3, (ll_hat, -res.fun)
    assert dz.max() < 0.05, dz

    # ---- 5 ---------------------------------------------------------------
    path = os.path.join(a.ref_dir, "Function_MCREML.py")
    if os.path.exists(path):
        spec = importlib.util.spec_from_file_location("ref_mc", path)
        ref = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(ref)
        s_mc, _ = ref.mc_reml(Z, Zd, W, y, Nmc=a.nmc, seed=3, w_psd=True)
        dmc = np.abs(s_mc - s_hat) / se
        print(f"[5] MC sibling (Nmc={a.nmc}) {np.round(s_mc, 5)}  "
              f"|MC - exact|/SE {np.round(dmc, 3)}")
        assert dmc.max() < 1.0, dmc
    else:
        print(f"[5] skipped: no MC sibling at {path}")
    print("all checks passed")


if __name__ == "__main__":
    main()
