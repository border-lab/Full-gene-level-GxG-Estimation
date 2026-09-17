# -*- coding: utf-8 -*-
"""Check the PRECOMPUTED-W pipeline against its _matfree_ sibling.

    python verify_preW.py [--n 300] [--m 120] [--G 4] [--r 5]
                          [--ref_dir ../Simulation_code_MCREML_pooled_matfree_Lowrank_Wu_std_4VC]

This directory claims to be the _matfree_ sibling with ONE change: REML
applies the epistasis kernel as a precomputed dense matrix instead of the
low-rank linear operator.  Everything below checks that claim at small n, m,
where both can be run side by side.  The sibling module is imported BY PATH
(it has the same module name), so run this where both directories exist.

  1. build_W_lowrank_dense at FULL rank equals build_W_pooled's exact W, and
     the exact W is PSD (lam_min >= -round-off).

  2. build_W_lowrank_dense at rank r equals the sibling's OPERATOR at rank r,
     column by column (W-hat @ I vs compute_WU_pooled(I)) -- so r > 0 here is
     literally the sibling's matrix.

  3. simulate_Cholesky_4vc hands save_W_est the same W it factorized (r = 0),
     and the save/load cache round-trips bit for bit; load_W_cache refuses a
     wrong r and a wrong c-hat.

  4. mc_reml with the dense W reproduces the sibling's mc_reml on the same y
     and seed: at full rank against the exact precomputed W, and at rank r
     against the dense W-hat.  The optimizer is the same code, so the
     estimates must agree to CG tolerance and the V/K/W apply COUNTS should
     match.
"""
import argparse
import importlib.util
import os
import tempfile

import numpy as np

import Function_MCREML as F

HERE = os.path.dirname(os.path.abspath(__file__))


def _rel(a, b):
    return np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300)


def load_ref(ref_dir):
    path = os.path.join(ref_dir, "Function_MCREML.py")
    if not os.path.exists(path):
        raise SystemExit(f"sibling module not found at {path}; pass --ref_dir")
    spec = importlib.util.spec_from_file_location("ref_matfree", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def genotype(n, m, seed=1):
    rng = np.random.default_rng(seed)
    p = rng.uniform(0.1, 0.5, size=m)
    return rng.binomial(2, p, size=(n, m)).astype(float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=300)
    ap.add_argument('--m', type=int, default=120)
    ap.add_argument('--G', type=int, default=4)
    ap.add_argument('--r', type=int, default=5)
    ap.add_argument('--nmc', type=int, default=30)
    ap.add_argument('--ref_dir', default=os.path.join(
        HERE, "..", "Simulation_code_MCREML_pooled_matfree_Lowrank_Wu_std_4VC"))
    a = ap.parse_args()
    ref = load_ref(a.ref_dir)
    n, m, G, r = a.n, a.m, a.G, a.r

    SNP = genotype(n, m)
    Z = F.additive_design(SNP)
    Zd = F.dominance_design(SNP)
    genes = F.split_into_genes(Z, G)
    full = max(min(Zg.shape) for Zg in genes)

    # ---- 1 ---------------------------------------------------------------
    W, c = F.build_W_pooled(genes, return_c=True)
    W_full = F.build_W_lowrank_dense(genes, full)
    e1 = _rel(W_full, W)
    lmin = np.linalg.eigvalsh(W)[0]
    print(f"[1] full-rank dense W-hat vs exact W: rel err {e1:.2e}; "
          f"lam_min(W) = {lmin:.2e}")
    assert e1 < 1e-12, e1
    assert lmin > -1e-10 * np.abs(W).max(), lmin

    # ---- 2 ---------------------------------------------------------------
    W_r = F.build_W_lowrank_dense(genes, r)
    F_list, P, c_ref = ref.setup_pooled(genes, r=r)
    Wop = ref.compute_WU_pooled(genes, F_list, P, c_ref, np.eye(n))
    e2 = _rel(W_r, Wop)
    asym = np.abs(Wop - Wop.T).max() / np.abs(Wop).max()
    print(f"[2] dense rank-{r} W-hat vs sibling operator: rel err {e2:.2e} "
          f"(operator's own asymmetry {asym:.1e}); c-hat equal: {c == c_ref}; "
          f"lam_min(W-hat) = {np.linalg.eigvalsh(W_r)[0]:+.3e}")
    assert e2 < 1e-12, e2
    assert c == c_ref

    # ---- 3 ---------------------------------------------------------------
    got = {}
    La, Ld, Lgxg, _, c_sim, _ = F.simulate_Cholesky_4vc(
        SNP, G, r=0,
        save_W_est=lambda We, cc, psd, t: got.update(W=We.copy(), c=cc, psd=psd))
    assert got['psd'] and got['c'] == c and np.array_equal(got['W'], W)
    tmp = tempfile.mkdtemp(prefix="vpreW_")
    npy, js = F.w_cache_paths(tmp, "Test", n, m, G, 0)
    meta = {'kernel': F.W_CACHE_KERNEL, 'mode': "Test", 'n': n, 'm': m, 'G': G,
            'r': 0, 'exact': True, 'psd': True, 'c_norm': c,
            'c_method': F.C_METHOD, 'build_time': 0.0}
    F.save_W_cache(got['W'], meta, npy, js)
    W_loaded, _ = F.load_W_cache(npy, js, "Test", n, m, G, 0, c)
    assert np.array_equal(W_loaded, W)
    for bad in (dict(r=5), dict(c_expected=c * 1.001)):
        kw = dict(mode="Test", n=n, m=m, G=G, r=0, c_expected=c)
        kw.update(bad)
        try:
            F.load_W_cache(npy, js, **kw)
        except (ValueError, FileNotFoundError):
            continue
        raise AssertionError(f"load_W_cache accepted a mismatched cache {bad}")
    print(f"[3] Cholesky hands over the factorized W; cache round-trips "
          f"bit-exact; mismatched r / c-hat refused  ({npy})")

    # ---- 4 ---------------------------------------------------------------
    np.random.seed(7)
    y = F.simulate_phenotype(La, Ld, Lgxg, n)
    for label, Wd, psd, r_ref in (("exact", W, True, full),
                                  (f"r={r}", W_r, False, r)):
        s_new, _ = F.mc_reml(Z, Zd, Wd, y, Nmc=a.nmc, seed=3, w_psd=psd)
        cnt_new = F.get_op_counts()
        s_ref, _ = ref.mc_reml(Z, Zd, y, G, Nmc=a.nmc, seed=3, r=r_ref)
        cnt_ref = ref.get_op_counts()
        d = np.abs(s_new - s_ref).max()
        same = {k: (cnt_new.get(k), cnt_ref.get(k)) for k in
                ("reml_iters", "cg_iters", "V_columns", "W_columns")}
        print(f"[4] mc_reml {label:>6}: preW {np.round(s_new, 5)}  "
              f"matfree {np.round(s_ref, 5)}  max|diff| {d:.1e}  "
              f"counts (preW, matfree) {same}")
        assert d < 1e-4, d
    print("all checks passed")


if __name__ == "__main__":
    main()
