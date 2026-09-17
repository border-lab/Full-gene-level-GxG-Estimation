# -*- coding: utf-8 -*-
"""
estimate_four_components_plot.py -- the four ESTIMATED variance components of the
nine-column result files as four stacked box panels, in the style of the
single-panel Vl_fixed_* figures.

Usage
-----
    python estimate_four_components_plot.py <result_dir> [--exact] [--paired]
                                            [--drop-boundary] [--ylim LO HI]
                                            [--outdir DIR] [--png]

<result_dir> holds the row-per-replicate files described in calc_stats.py,

    (s2a, s2d, s2gxg, s2e, Vl_hat, vell_gxg, vell_a, vell_d, vell_e)
     1    2    3      4    5       6         7        8        9

Plotted are columns 1, 2, 4 and 5 -- the four things the estimator produces that
are on a variance scale.  Column 3, the raw s2gxg_hat, is deliberately NOT one
of them: it lives on the un-corrected H-component scale, and only
Vl_hat = c_hat * s2gxg_hat is comparable to a realized variance (the same reason
calc_stats.py has no paired difference for gxg).  So the epistasis panel shows
column 5, not column 3.

Each panel is a component's deviation from its reference, one box per sample
size, with a one-sample t-test against zero shown as ns / * / ** / *** above the
box together with its mean and SD -- the layout of
Vl_fixed_m1000_G10_cexact_low_rank.pdf, repeated four times down the page.

References (what zero means)
---------------------------
By default each estimate is centred on the EXPECTED realized variance of its
component, which is what the single-panel figures use:

    col 1  s2a_hat   - s2a
    col 2  s2d_hat   - s2d
    col 4  s2e_hat   - s2e
    col 5  Vl_hat    - c * s2gxg      (c_hat_hwe, or c_exact with --exact)

--paired centres each estimate on THAT REPLICATE's own realized variance
instead -- columns 7, 8, 9 and 6 respectively -- which is calc_stats.py's
paired difference and the sharper target: the draws scatter around the nominal
sigma^2 by O(30%) at m = 1000 and that scatter does not shrink with n, so it
otherwise lands in the estimator's error budget.

--drop-boundary excludes replicates with a fitted component at mc_reml's lower
clamp, exactly as calc_stats.py does.  They are kept by default, for the same
reason: they are a real property of the estimator at these sample sizes.  The
count is reported under the figure either way.

The figure is written to <result_dir> (override with --outdir); --png writes a
PNG next to the PDF for quick viewing.  A ready-to-paste Typst #figure block is
printed at the end.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")

from plot_script import plot_four_effects_single_m

from realized_var_four_components_plot import (
    _parse_truth, _discover, _tag, _read, _expected, COL_VL)

# Fitted columns mc_reml can clamp to its lower bound, and the threshold
# calc_stats.py uses -- a decade above the 1e-9 clamp, so a component that
# merely converged very small is not mistaken for a clamped one.
FITTED_COLS = {"a": 0, "d": 1, "gxg": 2, "e": 3}
CLAMP_LO = 1e-8

# (key, estimate column, realized column for --paired, y-axis label)
PANELS = [
    ("a",   0, 6, r"$\hat{\sigma}^2_a - \sigma^2_a$"),
    ("d",   1, 7, r"$\hat{\sigma}^2_d - \sigma^2_d$"),
    ("e",   3, 8, r"$\hat{\sigma}^2_e - \sigma^2_e$"),
    ("gxg", 4, COL_VL,
     r"$\hat{c}\,\hat{V}_{\gamma} - \mathbb{E}[V_{\ell}^{\mathrm{realized}}]$"),
]

PAIRED_LABELS = {
    "a":   r"$\hat{\sigma}^2_a - \hat{V}_a$",
    "d":   r"$\hat{\sigma}^2_d - \hat{V}_d$",
    "e":   r"$\hat{\sigma}^2_e - \hat{V}_e$",
    "gxg": r"$\hat{c}\,\hat{V}_{\gamma} - \hat{V}_{\ell}$",
}


def _boundary_mask(arr):
    """Boolean mask of replicates with any fitted component at the lower clamp."""
    hit = np.zeros(arr.shape[0], dtype=bool)
    for col in FITTED_COLS.values():
        hit |= arr[:, col] <= CLAMP_LO
    return hit


def _print_caption(pdf_name, m_val, g_val, ns, kept, panel_stats, use_exact,
                   paired, dropped, normalized=False):
    """Emit the Typst #figure block, numbers filled in from this run."""
    gtag = "" if g_val is None else f", $G = {g_val}$"
    ref = ("each replicate's own realized variance (columns 6--9)" if paired
           else "the expected realized variance of the component")
    ctag = "c" if use_exact else "hat(c)"
    if normalized:
        # Column 5 is column 3 here; there is no post-fit c to name.
        gxg_sym = "hat(sigma)^2_(g times g)"
        col3_note = ("Column 5 is identically column 3 in this "
                     "$c$-normalized pipeline, so the epistasis panel is the "
                     "plain $hat(sigma)^2_(g times g)$")
    else:
        gxg_sym = f"{ctag} dot hat(V)_gamma"
        col3_note = ("The raw $hat(sigma)^2_(g times g)$ (column 3) is not "
                     "shown: it lives on the uncorrected $H$-component scale, "
                     f"and only ${gxg_sym}$ is comparable to a variance")
    bits = []
    for key, name in (("a", "additive"), ("d", "dominance"),
                      ("e", "environment"), ("gxg", "epistasis")):
        lo_sd, hi_sd = panel_stats[key]["sd_range"]
        bits.append(f"{name} ${lo_sd:.3f}$--${hi_sd:.3f}$")
    sig = [f"{panel_stats[k]['n_sig']}/{len(ns)}" for k in
           ("a", "d", "e", "gxg")]
    drop_txt = ("" if not dropped else
                f" {dropped} boundary-clamped replicate(s) were excluded.")

    print("\nTypst caption for Simulation.typ:\n")
    print("#figure(")
    print(f'image("{pdf_name}", width: 100%),caption: [The four estimated '
          f'variance components ($m = {m_val}${gtag}, $R = {kept}$ replicates '
          f'per box). Panels, top to bottom: '
          f'$hat(sigma)^2_a$, $hat(sigma)^2_d$, $hat(sigma)^2_e$ and '
          f'${gxg_sym}$ --- columns 1, 2, 4 and 5 of the result '
          f'file. {col3_note}. Each box is '
          f'the deviation from {ref}; the dashed line is zero, and the label '
          f'above each box is a one-sample $t$-test of that deviation against '
          f'zero (`ns` / `*` / `**` / `***`) with its mean and SD. SDs '
          f'contract with $n$ in every panel ({", ".join(bits)}). Significant '
          f'boxes: {sig[0]}, {sig[1]}, {sig[2]}, {sig[3]} of {len(ns)} for '
          f'additive, dominance, environment and epistasis respectively.'
          f'{drop_txt}]')
    print(")")


def main():
    _usage = ("Usage: python estimate_four_components_plot.py <result_dir> "
              "[--exact] [--paired] [--drop-boundary] [--ylim LO HI] "
              "[--outdir DIR] [--png]")
    if len(sys.argv) < 2:
        print(_usage)
        sys.exit(1)

    dir_path = os.path.abspath(sys.argv[1])
    if not os.path.isdir(dir_path):
        print(f"Error: '{dir_path}' is not a directory.")
        sys.exit(1)

    rest = sys.argv[2:]
    use_exact = "--exact" in rest
    paired = "--paired" in rest
    drop_boundary = "--drop-boundary" in rest
    want_png = "--png" in rest
    rest = [t for t in rest if t not in
            ("--exact", "--paired", "--drop-boundary", "--png")]

    shared_ylim = None
    if "--ylim" in rest:
        i = rest.index("--ylim")
        if i + 2 >= len(rest):
            print("Error: --ylim needs two numbers (LO HI).")
            print(_usage)
            sys.exit(1)
        shared_ylim = (float(rest[i + 1]), float(rest[i + 2]))
        if shared_ylim[0] >= shared_ylim[1]:
            print(f"Error: --ylim LO ({shared_ylim[0]}) must be less than "
                  f"HI ({shared_ylim[1]}).")
            sys.exit(1)
        rest = rest[:i] + rest[i + 3:]

    out_dir = dir_path
    if "--outdir" in rest:
        i = rest.index("--outdir")
        if i + 1 >= len(rest):
            print("Error: --outdir needs a directory.")
            print(_usage)
            sys.exit(1)
        out_dir = os.path.abspath(rest[i + 1])
        rest = rest[:i] + rest[i + 2:]
    if rest:
        print(f"Error: unrecognised argument(s) {rest}.")
        print(_usage)
        sys.exit(1)

    folder = os.path.basename(dir_path)
    truth = _parse_truth(folder)
    stem, ns, m_val, g_val = _discover(dir_path)

    print(f"Directory : {dir_path}")
    print(f"Truth     : s2a={truth['s2a']}, s2d={truth['s2d']}, "
          f"s2gxg={truth['s2gxg']}, s2e={truth['s2e']}")
    print(f"Fixed m   : {m_val}" + ("" if g_val is None else f", G = {g_val}"))
    from realized_var_four_components_plot import _read_c
    _probe = _read_c(dir_path, _tag(stem, ns[0], m_val, g_val)) or {}
    if "expected_realized_variance_gxg" in _probe:
        print("Kernel    : c-NORMALIZED -- column 5 is identically column 3, so "
              "the epistasis panel is the plain s2gxg_hat")
        if use_exact:
            print("            (--exact ignored: there is no post-fit c to choose)")
        print("Plotted   : columns 1, 2, 4, 5 -- the estimated components "
              "(all four already on the realized-variance scale)")
    else:
        print(f"c used    : {'c_exact' if use_exact else 'c_hat_hwe'} (epistasis panel)")
        print("Plotted   : columns 1, 2, 4, 5 -- the estimated components "
              "(column 3, raw s2gxg_hat, is off-scale by construction)")
    print(f"Reference : {'per-replicate realized variance (paired)' if paired else 'expected realized variance'}")

    dev = {key: {} for key, _c, _r, _l in PANELS}
    n_dropped, n_boundary, kept_total = 0, 0, None
    for n in ns:
        tag = _tag(stem, n, m_val, g_val)
        arr = _read(os.path.join(dir_path, tag + ".txt"))
        exp_n, src, normalized = _expected(dir_path, tag, truth, arr, use_exact)

        hit = _boundary_mask(arr)
        n_boundary += int(hit.sum())
        if drop_boundary and hit.any():
            arr = arr[~hit]
            n_dropped += int(hit.sum())
            if arr.shape[0] == 0:
                raise SystemExit(f"{tag}: every replicate is boundary-clamped.")
        kept_total = arr.shape[0]

        print(f"  n={n:<6d} R={arr.shape[0]:<4d} [{src}]"
              + (f"  boundary={int(hit.sum())}" if hit.any() else ""))
        for key, est_col, real_col, _lab in PANELS:
            ref = arr[:, real_col] if paired else exp_n[key]
            dev[key][n] = arr[:, est_col] - ref
            v = dev[key][n]
            print(f"      {key:<4s} mean={v.mean():+.6f}  SD={v.std(ddof=1):.6f}")

    if n_boundary:
        print(f"  Boundary: {n_boundary} replicate(s) across the grid have a "
              f"component at the lower clamp"
              + (" -- dropped." if drop_boundary else
                 " -- KEPT (their other components are biased too; "
                 "rerun with --drop-boundary to exclude them)."))

    # A c-normalized kernel makes column 5 identically column 3, so the
    # epistasis panel is the plain s2gxg estimate and must not be labelled as a
    # c-corrected one.
    labels = dict(PAIRED_LABELS) if paired else {k: lab for k, _c, _r, lab in PANELS}
    if normalized:
        labels["gxg"] = (r"$\hat{\sigma}^2_{g \times g} - \hat{V}_{\ell}$"
                         if paired else
                         r"$\hat{\sigma}^2_{g \times g} - \sigma^2_{g \times g}$")
    panels = [(labels[key], dev[key]) for key, _c, _r, _l in PANELS]

    note = None
    if n_boundary:
        note = (f"{n_dropped} boundary-clamped replicate(s) excluded."
                if drop_boundary else
                f"{n_boundary} replicate(s) across the grid sit at the lower "
                f"clamp and are retained.")

    os.makedirs(out_dir, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(out_dir)      # plot_four_effects_single_m writes to the CWD
    try:
        gsuf = "" if g_val is None else f"_G{g_val}"
        base = ("estimate_four_components"
                f"_m{m_val}{gsuf}"
                + ("_cexact" if use_exact else "")
                + ("_paired" if paired else "")
                + ("_nobound" if drop_boundary else ""))
        fig = plot_four_effects_single_m(
            panels,
            individual_sizes=ns,
            m_label=f"m = {m_val:,}" + ("" if g_val is None else f",  G = {g_val}"),
            save_name=base,
            shared_ylim=shared_ylim,
            boundary_note=note,
        )
        if want_png:
            fig.savefig(base + ".png", bbox_inches="tight", facecolor="white",
                        edgecolor="none", format="png", pad_inches=0.1, dpi=200)
            print(f"PNG saved to: {os.path.join(out_dir, base + '.png')}")
    finally:
        os.chdir(cwd)

    from scipy import stats
    panel_stats = {}
    for key, _c, _r, _l in PANELS:
        sds = [dev[key][n].std(ddof=1) for n in ns]
        n_sig = sum(1 for n in ns
                    if stats.ttest_1samp(dev[key][n], 0)[1] < 0.05)
        panel_stats[key] = {"sd_range": (min(sds), max(sds)), "n_sig": n_sig}

    _print_caption(base + ".pdf", m_val, g_val, ns, kept_total, panel_stats,
                   use_exact, paired, n_dropped, normalized)


if __name__ == "__main__":
    main()
