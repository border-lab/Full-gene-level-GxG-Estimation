# -*- coding: utf-8 -*-
"""
gxg_paired_error_plot.py -- the per-replicate epistasis error, c_hat*V_gamma_hat
minus that replicate's OWN realized variance V_ell_hat, as a single box panel in
the SE-figure style.

Usage
-----
    python gxg_paired_error_plot.py <result_dir> [--drop-boundary]
                                    [--ylim LO HI] [--outdir DIR]

<result_dir> holds the row-per-replicate files described in calc_stats.py,

    (s2a, s2d, s2gxg, s2e, Vl_hat, vell_gxg, vell_a, vell_d, vell_e)
     1    2    3      4    5       6         7        8        9

and this figure is columns 5 - 6, one box per sample size.  That difference is
calc_stats.py's paired difference for epistasis, and it is the sharper target
than col 5 - E[V_ell]: the epistasis draws scatter around the nominal
sigma^2_gxg by tens of percent at m = 1000 and that scatter does NOT shrink with
n, so centring on the expectation folds the simulator's own draw noise into the
estimator's error budget.  Pairing removes it row by row.

Column 3, the raw s2gxg_hat, is deliberately not used: it lives on the
uncorrected H-component scale, and only Vl_hat = c_hat * s2gxg_hat is on a
variance scale and hence comparable to a realized variance.  Nothing here reads
c_<tag>.txt for the same reason -- the c correction is already inside column 5,
so the choice of c_hat_hwe vs c_exact (which only sets the EXPECTED reference)
does not enter a paired figure at all.

The signed difference is plotted, not |difference|: the box then shows bias and
spread at once, the dashed line at zero is the unbiasedness target, and the
label above each box is a one-sample t-test of the paired difference against
zero -- the correct paired test, since the reference is this replicate's own
draw rather than something estimated across replicates.

--drop-boundary excludes replicates with a fitted component at mc_reml's lower
clamp, exactly as calc_stats.py does; they are kept by default, since they are a
real property of the estimator at these sample sizes.  The count is reported
either way.

The figure is written as a PDF to <result_dir> (override with --outdir).
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")

from plot_script import plot_relative_error_across_groups_combined
from realized_var_four_components_plot import (
    _parse_truth, _discover, _tag, _read, COL_VL)

COL_VL_HAT = 4          # column 5, c_hat * V_gamma_hat
FITTED_COLS = (0, 1, 2, 3)
CLAMP_LO = 1e-8

YLABEL = r"$\hat{c}\,\hat{V}_{\gamma} - \hat{V}_{\ell}$"


def _boundary_mask(arr):
    """Replicates with any fitted component at mc_reml's lower clamp."""
    hit = np.zeros(arr.shape[0], dtype=bool)
    for col in FITTED_COLS:
        hit |= arr[:, col] <= CLAMP_LO
    return hit


def _auto_ylim(dev, ns, pad=0.15):
    """Limits covering the box-whisker range, symmetric about zero.

    Whiskers follow matplotlib's default 1.5*IQR rule and the helper hides
    outliers (showfliers=False), so the range is taken from the whiskers and
    not from the data extremes -- min/max would leave a band of empty axis
    below the boxes.  Symmetric about zero because zero is the reference line
    this figure is read against; fixed_m.py's _auto_ylimits, whose rule this
    is, does not centre because its panels have no such line.
    """
    half = 0.0
    for n in ns:
        v = dev[n]
        q1, q3 = np.percentile(v, [25, 75])
        iqr = q3 - q1
        w_lo = v[v >= q1 - 1.5 * iqr].min()
        w_hi = v[v <= q3 + 1.5 * iqr].max()
        half = max(half, abs(w_lo), abs(w_hi))
    if half <= 0:
        half = 1.0
    return -half * (1.0 + pad), half * (1.0 + pad)


def _print_caption(pdf_name, m_val, g_val, ns, reps, dev, n_boundary,
                   n_dropped, drop_boundary):
    """Emit the Typst #figure block, numbers filled in from this run."""
    from scipy import stats
    gtag = "" if g_val is None else f", $G = {g_val}$"
    sds = [dev[n].std(ddof=1) for n in ns]
    n_sig = sum(1 for n in ns if stats.ttest_1samp(dev[n], 0)[1] < 0.05)
    drop_txt = ""
    if n_boundary:
        drop_txt = (f" {n_dropped} boundary-clamped replicate(s) were excluded."
                    if drop_boundary else
                    f" {n_boundary} replicate(s) across the grid sit at "
                    f"`mc_reml`'s lower clamp and are retained, as in "
                    f"`calc_stats.py`.")
    print("\nTypst caption for Simulation.typ:\n")
    print("#figure(")
    print(f'image("{pdf_name}", width: 100%),caption: [Per-replicate epistasis '
          f'error ($m = {m_val}${gtag}, $R = {reps}$ replicates per box). Each '
          f'value is $hat(c) dot hat(V)_gamma - hat(V)_ell$ for one replicate '
          f'--- columns 5 and 6 of the result file --- so the reference is the '
          f'realized epistasis variance that *that* replicate\'s own effect '
          f'draw produced, not the nominal $sigma^2_(g times g)$. Pairing '
          f'matters here: the draws scatter around the nominal value by tens '
          f'of percent at $m = {m_val}$ and that scatter does not shrink with '
          f'$n$, so centring on the expectation would fold the simulator\'s '
          f'draw noise into the estimator\'s error. The dashed line is zero '
          f'and the label above each box is a one-sample $t$-test of the '
          f'paired difference against zero (`ns` / `*` / `**` / `***`) with '
          f'its mean and SD; {n_sig} of the {len(ns)} boxes reject. The SD '
          f'goes from ${sds[0]:.3f}$ at $n = {ns[0]}$ to ${sds[-1]:.3f}$ at '
          f'$n = {ns[-1]}$.{drop_txt}]')
    print(")")


def main():
    _usage = ("Usage: python gxg_paired_error_plot.py <result_dir> "
              "[--drop-boundary] [--ylim LO HI] [--outdir DIR]")
    if len(sys.argv) < 2:
        print(_usage)
        sys.exit(1)

    dir_path = os.path.abspath(sys.argv[1])
    if not os.path.isdir(dir_path):
        print(f"Error: '{dir_path}' is not a directory.")
        sys.exit(1)

    rest = sys.argv[2:]
    drop_boundary = "--drop-boundary" in rest
    rest = [t for t in rest if t != "--drop-boundary"]

    ylim = None
    if "--ylim" in rest:
        i = rest.index("--ylim")
        if i + 2 >= len(rest):
            print("Error: --ylim needs two numbers (LO HI).")
            sys.exit(1)
        ylim = (float(rest[i + 1]), float(rest[i + 2]))
        if ylim[0] >= ylim[1]:
            print(f"Error: --ylim LO ({ylim[0]}) must be less than HI ({ylim[1]}).")
            sys.exit(1)
        rest = rest[:i] + rest[i + 3:]

    out_dir = dir_path
    if "--outdir" in rest:
        i = rest.index("--outdir")
        if i + 1 >= len(rest):
            print("Error: --outdir needs a directory.")
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
    print("Plotted   : column 5 - column 6, per replicate "
          "(c_hat*V_gamma_hat - V_ell_hat)")

    dev, n_boundary, n_dropped, reps = {}, 0, 0, None
    for n in ns:
        tag = _tag(stem, n, m_val, g_val)
        arr = _read(os.path.join(dir_path, tag + ".txt"))

        hit = _boundary_mask(arr)
        n_boundary += int(hit.sum())
        if drop_boundary and hit.any():
            arr = arr[~hit]
            n_dropped += int(hit.sum())
            if arr.shape[0] == 0:
                raise SystemExit(f"{tag}: every replicate is boundary-clamped.")
        reps = arr.shape[0]

        d = arr[:, COL_VL_HAT] - arr[:, COL_VL]
        dev[n] = d
        se = d.std(ddof=1) / np.sqrt(d.size)
        print(f"  n={n:<6d} R={d.size:<4d} mean={d.mean():+.6f} "
              f"SD={d.std(ddof=1):.6f} SE={se:.6f} "
              f"median={np.median(d):+.6f} "
              f"mean|err|={np.abs(d).mean():.6f}"
              + (f"  boundary={int(hit.sum())}" if hit.any() else ""))

    if n_boundary:
        print(f"  Boundary: {n_boundary} replicate(s) across the grid have a "
              f"component at the lower clamp"
              + (" -- dropped." if drop_boundary else
                 " -- KEPT (rerun with --drop-boundary to exclude them)."))

    if ylim is None:
        ylim = _auto_ylim(dev, ns)

    # plot_relative_error_across_groups_combined subtracts a single scalar
    # reference; the reference here is per-replicate, so it is already inside
    # dev and the scalar is 0.
    data_dict = {n: pd.DataFrame({"paired": dev[n]}) for n in ns}

    gsuf = "" if g_val is None else f"_G{g_val}"
    base = (f"gxg_paired_error_m{m_val}{gsuf}"
            + ("_nobound" if drop_boundary else ""))

    os.makedirs(out_dir, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(out_dir)      # the plot helper writes to the CWD
    try:
        plot_relative_error_across_groups_combined(
            data_dict,
            x_labels=[f"m = {m_val:,}"
                      + ("" if g_val is None else f",  G = {g_val}")],
            individual_sizes=ns,
            col_num=0,
            real_value=0.0,
            ymin=ylim[0],
            ymax=ylim[1],
            x_axis_name="Sample size (n)",
            y_axis_name=YLABEL,
            save_name=base,
        )
    finally:
        os.chdir(cwd)

    print(f"Saved     : {os.path.join(out_dir, base + '.pdf')}  "
          f"(y in [{ylim[0]:.3f}, {ylim[1]:.3f}])")

    _print_caption(base + ".pdf", m_val, g_val, ns, reps, dev,
                   n_boundary, n_dropped, drop_boundary)


if __name__ == "__main__":
    main()
