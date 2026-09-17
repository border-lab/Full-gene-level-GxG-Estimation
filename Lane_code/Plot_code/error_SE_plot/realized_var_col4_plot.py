# -*- coding: utf-8 -*-
"""
realized_var_col4_plot.py -- fixed-m, increasing-n figure for the FOURTH column
alone: the simulated realized variance V_ell^realized = Var(H gamma) of each
replicate's own gamma draw.

Usage
-----
    python realized_var_col4_plot.py <result_dir> [ymin ymax] [--outdir DIR]

This is the diagnostic companion to realized_var_perrep_plot.py.  That figure
plots col3 - col4 and is dominated by the spread of col4, so this one shows col4
on its own to make its n-behaviour visible.

col4 is not an estimate -- it is the truth of its replicate -- so it is centred
on its own expectation, E[V_ell^realized] = c_n * s2gxg, where c_n is the
genotype-only constant the estimator applied at that n (recovered exactly as
col3 / col1, which is constant within a file).  Centring makes the reference
line 0 and the one-sample t-test meaningful: it asks whether the simulated
realized variance sits on its theoretical mean.  The raw per-n means and c_n are
printed to the console so the absolute scale is not hidden.

What to read off the figure: the box HEIGHTS.  Var(V_ell^realized) is set by the
number of effective interaction dimensions of H, not by n, so the boxes should
be flat across n -- that flat height is the irreducible floor under
SD(col3 - col4).

Default output directory is <result_dir> itself (override with --outdir).
"""
import os
import re
import sys
import glob

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")

from plot_script import plot_relative_error_across_groups_combined

from realized_var_perrep_plot import _parse_truth, _discover, _tag, _read


def main():
    _usage = ("Usage: python realized_var_col4_plot.py <result_dir> "
              "[ymin ymax] [--outdir DIR]")
    if len(sys.argv) < 2:
        print(_usage)
        sys.exit(1)

    dir_path = os.path.abspath(sys.argv[1])
    if not os.path.isdir(dir_path):
        print(f"Error: '{dir_path}' is not a directory.")
        sys.exit(1)

    out_dir = dir_path          # store next to the data by default
    rest = sys.argv[2:]

    if "--outdir" in rest:
        i = rest.index("--outdir")
        if i + 1 >= len(rest):
            print("Error: --outdir needs a directory.")
            print(_usage)
            sys.exit(1)
        out_dir = os.path.abspath(rest[i + 1])
        rest = rest[:i] + rest[i + 2:]

    nums = [float(t) for t in rest]
    if len(nums) == 2:
        ymin, ymax = nums
        if ymin >= ymax:
            print(f"Error: ymin ({ymin}) must be less than ymax ({ymax}).")
            sys.exit(1)
    elif len(nums) == 0:
        ymin, ymax = None, None
    else:
        print("Error: need exactly two numbers (ymin ymax).")
        print(_usage)
        sys.exit(1)

    folder = os.path.basename(dir_path)
    real_gxg, real_e = _parse_truth(folder)
    stem, ns, m_val, g_val = _discover(dir_path)
    data = _read(dir_path, stem, ns, m_val, g_val)

    print(f"Directory : {dir_path}")
    print(f"Truth     : s2gxg={real_gxg}, s2e={real_e}")
    print(f"Fixed m   : {m_val}" + ("" if g_val is None else f", G = {g_val}"))
    print("Plotted   : col4 (V_realized) - c_n*s2gxg, its own expectation")

    centred = {}
    for n in ns:
        df = data[n]

        # c_n the estimator applied; col3 = c_n * col1 exactly within a file.
        ratio = (df[2] / df[0]).values
        c_n = float(np.median(ratio))
        spread = float(np.max(np.abs(ratio - c_n)))
        if spread > 1e-8 * max(1.0, abs(c_n)):
            raise ValueError(f"n={n}: col3/col1 is not constant (spread {spread:.2e}); "
                             "c_n cannot be recovered from the results.")

        truth_n = c_n * real_gxg
        v = df[3].values
        centred[n] = pd.DataFrame({0: v - truth_n})
        print(f"  n={n:<6d} R={len(v):<4d} c_n={c_n:.4f}  E[V_real]=c_n*s2gxg={truth_n:.6f}")
        print(f"           mean(V_real)={v.mean():.6f}  SD={v.std(ddof=1):.6f}  "
              f"deviation from E = {v.mean() - truth_n:+.6f}")

    if ymin is None:
        lo, hi = np.inf, -np.inf
        for n in ns:
            x = centred[n][0].values
            q1, q3 = np.percentile(x, [25, 75])
            iqr = q3 - q1
            lo = min(lo, x[x >= q1 - 1.5 * iqr].min())
            hi = max(hi, x[x <= q3 + 1.5 * iqr].max())
        span = (hi - lo) or 1.0
        ymin, ymax = lo - 0.15 * span, hi + 0.15 * span
    print(f"y-axis    : [{ymin}, {ymax}]")

    # plot_relative_error_across_groups_combined writes its PDF to the CWD.
    os.makedirs(out_dir, exist_ok=True)
    os.chdir(out_dir)

    gtag = "" if g_val is None else f"_G{g_val}"
    # Writing beside the data: the folder already names the setting, and
    # repeating it in the filename overruns the Windows 260-char path limit.
    prefix = "" if os.path.normcase(out_dir) == os.path.normcase(dir_path) else folder + "_"
    plot_relative_error_across_groups_combined(
        centred,
        x_labels=[f"m = {m_val:,}"],
        individual_sizes=ns,
        col_num=0,
        real_value=0.0,          # already centred on the per-n expectation
        ymin=ymin,
        ymax=ymax,
        x_axis_name="Sample size (n)",
        y_axis_name=r"$V_{\ell}^{\mathrm{realized}} - c_n\,\sigma^2_{g \times g}$",
        save_name=f"{prefix}Vrealized_fixed_m{m_val}{gtag}",
    )


if __name__ == "__main__":
    main()
