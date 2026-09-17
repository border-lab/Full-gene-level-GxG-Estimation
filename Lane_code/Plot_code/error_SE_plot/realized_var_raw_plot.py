# -*- coding: utf-8 -*-
"""
realized_var_raw_plot.py -- fixed-m, increasing-n figure of the THIRD column on
its own: the estimated realized variance V_ell-hat = c-hat * V_gamma-hat, in
absolute units, with no reference and no comparison.

Usage
-----
    python realized_var_raw_plot.py <result_dir> [ymin ymax] [--outdir DIR]

Every other script in this family centres col3 on something -- c_hat * s2gxg
(realized_var_plot.py), the per-replicate col4 (realized_var_perrep_plot.py), or
mean(col4) (realized_var_vs_empirical_plot.py).  This one does not: it shows the
raw distribution of V_ell-hat at each n so the absolute scale and the shape of
the estimator's sampling distribution are visible directly.

Because there is no reference point there is no null hypothesis, so the
significance stars and the dashed zero line are both switched off -- a t-test of
a variance against 0 would print *** everywhere and mean nothing.  The mean and
SD annotations are kept; they are descriptive, not comparative.

Default output directory is <result_dir> itself (override with --outdir).
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")

from plot_script import plot_relative_error_across_groups_combined

from realized_var_perrep_plot import _parse_truth, _discover, _tag, _read


def main():
    _usage = ("Usage: python realized_var_raw_plot.py <result_dir> "
              "[ymin ymax] [--outdir DIR]")
    if len(sys.argv) < 2:
        print(_usage)
        sys.exit(1)

    dir_path = os.path.abspath(sys.argv[1])
    if not os.path.isdir(dir_path):
        print(f"Error: '{dir_path}' is not a directory.")
        sys.exit(1)

    out_dir = dir_path
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
    print("Plotted   : col3 (c_hat * V_gamma_hat), raw -- no reference")

    raw = {}
    for n in ns:
        v = data[n][2].values
        raw[n] = pd.DataFrame({0: v})
        q = np.percentile(v, [25, 50, 75])
        print(f"  n={n:<6d} R={len(v):<4d} mean={v.mean():.6f}  SD={v.std(ddof=1):.6f}  "
              f"min={v.min():.4f}  q1={q[0]:.4f}  med={q[1]:.4f}  q3={q[2]:.4f}  max={v.max():.4f}")

    if ymin is None:
        lo, hi = np.inf, -np.inf
        for n in ns:
            v = raw[n][0].values
            q1, q3 = np.percentile(v, [25, 75])
            iqr = q3 - q1
            lo = min(lo, v[v >= q1 - 1.5 * iqr].min())
            hi = max(hi, v[v <= q3 + 1.5 * iqr].max())
        span = (hi - lo) or 1.0
        ymin, ymax = lo - 0.15 * span, hi + 0.15 * span
    print(f"y-axis    : [{ymin}, {ymax}]")

    os.makedirs(out_dir, exist_ok=True)
    os.chdir(out_dir)

    gtag = "" if g_val is None else f"_G{g_val}"
    prefix = "" if os.path.normcase(out_dir) == os.path.normcase(dir_path) else folder + "_"
    plot_relative_error_across_groups_combined(
        raw,
        x_labels=[f"m = {m_val:,}"],
        individual_sizes=ns,
        col_num=0,
        real_value=0.0,              # raw values, nothing subtracted
        ymin=ymin,
        ymax=ymax,
        x_axis_name="Sample size (n)",
        y_axis_name=r"$\hat{c}\,\hat{V}_{\gamma}$",
        show_significance=False,     # no reference -> no null hypothesis
        ref_line=None,               # and no dashed zero line
        save_name=f"{prefix}Vl_fixed_m{m_val}{gtag}_raw",
    )


if __name__ == "__main__":
    main()
