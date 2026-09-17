# -*- coding: utf-8 -*-
"""
realized_var_vs_empirical_plot.py -- fixed-m, increasing-n figure validating the
c construction: col3 (V_ell-hat = c-hat * V_gamma-hat) referenced against the
EMPIRICAL mean of col4, the simulated realized variance.

Usage
-----
    python realized_var_vs_empirical_plot.py <result_dir> [ymin ymax] [--outdir DIR]

Plotted per n:   V_ell-hat_r  -  mean_r(V_ell^realized)

Why the empirical mean and not c_exact * s2gxg
----------------------------------------------
realized_var_plot.py --exact centres on the theoretical E[V_ell^realized] =
c_exact * s2gxg.  That reference ASSUMES the c construction is right, so it can
only ever measure the estimator given c.  Referencing the simulated col4 mean
instead makes no such assumption: the whole chain -- the HWE closed form for
c-hat, the theory behind it, and the variance-component estimate -- is tested
against the realized variance the simulator actually produced.  That is the
figure to show when the claim is "our c is correct".

The price is that the reference is estimated, with SE = SD(col4)/sqrt(R), and
here that SE (~0.008-0.010) is comparable to the mean errors being judged.
Two consequences, both handled:

  * The BOX HEIGHTS are SD(col3), the spread of the estimator about a fixed
    reference -- not the SD of the error.  They shrink with n.  For the honest
    per-replicate error spread see realized_var_perrep_plot.py.
  * The STARS cannot come from a one-sample t-test of the plotted values: that
    would use SD(col3) as the standard error and ignore the sampling error of
    the reference, over-reporting significance.  col3 and col4 share the
    replicate, so the correct test of mean(col3) == mean(col4) is the PAIRED
    t-test on col3 - col4; those p-values are computed here and passed to the
    plotter.  Note mean(col3 - mean(col4)) == mean(col3 - col4) exactly, so the
    means annotated on this figure match realized_var_perrep_plot.py.

Default output directory is <result_dir> itself (override with --outdir).
"""
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
from scipy import stats

from plot_script import plot_relative_error_across_groups_combined

from realized_var_perrep_plot import _parse_truth, _discover, _tag, _read


def main():
    _usage = ("Usage: python realized_var_vs_empirical_plot.py <result_dir> "
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
    print("Reference : mean(col4), the empirical mean simulated realized variance")
    print("Stars     : paired t-test of col3 - col4 (reference is estimated)")

    centred, pvals = {}, []
    for n in ns:
        df = data[n]
        vh, vr = df[2].values, df[3].values
        ref = vr.mean()
        se_ref = vr.std(ddof=1) / np.sqrt(len(vr))
        centred[n] = pd.DataFrame({0: vh - ref})

        # Reference estimated from the same replicates -> paired test.
        t_stat, p_paired = stats.ttest_rel(vh, vr)
        pvals.append(p_paired)

        print(f"  n={n:<6d} R={len(vh):<4d} mean(col4)={ref:.6f} (SE {se_ref:.6f})  "
              f"mean(col3)={vh.mean():.6f}")
        print(f"           offset={vh.mean() - ref:+.6f}  paired p={p_paired:.4f}  "
              f"| box SD=SD(col3)={vh.std(ddof=1):.6f}  "
              f"err SD=SD(col3-col4)={(vh - vr).std(ddof=1):.6f}")

    if ymin is None:
        lo, hi = np.inf, -np.inf
        for n in ns:
            v = centred[n][0].values
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
        centred,
        x_labels=[f"m = {m_val:,}"],
        individual_sizes=ns,
        col_num=0,
        real_value=0.0,
        ymin=ymin,
        ymax=ymax,
        x_axis_name="Sample size (n)",
        y_axis_name=r"$\hat{c}\,\hat{V}_{\gamma} - \overline{V_{\ell}^{\mathrm{realized}}}$",
        p_values=pvals,
        save_name=f"{prefix}Vl_fixed_m{m_val}{gtag}_vs_meanVreal",
    )


if __name__ == "__main__":
    main()
