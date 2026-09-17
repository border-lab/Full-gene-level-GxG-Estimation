# -*- coding: utf-8 -*-
"""
realized_var_perrep_plot.py -- fixed-m, increasing-n figure for the FOUR-column
result files, where the realized variance is known per replicate.

Usage
-----
    python realized_var_perrep_plot.py <result_dir> [ymin ymax] [--outdir DIR]

<result_dir> holds files named <basename>_n<N>m<M>_G<G>.txt whose rows are
"(s2gxg_hat, s2e_hat, Vl_hat, Vl_realized)":

    col 3  Vl_hat       = c-hat * s2gxg-hat, the estimated realized variance
    col 4  Vl_realized  = the simulated realized variance Var(H gamma) of THAT
                          replicate's own gamma draw

Unlike realized_var_plot.py -- which centres column 3 on the per-n constant
E[Vl] = c * s2gxg because no per-replicate truth was recorded -- here the truth
is available replicate by replicate, so the error plotted is simply

    error_r = Vl_hat_r - Vl_realized_r

paired within each replicate.  This removes the Monte-Carlo scatter of gamma
itself from the error, so the boxes show only the estimator's error and not the
spread of the realized variance around its expectation.  The reference line is
therefore 0 for every n and the boxes are directly comparable across n.

The PDF goes to result_figure/ next to this script (override with --outdir).
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


def _parse_truth(basename):
    m_gxg = re.search(r"s2gxg([0-9]*\.?[0-9]+)", basename)
    m_e = re.search(r"s2e([0-9]*\.?[0-9]+)", basename)
    if not (m_gxg and m_e):
        raise ValueError(f"Cannot parse s2gxg / s2e from '{basename}'.")
    return float(m_gxg.group(1)), float(m_e.group(1))


def _discover(dir_path):
    """Find the file stem and every (n, m, G) present in dir_path."""
    grid = {}
    stem = None
    for f in sorted(glob.glob(os.path.join(dir_path, "*_n*m*.txt"))):
        b = os.path.basename(f)
        if b.startswith("c_"):
            continue
        mobj = re.search(r"^(.*)_n(\d+)m(\d+)(?:_G(\d+))?\.txt$", b)
        if not mobj:
            continue
        stem = mobj.group(1)
        grid.setdefault((int(mobj.group(3)), mobj.group(4)), []).append(int(mobj.group(2)))
    if not grid:
        raise FileNotFoundError(f"No result files found in {dir_path}")
    if len(grid) != 1:
        raise ValueError(f"Expected one (m, G) setting, found {sorted(grid)}.")
    (m_val, g_val), ns = next(iter(grid.items()))
    return stem, sorted(ns), m_val, g_val


def _tag(stem, n, m_val, g_val):
    suffix = "" if g_val is None else f"_G{g_val}"
    return f"{stem}_n{n}m{m_val}{suffix}"


def _read(dir_path, stem, ns, m_val, g_val):
    """Read the four-column result files into {n: DataFrame}."""
    out = {}
    for n in ns:
        path = os.path.join(dir_path, _tag(stem, n, m_val, g_val) + ".txt")
        df = pd.read_csv(path, header=None)
        if df.shape[1] != 4:
            raise ValueError(
                f"{path} has {df.shape[1]} columns, need 4 "
                "(s2gxg_hat, s2e_hat, Vl_hat, Vl_realized). "
                "For three-column files use realized_var_plot.py instead.")
        df[0] = df[0].astype(str).str.replace("(", "", regex=False).astype(float)
        df[1] = df[1].astype(float)
        df[2] = df[2].astype(float)
        df[3] = df[3].astype(str).str.replace(")", "", regex=False).astype(float)
        out[n] = df
    return out


def main():
    _usage = ("Usage: python realized_var_perrep_plot.py <result_dir> "
              "[ymin ymax] [--outdir DIR]")
    if len(sys.argv) < 2:
        print(_usage)
        sys.exit(1)

    dir_path = os.path.abspath(sys.argv[1])
    if not os.path.isdir(dir_path):
        print(f"Error: '{dir_path}' is not a directory.")
        sys.exit(1)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(script_dir, "result_figure")
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
    print("Error     : col3 (Vl_hat) - col4 (Vl_realized), paired per replicate")

    centred = {}
    for n in ns:
        df = data[n]
        err = df[2].values - df[3].values
        centred[n] = pd.DataFrame({0: err})
        print(f"  n={n:<6d} R={len(err):<4d} "
              f"mean(Vl_hat)={df[2].mean():.6f}  mean(Vl_realized)={df[3].mean():.6f}")
        print(f"           error mean={err.mean():+.6f}  SD={err.std(ddof=1):.6f}  "
              f"RMSE={np.sqrt(np.mean(err ** 2)):.6f}")

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

    # plot_relative_error_across_groups_combined writes its PDF to the CWD.
    os.makedirs(out_dir, exist_ok=True)
    os.chdir(out_dir)

    gtag = "" if g_val is None else f"_G{g_val}"
    plot_relative_error_across_groups_combined(
        centred,
        x_labels=[f"m = {m_val:,}"],
        individual_sizes=ns,
        col_num=0,
        real_value=0.0,          # errors are already differences against the truth
        ymin=ymin,
        ymax=ymax,
        x_axis_name="Sample size (n)",
        y_axis_name=r"Relative error ($\hat{V}_{\ell} - V_{\ell}^{\mathrm{realized}}$)",
        save_name=f"{folder}_fixed_m{m_val}{gtag}_Vl_perrep",
    )


if __name__ == "__main__":
    main()
