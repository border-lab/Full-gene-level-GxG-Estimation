# -*- coding: utf-8 -*-
"""
fixed_m.py -- fixed-m, increasing-n figure for the epistasis-only MCREML results.

Usage
-----
    python fixed_m.py <result_dir> [component] [ymin ymax]

<result_dir> is a folder of two-column result files (rows "(s2gxg, s2e)") named
    <basename>_n<N>m<M>.txt
where <basename> encodes the truth, e.g.  Random_s2gxg0.05_s2e0.95.

[component] is optional and selects which variance to plot:
    both  (default) -- stacked two-panel figure (sigma^2_gxg on top, sigma^2_e below)
    gxg             -- single panel for sigma^2_gxg only
    e               -- single panel for sigma^2_e only

[ymin ymax] is an optional pair of numbers fixing the y-axis range of the
relative-error axis (e.g. -0.5 0.5).  If omitted, the range is auto-computed
from the data.  component and the y-limits may be given in any order.

Examples

    python fixed_m.py Random_s2gxg0.05_s2e0.95              # both, auto y-axis
    python fixed_m.py Random_s2gxg0.05_s2e0.95 gxg          # gxg only, auto y-axis
    python fixed_m.py Random_s2gxg0.05_s2e0.95 gxg -0.5 0.5 # gxg only, y in [-0.5, 0.5]
    python fixed_m.py Random_s2gxg0.05_s2e0.95 -0.5 0.5     # both, y in [-0.5, 0.5]

reads every n at the (single) fixed m in that folder and writes a box-plot with
sample size increasing left-to-right.  The true (s2gxg, s2e) used to centre the
relative errors is parsed from the folder name.  The PDF is written to the
current working directory as <basename>_fixed_m<M>.pdf (or, for a single
component, <basename>_fixed_m<M>_<component>.pdf).
"""
import os
import re
import sys
import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")            # headless: save PDF without a display

from plot_script import (read_MoM_results, plot_two_effects_single_m,
                         plot_relative_error_across_groups_combined)

# component name -> (column index, PDF suffix)
_COMPONENTS = {"gxg": 0, "e": 1}


def _parse_truth(basename):
    """Extract (s2gxg, s2e) truth from a folder name like Random_s2gxg0.05_s2e0.95."""
    m_gxg = re.search(r"s2gxg([0-9]*\.?[0-9]+)", basename)
    m_e = re.search(r"s2e([0-9]*\.?[0-9]+)", basename)
    if not (m_gxg and m_e):
        raise ValueError(
            f"Cannot parse true s2gxg / s2e from folder name '{basename}'. "
            "Expected it to contain 's2gxg<value>' and 's2e<value>'."
        )
    return float(m_gxg.group(1)), float(m_e.group(1))


def _discover_n_and_m(dir_path, basename):
    """Find every (n, m) present as <basename>_n<N>m<M>.txt in dir_path.

    Returns (sorted list of n, the single fixed m).  Errors if more than one m
    is present, since this figure is defined for a fixed m.
    """
    pattern = os.path.join(dir_path, f"{basename}_n*m*.txt")
    n_by_m = {}
    for f in glob.glob(pattern):
        mobj = re.search(r"_n(\d+)m(\d+)\.txt$", os.path.basename(f))
        if not mobj:
            continue
        n_val, m_val = int(mobj.group(1)), int(mobj.group(2))
        n_by_m.setdefault(m_val, []).append(n_val)

    if not n_by_m:
        raise FileNotFoundError(
            f"No files matching '{basename}_n*m*.txt' found in {dir_path}"
        )
    if len(n_by_m) > 1:
        raise ValueError(
            f"Multiple m values found in {dir_path}: {sorted(n_by_m)}. "
            "fixed_m.py expects a single (fixed) m; please separate them."
        )

    m_val = next(iter(n_by_m))
    return sorted(n_by_m[m_val]), m_val


def _auto_ylimits(data_dict, individual_sizes, real_values, cols=(0, 1), pad=0.15):
    """y-limits covering the box-whisker range of the requested component(s).

    Whiskers follow matplotlib's default 1.5*IQR rule (outliers are hidden in
    the plot), so the limits are derived from the same range and padded.
    """
    lo, hi = np.inf, -np.inf
    for col_num in cols:
        real_value = real_values[col_num]
        for n in individual_sizes:
            v = data_dict[n].iloc[:, col_num].values - real_value
            q1, q3 = np.percentile(v, [25, 75])
            iqr = q3 - q1
            w_lo = v[v >= q1 - 1.5 * iqr].min()
            w_hi = v[v <= q3 + 1.5 * iqr].max()
            lo = min(lo, w_lo)
            hi = max(hi, w_hi)
    span = hi - lo
    if span <= 0:
        span = 1.0
    return lo - pad * span, hi + pad * span


def main():
    _usage = "Usage: python fixed_m.py <result_dir> [gxg|e|both] [ymin ymax]"
    if len(sys.argv) < 2:
        print(_usage)
        sys.exit(1)

    dir_path = os.path.normpath(sys.argv[1])
    if not os.path.isdir(dir_path):
        print(f"Error: '{dir_path}' is not a directory.")
        sys.exit(1)

    # Parse the trailing tokens: a component keyword and/or a pair of numbers,
    # in any order.  Numbers are the y-axis limits (ymin ymax).
    component = "both"
    nums = []
    for tok in sys.argv[2:]:
        if tok.lower() in ("both", "gxg", "e"):
            component = tok.lower()
            continue
        try:
            nums.append(float(tok))
        except ValueError:
            print(f"Error: unrecognized argument '{tok}'.\n{_usage}")
            sys.exit(1)

    ylim = None
    if len(nums) == 2:
        ylim = (nums[0], nums[1])
        if ylim[0] >= ylim[1]:
            print(f"Error: ymin ({ylim[0]}) must be less than ymax ({ylim[1]}).")
            sys.exit(1)
    elif len(nums) != 0:
        print(f"Error: y-axis limits need exactly two numbers (ymin ymax), "
              f"got {len(nums)}.\n{_usage}")
        sys.exit(1)

    basename = os.path.basename(dir_path)
    real_gxg, real_e = _parse_truth(basename)
    individual_sizes, m_val = _discover_n_and_m(dir_path, basename)

    print(f"Directory : {dir_path}")
    print(f"Truth     : s2gxg={real_gxg}, s2e={real_e}")
    print(f"Fixed m   : {m_val}")
    print(f"n values  : {individual_sizes}")
    print(f"Component : {component}")
    print(f"y-axis    : {'auto' if ylim is None else f'[{ylim[0]}, {ylim[1]}]'}")

    # read_MoM_results reads two-column files (col 0 = gxg, col 1 = e).
    data_dict = read_MoM_results(individual_sizes, dir_path, basename, m_val)
    real_values = [real_gxg, real_e]

    if component == "both":
        if ylim is not None:
            ymin, ymax = ylim
        else:
            ymin, ymax = _auto_ylimits(data_dict, individual_sizes, real_values, cols=(0, 1))
        plot_two_effects_single_m(
            data_dict,
            individual_sizes,
            real_values=real_values,
            m_label=f"m = {m_val:,}",
            ymin=ymin,
            ymax=ymax,
            save_name=f"{basename}_fixed_m{m_val}",
        )
    else:
        # Single variance component -> single-panel box plot.
        col_num = _COMPONENTS[component]
        if ylim is not None:
            ymin, ymax = ylim
        else:
            ymin, ymax = _auto_ylimits(data_dict, individual_sizes, real_values, cols=(col_num,))
        plot_relative_error_across_groups_combined(
            data_dict,
            x_labels=[f"m = {m_val:,}"],
            individual_sizes=individual_sizes,
            col_num=col_num,
            real_value=real_values[col_num],
            ymin=ymin,
            ymax=ymax,
            x_axis_name="Sample size (n)",
            save_name=f"{basename}_fixed_m{m_val}_{component}",
        )


if __name__ == "__main__":
    main()
