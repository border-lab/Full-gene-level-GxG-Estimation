# -*- coding: utf-8 -*-
"""
gxg_col3_error_plot.py -- SIGNED error of the epistasis component, column 3
of a result file, as a single box panel in the SE-figure style.

Usage
-----
    python gxg_col3_error_plot.py <result_dir> [--ylim LO HI] [--outdir DIR]

<result_dir> holds the row-per-replicate files written by the c-normalized
4VC pipelines (see that pipeline's calc_stats.py),

    (s2a, s2d, s2gxg, s2e [, ...])
     1    2    3      4

and this figure is column 3 - sigma^2_gxg, the truth parsed from the folder
name, one box per sample size.  The signed error is used rather than a relative
one because the truth may be 0 (a null run), where dividing by it is undefined.
Keeping the sign shows bias and spread at once: the dashed line at zero is the
unbiasedness target, and the label above each box is a one-sample t-test of the
signed error against zero.

In these pipelines the lower box for s2gxg is -0.2 * the upper bound, not 0, so
negative estimates are ordinary REML output and nothing is dropped.

The figure is written as a PDF to <result_dir> (override with --outdir).
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")

from plot_script import plot_relative_error_across_groups_combined
from realized_var_four_components_plot import _parse_truth, _discover, _tag
from gxg_paired_error_plot import _auto_ylim

COL_GXG = 2             # column 3, s2gxg_hat

YLABEL = (r"Signed error "
          r"($\hat{\sigma}^2_{g \times g} - \sigma^2_{g \times g}$)")


def _read(path):
    """Parse "(a,b,...)" rows into an (R, ncol) array with at least 3 columns."""
    rows = []
    with open(path, "r", encoding="utf-8-sig") as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                rows.append([float(v) for v in line.strip("()[] ").split(",")])
            except ValueError:
                raise SystemExit(f"{path}:{lineno}: cannot parse '{line}' as numbers.")
    if not rows:
        raise SystemExit(f"{path}: no data rows.")
    widths = {len(r) for r in rows}
    if len(widths) != 1:
        raise SystemExit(f"{path}: ragged file, row widths {sorted(widths)}.")
    arr = np.array(rows)
    if arr.shape[1] < 3:
        raise SystemExit(f"{path} has {arr.shape[1]} columns, need at least 3.")
    return arr


def main():
    _usage = ("Usage: python gxg_col3_error_plot.py <result_dir> "
              "[--ylim LO HI] [--outdir DIR]")
    if len(sys.argv) < 2:
        print(_usage)
        sys.exit(1)

    dir_path = os.path.abspath(sys.argv[1])
    if not os.path.isdir(dir_path):
        print(f"Error: '{dir_path}' is not a directory.")
        sys.exit(1)

    rest = sys.argv[2:]
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

    truth = _parse_truth(os.path.basename(dir_path))["s2gxg"]
    stem, ns, m_val, g_val = _discover(dir_path)

    print(f"Directory : {dir_path}")
    print(f"Truth     : s2gxg={truth}")
    print(f"Fixed m   : {m_val}" + ("" if g_val is None else f", G = {g_val}"))
    print("Plotted   : column 3 - s2gxg (signed), per replicate")

    err = {}
    for n in ns:
        arr = _read(os.path.join(dir_path, _tag(stem, n, m_val, g_val) + ".txt"))
        e = arr[:, COL_GXG] - truth
        err[n] = e
        print(f"  n={n:<6d} R={e.size:<4d} mean={e.mean():+.6f} "
              f"SD={e.std(ddof=1):.6f} median={np.median(e):+.6f}")

    if ylim is None:
        ylim = _auto_ylim(err, ns)

    # The truth is already subtracted in err, so the helper's scalar is 0.
    data_dict = {n: pd.DataFrame({"gxg": err[n]}) for n in ns}
    gsuf = "" if g_val is None else f"_G{g_val}"
    base = f"gxg_col3_signed_error_m{m_val}{gsuf}"

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


if __name__ == "__main__":
    main()
