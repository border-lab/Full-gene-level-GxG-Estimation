# -*- coding: utf-8 -*-
"""
gxg_col3_error_scaling_plot.py -- SIGNED error of the epistasis component
(column 3 - sigma^2_gxg) when n, m and G are scaled TOGETHER.

The sibling gxg_col3_error_plot.py requires one (m, G) setting and a sweep
over n.  This script is for folders where every file has its own (n, m, G),
e.g. n/m fixed at 2 with G = m/100:

    ..._n2000m1000_G10.txt   ..._n4000m2000_G20.txt   ...

Usage
-----
    python gxg_col3_error_scaling_plot.py <result_dir> [--component gxg|a|d|e]
                                          [--relative] [--ylim LO HI] [--outdir DIR]

One box per file, ordered by n, labelled with its n, m and G.  By default
column 3 (s2gxg) is plotted as a SIGNED error, as in gxg_col3_error_plot.py.
--component picks another variance component (columns 1, 2, 4 for a, d, e),
and --relative divides the error by the truth, which must then be non-zero.
"""
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")

from plot_script import plot_relative_error_across_groups_combined
from realized_var_four_components_plot import _parse_truth
from gxg_paired_error_plot import _auto_ylim
from gxg_col3_error_plot import _read

# column index and label subscript for each variance component
COMPONENTS = {"a": (0, "a"), "d": (1, "d"), "gxg": (2, r"g \times g"), "e": (3, "e")}


def _discover_triples(dir_path, gp=None):
    """Every (n, m, G) present, one file each, sorted by n.  gp=None keeps only
    files without a _gp<pct> suffix; gp="3.125" keeps only that gene_pct."""
    found = []
    for f in sorted(glob.glob(os.path.join(dir_path, "*_n*m*.txt"))):
        b = os.path.basename(f)
        if b.startswith("c_"):
            continue
        mobj = re.search(r"_n(\d+)m(\d+)(?:_G(\d+))?(?:_gp([\d.]+))?\.txt$", b)
        if not mobj or mobj.group(4) != gp:
            continue
        g = None if mobj.group(3) is None else int(mobj.group(3))
        found.append((int(mobj.group(1)), int(mobj.group(2)), g, f))
    if not found:
        raise FileNotFoundError(f"No result files found in {dir_path}")
    found.sort()
    ns = [t[0] for t in found]
    if len(set(ns)) != len(ns):
        raise ValueError(f"Duplicate sample sizes in {dir_path}: {ns}")
    return found


def main():
    _usage = ("Usage: python gxg_col3_error_scaling_plot.py <result_dir> "
              "[--component gxg|a|d|e] [--relative] [--ylim LO HI] [--outdir DIR]")
    if len(sys.argv) < 2:
        print(_usage)
        sys.exit(1)
    dir_path = os.path.abspath(sys.argv[1])
    if not os.path.isdir(dir_path):
        print(f"Error: '{dir_path}' is not a directory.")
        sys.exit(1)

    rest = sys.argv[2:]
    comp = "gxg"
    if "--component" in rest:
        i = rest.index("--component")
        if i + 1 >= len(rest) or rest[i + 1] not in COMPONENTS:
            print(f"Error: --component needs one of {sorted(COMPONENTS)}.")
            sys.exit(1)
        comp = rest[i + 1]
        rest = rest[:i] + rest[i + 2:]
    relative = "--relative" in rest
    if relative:
        rest.remove("--relative")
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

    truth = _parse_truth(os.path.basename(dir_path))["s2" + comp]
    if relative and truth == 0:
        raise SystemExit(f"s2{comp} = 0 in this folder; a relative error is undefined. "
                         "Drop --relative.")
    col, sub = COMPONENTS[comp]
    triples = _discover_triples(dir_path)
    ns = [t[0] for t in triples]
    kind = "relative" if relative else "signed"
    hat, sig = rf"\hat{{\sigma}}^2_{{{sub}}}", rf"\sigma^2_{{{sub}}}"
    ylabel = (rf"Relative error $({hat} - {sig}) / {sig}$" if relative
              else rf"Signed error $({hat} - {sig})$")

    print(f"Directory : {dir_path}")
    print(f"Truth     : s2{comp}={truth}")
    print(f"Plotted   : column {col + 1} - s2{comp} ({kind}), per replicate; "
          "n, m, G scaled together")

    err, labels = {}, []
    for n, m, g, path in triples:
        e = _read(path)[:, col] - truth
        if relative:
            e = e / truth
        err[n] = e
        labels.append(f"n = {n:,}\nm = {m:,}" + ("" if g is None else f"\nG = {g}"))
        print(f"  n={n:<6d} m={m:<6d} G={g!s:<4} R={e.size:<4d} mean={e.mean():+.6f} "
              f"SD={e.std(ddof=1):.6f} median={np.median(e):+.6f}")

    if ylim is None:
        ylim = _auto_ylim(err, ns)

    ratios = {round(n / m, 3) for n, m, _, _ in triples}
    top = (f"n / m = {ratios.pop():g}" if len(ratios) == 1 else "n, m, G scaled together")
    data_dict = {n: pd.DataFrame({"gxg": err[n]}) for n in ns}
    base = f"{comp}_col{col + 1}_{kind}_error_nmG_scaling"

    os.makedirs(out_dir, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(out_dir)      # the plot helper writes to the CWD
    try:
        plot_relative_error_across_groups_combined(
            data_dict,
            x_labels=[top],
            individual_sizes=ns,
            col_num=0,
            real_value=0.0,
            ymin=ylim[0],
            ymax=ylim[1],
            x_axis_name="Sample size (n), SNPs (m), genes (G)",
            y_axis_name=ylabel,
            custom_bottom_labels=labels,
            save_name=base,
        )
    finally:
        os.chdir(cwd)
    print(f"Saved     : {os.path.join(out_dir, base + '.pdf')}  "
          f"(y in [{ylim[0]:.3f}, {ylim[1]:.3f}])")


if __name__ == "__main__":
    main()
