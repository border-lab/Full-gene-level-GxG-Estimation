# -*- coding: utf-8 -*-
"""
realized_var_plot.py -- fixed-m, increasing-n figure for the THIRD column of the
three-column result files, V_ell-hat = c-hat * s2gxg-hat, the c-corrected
estimate of the realized variance Var-hat(H gamma) (realized_variance.pdf).

Usage
-----
    python realized_var_plot.py <result_dir> [ymin ymax] [--exact] [--outdir DIR]

<result_dir> holds files named <basename>_n<N>m<M>_G<G>.txt whose rows are
"(s2gxg_hat, s2e_hat, Vl_hat)".  The truth this column is centred on is
E[Var-hat(H gamma)] = c * s2gxg, and c is genotype-only -- one constant per n,
so the data are centred per n and the shared reference line is 0.

c per n comes from the Cholesky job c_<basename>_n<N>m<M>_G<G>.txt file when
that is present: c_hat_hwe (the O(nm) HWE closed form the estimator actually
applied) by default, c_exact (mean sample variance of the P interaction
columns) with --exact.  Both are printed either way, so the HWE approximation
error is visible next to the truncation bias.  Without a c file, c is recovered
from the results themselves: Vl_hat = c * s2gxg_hat exactly, so the ratio
Vl_hat / s2gxg_hat is constant within a file and equals c_hat_hwe.

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
    """Find the file stem and every (n, m, G) present in dir_path.

    c_<...>.txt files sit in the same folder and are skipped here; they are
    read by _read_c.
    """
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
    """Read the three-column result files into {n: DataFrame}."""
    out = {}
    for n in ns:
        path = os.path.join(dir_path, _tag(stem, n, m_val, g_val) + ".txt")
        df = pd.read_csv(path, header=None)
        if df.shape[1] < 3:
            raise ValueError(f"{path} has {df.shape[1]} columns, need 3.")
        df[0] = df[0].astype(str).str.replace("(", "", regex=False).astype(float)
        df[1] = df[1].astype(float)
        df[2] = df[2].astype(str).str.replace(")", "", regex=False).astype(float)
        out[n] = df
    return out


def _read_c(dir_path, stem, n, m_val, g_val):
    """Return {key: value} from c_<tag>.txt, or None if the file is absent."""
    path = os.path.join(dir_path, "c_" + _tag(stem, n, m_val, g_val) + ".txt")
    if not os.path.isfile(path):
        return None
    vals = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) == 2:
                vals[parts[0]] = float(parts[1])
    return vals


def main():
    _usage = ("Usage: python realized_var_plot.py <result_dir> [ymin ymax] "
              "[--exact] [--outdir DIR]")
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

    use_exact = "--exact" in rest
    rest = [t for t in rest if t != "--exact"]

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

    which = "c_exact" if use_exact else "c_hat_hwe"
    print(f"Directory : {dir_path}")
    print(f"Truth     : s2gxg={real_gxg}, s2e={real_e}")
    print(f"Fixed m   : {m_val}" + ("" if g_val is None else f", G = {g_val}"))
    print(f"c used    : {which}")

    centred = {}
    for n in ns:
        df = data[n]

        # c the estimator applied, recoverable from the results alone.
        ratio = (df[2] / df[0]).values
        c_ratio = float(np.median(ratio))
        spread = float(np.max(np.abs(ratio - c_ratio)))
        if spread > 1e-8 * max(1.0, abs(c_ratio)):
            raise ValueError(f"n={n}: Vl/s2gxg is not constant (spread {spread:.2e}).")

        cvals = _read_c(dir_path, stem, n, m_val, g_val)
        if cvals is None:
            c_hwe, c_exact, src = c_ratio, None, "from results"
        else:
            c_hwe = cvals.get("c_hat_hwe", c_ratio)
            c_exact = cvals.get("c_exact")
            src = "c file"
            if abs(c_hwe - c_ratio) > 1e-6 * max(1.0, abs(c_hwe)):
                raise ValueError(
                    f"n={n}: c_hat_hwe in the c file ({c_hwe:.10f}) disagrees with "
                    f"Vl_hat/s2gxg_hat in the results ({c_ratio:.10f}); the c file "
                    "does not belong to these results.")

        if use_exact:
            if c_exact is None:
                raise ValueError(f"n={n}: --exact needs c_exact, but no c file was found.")
            c_n = c_exact
        else:
            c_n = c_hwe

        truth_n = c_n * real_gxg
        centred[n] = pd.DataFrame({0: df[2].values - truth_n})
        extra = "" if c_exact is None else f"  c_exact={c_exact:.6f} (truth {c_exact * real_gxg:.6f})"
        print(f"  n={n:<6d} [{src}] c_hat_hwe={c_hwe:.6f}{extra}")
        print(f"           truth = c*s2gxg = {truth_n:.6f}  mean(Vl_hat)={df[2].mean():.6f}  "
              f"error mean={centred[n][0].mean():+.6f}  SD={centred[n][0].std():.6f}")

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
    etag = "_cexact" if use_exact else ""
    # Writing beside the data: the folder already names the setting, and
    # repeating it in the filename overruns the Windows 260-char path limit.
    prefix = "" if os.path.normcase(out_dir) == os.path.normcase(dir_path) else folder + "_"
    plot_relative_error_across_groups_combined(
        centred,
        x_labels=[f"m = {m_val:,}"],
        individual_sizes=ns,
        col_num=0,
        real_value=0.0,          # data are already centred on the per-n truth
        ymin=ymin,
        ymax=ymax,
        x_axis_name="Sample size (n)",
        y_axis_name=r"$\hat{c}\,\hat{V}_{\gamma} - \mathbb{E}[V_{\ell}^{\mathrm{realized}}]$",
        save_name=f"{prefix}Vl_fixed_m{m_val}{gtag}{etag}",
    )


if __name__ == "__main__":
    main()
