# -*- coding: utf-8 -*-
"""
estimate_plot.py -- fixed-m, increasing-n box plot of the *raw* sigma^2_gxg
ESTIMATE (not the relative error) for a folder of two-column MCREML result
files, written into result_figure/.

Same journal box-plot style as fixed_m.py / se_plot.py, but the y-axis is the
estimate sigma^2_gxg itself instead of (estimate - truth).  A dashed reference
line is drawn at the true sigma^2_gxg (parsed from the folder name), and the
per-box "Mean"/"SD" annotations plus a one-sample t-test against the truth
(bias test) are shown.

The result files in these folders carry the n/m tokens in the MIDDLE of the
name, e.g.

    RandomSNP_s2gxg0.2_s2e0.8_n1000m1000_G2_r0.5-0.5_est2.txt

(prefix "RandomSNP_s2gxg0.2_s2e0.8", suffix "G2_r0.5-0.5_est2"), so the file
grid is auto-detected rather than assuming "<foldername>_n*m*.txt".

    python estimate_plot.py RandomSNP_s2gxg0.2_s2e0.8_G2_r0.5-0.5_est2 [more folders...]
"""
import os
import re
import sys
import glob

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")            # headless: save PDF without a display
import matplotlib.pyplot as plt
from scipy import stats

from fixed_m import _parse_truth

HERE = os.path.dirname(os.path.abspath(__file__))
RESULT_ROOT = os.path.join(HERE, "estimation_result")
OUT_FOLDER = os.path.join(HERE, "result_figure")


def discover_grid(folder_path):
    """Return (prefix, suffix, sorted n list, single m) from the actual result
    files, where a file is named <prefix>_n<N>m<M>[_<suffix>].txt."""
    n_by_m = {}
    keys = set()
    for f in glob.glob(os.path.join(folder_path, "*_n*m*.txt")):
        b = os.path.basename(f)
        mobj = re.search(r"^(.*)_n(\d+)m(\d+)(.*)\.txt$", b)
        if not mobj:
            continue
        prefix, n_val, m_val, suffix = (
            mobj.group(1), int(mobj.group(2)), int(mobj.group(3)), mobj.group(4)
        )
        keys.add((prefix, suffix))
        n_by_m.setdefault(m_val, []).append(n_val)

    if not n_by_m:
        raise FileNotFoundError(f"No *_n<N>m<M>*.txt result files in {folder_path}")
    if len(keys) != 1:
        raise ValueError(f"Expected one (prefix, suffix) in {folder_path}, got {keys}")
    if len(n_by_m) != 1:
        raise ValueError(f"Expected a single (fixed) m in {folder_path}, got {sorted(n_by_m)}")

    prefix, suffix = keys.pop()
    m_val = next(iter(n_by_m))
    return prefix, suffix, sorted(n_by_m[m_val]), m_val


def read_estimates(folder_path, prefix, suffix, ns, m_val):
    """{n -> 1-D array of sigma^2_gxg estimates} read from column 0."""
    out = {}
    for n in ns:
        fp = os.path.join(folder_path, f"{prefix}_n{n}m{m_val}{suffix}.txt")
        df = pd.read_csv(fp, header=None)
        out[n] = df[0].astype(str).str.replace("(", "", regex=False).astype(float).values
    return out


def _auto_ylimits(estimates, ns, ref, pad=0.15):
    """y-limits covering the box-whisker (1.5*IQR) range of the raw estimates and
    always including the reference line."""
    lo, hi = ref, ref
    for n in ns:
        v = estimates[n]
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


def plot_estimate(estimates, ns, m_val, truth, save_name, title=None, ref_line=None):
    # Dashed reference line; defaults to the truth but may be set independently.
    if ref_line is None:
        ref_line = truth
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 10,
        'axes.linewidth': 1,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.major.width': 1,
        'ytick.major.width': 1,
        'xtick.major.size': 4,
        'ytick.major.size': 4,
        'figure.dpi': 150,
        'savefig.dpi': 600,
    })

    labels = [f"n = {n:,}" for n in ns]
    box_data = [estimates[n] for n in ns]
    x_positions = np.arange(len(ns))

    ymin, ymax = _auto_ylimits(estimates, ns, ref_line)
    ymax_extended = ymax + 0.28 * (ymax - ymin)

    fig, ax = plt.subplots(figsize=(3 + len(ns) * 0.8, 5))
    ax.set_ylim(ymin, ymax_extended)
    ax.set_xlim(-0.6, len(ns) - 0.4)

    # Alternating white / gray background for each box
    for i in range(len(ns)):
        ax.axvspan(i - 0.5, i + 0.5,
                   facecolor='white' if i % 2 == 0 else '#E8E8E8',
                   alpha=1.0 if i % 2 == 0 else 0.8, zorder=0)

    box_color = '#3274A1'
    median_color = '#CC0000'
    label_color = '#000000'

    ax.boxplot(
        box_data,
        positions=x_positions,
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        boxprops=dict(linewidth=1.5, edgecolor=box_color, facecolor='white'),
        whiskerprops=dict(linewidth=1.2, color=box_color),
        capprops=dict(linewidth=1.2, color=box_color),
        medianprops=dict(linewidth=2, color=median_color),
    )

    # Top group label (fixed m)
    ax.text((len(ns) - 1) / 2.0, ymax_extended - 0.01 * (ymax_extended - ymin),
            f"m = {m_val:,}", ha='center', va='top',
            fontsize=11, fontweight='bold', color=label_color)

    # Per-box significance (bias vs truth), Mean, SD
    for i, n in enumerate(ns):
        v = estimates[n]
        _, p_val = stats.ttest_1samp(v, truth)
        if p_val < 0.001:
            sig = '***'
        elif p_val < 0.01:
            sig = '**'
        elif p_val < 0.05:
            sig = '*'
        else:
            sig = 'ns'
        top = ymax_extended
        ax.text(i, top - 0.09 * (top - ymin), sig, ha='center', va='top',
                fontsize=11, fontweight='bold', color=label_color)
        ax.text(i, top - 0.16 * (top - ymin), f"Mean={v.mean():.3f}",
                ha='center', va='top', fontsize=9, fontweight='bold', color=label_color)
        ax.text(i, top - 0.23 * (top - ymin), f"SD={v.std(ddof=1):.3f}",
                ha='center', va='top', fontsize=9, fontweight='bold', color=label_color)

    # Dashed reference line (defaults to the truth; overridable via ref_line)
    ax.axhline(ref_line, color='#666666', linestyle='--', linewidth=0.8, zorder=1)

    ax.set_xlabel("Sample size (n)", fontsize=10, labelpad=8)
    ax.set_ylabel(r"Estimate ($\sigma^2_{g \times g}$)", fontsize=10, labelpad=8)

    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels, fontsize=9)

    if title:
        ax.set_title(title, fontsize=11, pad=10, fontweight='normal')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#333333')
    ax.spines['bottom'].set_color('#333333')
    ax.tick_params(axis='both', which='major', labelsize=9, colors='#333333')

    plt.tight_layout()
    os.makedirs(OUT_FOLDER, exist_ok=True)
    save_path = os.path.join(OUT_FOLDER, f"{save_name}.pdf")
    plt.savefig(save_path, bbox_inches='tight', facecolor='white',
                edgecolor='none', format='pdf', transparent=False, pad_inches=0.1)
    plt.rcParams.update(plt.rcParamsDefault)
    return save_path, (ymin, ymax)


def main():
    args = sys.argv[1:]
    ref_line = None                      # dashed line value; None -> use truth
    if "--ref" in args:
        i = args.index("--ref")
        ref_line = float(args[i + 1])
        del args[i:i + 2]

    if not args:
        print("Usage: python estimate_plot.py [--ref <value>] <result_dir> [<result_dir> ...]")
        sys.exit(1)

    for arg in args:
        folder = os.path.basename(os.path.normpath(arg))
        # Accept either a bare folder name or a full/relative path.
        candidates = [arg, os.path.join(RESULT_ROOT, folder), os.path.join(HERE, folder)]
        folder_path = next((c for c in candidates if os.path.isdir(c)), None)
        if folder_path is None:
            print(f"Error: could not locate folder '{arg}'.")
            sys.exit(1)

        real_gxg, real_e = _parse_truth(folder)
        prefix, suffix, ns, m_val = discover_grid(folder_path)

        print(f"\nFolder    : {folder_path}")
        print(f"Truth     : s2gxg={real_gxg}, s2e={real_e}")
        print(f"Fixed m   : {m_val}")
        print(f"n values  : {ns}")

        estimates = read_estimates(folder_path, prefix, suffix, ns, m_val)
        save_path, (ymin, ymax) = plot_estimate(
            estimates, ns, m_val, real_gxg,
            save_name=f"{folder}_fixed_m{m_val}_gxg_estimate",
            ref_line=ref_line,
        )
        print(f"Saved     : {save_path}  (y in [{ymin:.3f}, {ymax:.3f}])")


if __name__ == "__main__":
    main()
