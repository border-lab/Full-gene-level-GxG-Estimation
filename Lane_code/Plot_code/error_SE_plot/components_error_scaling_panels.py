# -*- coding: utf-8 -*-
"""
components_error_scaling_panels.py -- one figure, one panel per variance
component stacked in a column, for a folder where n, m and G are scaled together.

Usage
-----
    python components_error_scaling_panels.py <result_dir>
        [--panels a,d,gxg] [--relative a,d,gxg] [--gp PCT] [--outdir DIR]

Each panel is the box-per-setting figure of gxg_col3_error_scaling_plot.py:
error of the component over replicates, dashed line at zero, one-sample
t-test label and Mean / SD above each box.  --panels lists the components in
order (default a,d,gxg); --relative names those shown as a relative error
(default: all panels; pass e.g. --relative a,d in a null folder, where the
gxg truth is zero and its relative error is undefined, to show gxg signed).  Every panel
has its own y-axis, symmetric about zero and covering the whiskers.
--gp picks the *_gp<PCT>.txt files (W from the first PCT% of SNPs); without
it only files with no _gp suffix are used.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from plot_script import _nice_tick_step
from realized_var_four_components_plot import _parse_truth
from gxg_paired_error_plot import _auto_ylim
from gxg_col3_error_plot import _read
from gxg_col3_error_scaling_plot import _discover_triples, COMPONENTS

BOX_COLOR, MEDIAN_COLOR = "#3274A1", "#CC0000"
NAMES = {"a": "Additive", "d": "Dominance", "gxg": "Epistatic", "e": "Residual"}


def _stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def _panel(ax, err, ns, ylabel, title):
    ymin, ymax = _auto_ylim(err, ns)
    ytop = ymax + 0.25 * (ymax - ymin)
    ax.set_ylim(ymin, ytop)
    ax.set_xlim(-0.6, len(ns) - 0.4)
    for i in range(len(ns)):
        if i % 2:
            ax.axvspan(i - 0.5, i + 0.5, facecolor="#E8E8E8", alpha=0.8, zorder=0)
    ax.boxplot([err[n] for n in ns], positions=range(len(ns)), widths=0.5,
               patch_artist=True, showfliers=False,
               boxprops=dict(linewidth=1.5, edgecolor=BOX_COLOR, facecolor="white"),
               whiskerprops=dict(linewidth=1.2, color=BOX_COLOR),
               capprops=dict(linewidth=1.2, color=BOX_COLOR),
               medianprops=dict(linewidth=2, color=MEDIAN_COLOR))
    ax.axhline(0, color="#666666", linestyle="--", linewidth=0.8, zorder=1)
    span = ytop - ymin
    for i, n in enumerate(ns):
        v = err[n]
        p = stats.ttest_1samp(v, 0).pvalue
        ax.text(i, ytop - 0.08 * span, _stars(p), ha="center", va="top",
                fontsize=11, fontweight="bold")
        ax.text(i, ytop - 0.15 * span, f"Mean={v.mean():.3f}", ha="center",
                va="top", fontsize=8.5, fontweight="bold")
        ax.text(i, ytop - 0.22 * span, f"SD={v.std(ddof=1):.3f}", ha="center",
                va="top", fontsize=8.5, fontweight="bold")
    step = _nice_tick_step(ymax - ymin)
    k0, k1 = int(np.ceil(ymin / step)), int(np.floor(ymax / step))
    ax.set_yticks([round(k * step, 10) for k in range(k0, k1 + 1)])
    ax.set_ylabel(ylabel, fontsize=10, labelpad=6)
    ax.set_title(title, fontsize=10, loc="left", pad=6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="both", which="major", labelsize=9, colors="#333333")


def _opt(args, flag, default):
    if flag in args:
        i = args.index(flag)
        if i + 1 >= len(args):
            raise SystemExit(f"{flag} needs a value.")
        val = args[i + 1]
        del args[i:i + 2]
        return val
    return default


def main():
    args = sys.argv[1:]
    panels = _opt(args, "--panels", "a,d,gxg").split(",")
    relative = set(_opt(args, "--relative", ",".join(panels)).split(","))
    out_dir = _opt(args, "--outdir", None)
    gp = _opt(args, "--gp", None)
    if len(args) != 1:
        print(__doc__)
        sys.exit(1)
    dir_path = os.path.abspath(args[0])
    if not os.path.isdir(dir_path):
        raise SystemExit(f"'{dir_path}' is not a directory.")
    bad = [c for c in panels if c not in COMPONENTS]
    if bad:
        raise SystemExit(f"Unknown component(s) {bad}; choose from {sorted(COMPONENTS)}.")
    out_dir = os.path.abspath(out_dir) if out_dir else dir_path

    truths = _parse_truth(os.path.basename(dir_path))
    triples = _discover_triples(dir_path, gp)
    ns = [t[0] for t in triples]
    arrs = {n: _read(path) for n, _, _, path in triples}
    labels = [f"n = {n:,}\nm = {m:,}" + ("" if g is None else f"\nG = {g}")
              for n, m, g, _ in triples]
    ratios = {round(n / m, 3) for n, m, _, _ in triples}
    top = f"n / m = {ratios.pop():g}" if len(ratios) == 1 else "n, m, G scaled together"
    if gp is not None:
        top += f", W from first {gp}% of SNPs"

    plt.rcParams.update({"font.family": "Arial", "font.size": 10, "axes.linewidth": 1,
                         "xtick.major.width": 1, "ytick.major.width": 1,
                         "xtick.major.size": 4, "ytick.major.size": 4})
    # Panels stacked vertically: every panel gets the full width for its
    # x labels and the reader compares components down the page.
    fig, axes = plt.subplots(len(panels), 1,
                             figsize=(2.6 + len(ns) * 1.0, 4.2 * len(panels)))
    axes = np.atleast_1d(axes)
    print(f"Directory : {dir_path}")
    for k, (comp, ax) in enumerate(zip(panels, axes)):
        col, sub = COMPONENTS[comp]
        truth = truths["s2" + comp]
        rel = comp in relative
        if rel and truth == 0:
            raise SystemExit(f"s2{comp} = 0 here; drop it from --relative.")
        err = {n: (arrs[n][:, col] - truth) / (truth if rel else 1.0) for n in ns}
        hat, sig = rf"\hat{{\sigma}}^2_{{{sub}}}", rf"\sigma^2_{{{sub}}}"
        ylabel = (rf"Relative error $({hat} - {sig}) / {sig}$" if rel
                  else rf"Signed error $({hat} - {sig})$")
        title = f"({'abcd'[k]}) {NAMES[comp]}, ${sig} = {truth:g}$"
        _panel(ax, err, ns, ylabel, title)
        ax.set_xticks(range(len(ns)))
        ax.set_xticklabels(labels, fontsize=8.5)
        print(f"  {comp:<4} ({'relative' if rel else 'signed'}):")
        for n in ns:
            print(f"    n={n:<6d} mean={err[n].mean():+.4f} SD={err[n].std(ddof=1):.4f}")
    axes[-1].set_xlabel("Sample size (n), SNPs (m), genes (G)", fontsize=10, labelpad=8)
    fig.suptitle(top, fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.98))

    os.makedirs(out_dir, exist_ok=True)
    suffix = "" if gp is None else f"_gp{gp}"
    path = os.path.join(out_dir, f"components_error_nmG_scaling{suffix}.pdf")
    fig.savefig(path, bbox_inches="tight", facecolor="white", format="pdf", pad_inches=0.1)
    print(f"Saved: {path}")


if __name__ == "__main__":
    main()
