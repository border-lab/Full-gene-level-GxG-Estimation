# -*- coding: utf-8 -*-
"""
se_scaling_plot.py -- empirical SE of each variance-component estimate, one panel
per component stacked in a column, for a folder where n, m and G are scaled
together.  Companion to components_error_scaling_panels.py (same layout, same
x labels, same --panels / --relative / --gp options).

Usage
-----
    python se_scaling_plot.py <result_dir>
        [--panels a,d,gxg] [--relative a,d,gxg] [--gp PCT] [--outdir DIR]

Each panel shows, per setting, the SD (ddof=1) over replicates of the error of
that component -- exactly the "SD=" number printed above each box in the
components error figure -- with a 95% chi-square confidence interval.  The
dashed grey line is the 1/sqrt(n) rate anchored at the first setting.
Writes SE_nmG_scaling[_gp<PCT>].pdf into the result folder (or --outdir).
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
from gxg_col3_error_plot import _read
from gxg_col3_error_scaling_plot import _discover_triples, COMPONENTS
from components_error_scaling_panels import BOX_COLOR, NAMES, _opt


def _sd_ci(v, level=0.95):
    """Sample SD and its chi-square confidence interval."""
    k = len(v) - 1
    s = v.std(ddof=1)
    a = (1 - level) / 2
    return s, s * np.sqrt(k / stats.chi2.ppf(1 - a, k)), s * np.sqrt(k / stats.chi2.ppf(a, k))


def _panel(ax, err, ns, ylabel, title):
    sd, lo, hi = map(np.array, zip(*(_sd_ci(err[n]) for n in ns)))
    x = np.arange(len(ns))
    ymax = hi.max()
    ytop = ymax * 1.25
    ax.set_ylim(0, ytop)
    ax.set_xlim(-0.6, len(ns) - 0.4)
    for i in range(len(ns)):
        if i % 2:
            ax.axvspan(i - 0.5, i + 0.5, facecolor="#E8E8E8", alpha=0.8, zorder=0)
    ref = sd[0] * np.sqrt(ns[0] / np.asarray(ns, float))
    ax.plot(x, ref, color="#666666", linestyle="--", linewidth=0.8, zorder=1,
            label=r"$\propto 1/\sqrt{n}$")
    ax.errorbar(x, sd, yerr=[sd - lo, hi - sd], color=BOX_COLOR, linewidth=1.5,
                marker="o", markersize=6, markerfacecolor="white",
                markeredgewidth=1.5, capsize=4, elinewidth=1.2, zorder=3,
                label="Empirical SE (95% CI)")
    for i, s in enumerate(sd):
        ax.text(i, ytop - 0.04 * ytop, f"SE={s:.3f}", ha="center", va="top",
                fontsize=8.5, fontweight="bold")
    step = _nice_tick_step(ymax)
    ax.set_yticks([round(k * step, 10) for k in range(int(np.floor(ymax / step)) + 1)])
    ax.set_ylabel(ylabel, fontsize=10, labelpad=6)
    ax.set_title(title, fontsize=10, loc="left", pad=6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="both", which="major", labelsize=9, colors="#333333")
    return sd


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
    fig, axes = plt.subplots(len(panels), 1,
                             figsize=(2.6 + len(ns) * 1.0, 3.4 * len(panels)))
    axes = np.atleast_1d(axes)
    print(f"Directory : {dir_path}")
    print("Replicates: " + ", ".join(f"n={n}:{arrs[n].shape[0]}" for n in ns))
    for k, (comp, ax) in enumerate(zip(panels, axes)):
        col, sub = COMPONENTS[comp]
        truth = truths["s2" + comp]
        rel = comp in relative
        if rel and truth == 0:
            raise SystemExit(f"s2{comp} = 0 here; drop it from --relative.")
        err = {n: (arrs[n][:, col] - truth) / (truth if rel else 1.0) for n in ns}
        hat, sig = rf"\hat{{\sigma}}^2_{{{sub}}}", rf"\sigma^2_{{{sub}}}"
        ylabel = (rf"SE of relative error $({hat} - {sig}) / {sig}$" if rel
                  else rf"SE of ${hat}$")
        title = f"({'abcd'[k]}) {NAMES[comp]}, ${sig} = {truth:g}$"
        sd = _panel(ax, err, ns, ylabel, title)
        ax.set_xticks(range(len(ns)))
        ax.set_xticklabels(labels, fontsize=8.5)
        print(f"  {comp:<4} ({'relative' if rel else 'signed'}): "
              + ", ".join(f"n={n}:{s:.4f}" for n, s in zip(ns, sd)))
    axes[0].legend(loc="lower left", frameon=False, fontsize=8.5)
    axes[-1].set_xlabel("Sample size (n), SNPs (m), genes (G)", fontsize=10, labelpad=8)
    fig.suptitle(top, fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.98))

    os.makedirs(out_dir, exist_ok=True)
    suffix = "" if gp is None else f"_gp{gp}"
    path = os.path.join(out_dir, f"SE_nmG_scaling{suffix}.pdf")
    fig.savefig(path, bbox_inches="tight", facecolor="white", format="pdf", pad_inches=0.1)
    print(f"Saved: {path}")


if __name__ == "__main__":
    main()
