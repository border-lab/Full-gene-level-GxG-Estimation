# -*- coding: utf-8 -*-
"""
plot_time_buildW.py -- time comparison for three explicit conditions:

    1. RandomSNP_matfree_s2gxg0.3_s2e0.7          -> matrix-free estimation time
    2. RandomSNP_preW_s2gxg0.3_s2e0.7             -> pre-computed-W estimation time
    3. RandomSNP_preW_s2gxg0.3_s2e0.7_buildW      -> pre-computed-W kernel-build time

x-axis : individual size n
y-axis : wall-clock time (s)
legend : one line per condition.

Each folder holds one-number timing files named  <...>_n<N>m<M>.txt  (the single
average written by the pipeline's time step).  Files that instead contain
"(s2gxg,s2e)" result tuples -- e.g. a folder accidentally filled with estimation
RESULTS rather than TIMES -- are detected and skipped with a warning, so the plot
never mistakes a result file for a timing.

    python plot_time_buildW.py
writes time_comparison_buildW.pdf next to this script.
"""
import os
import re
import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")            # headless: save without a display
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))

# (folder name, pretty label, colour) -- explicit, in legend order.
_SERIES = [
    ("RandomSNP_matfree_s2gxg0.3_s2e0.7",       "Matrix-free (estimation)",     "#C03D3E"),
    ("RandomSNP_preW_s2gxg0.3_s2e0.7",          "Pre-computed W (estimation)",  "#3274A1"),
    ("RandomSNP_preW_s2gxg0.3_s2e0.7_buildW",   "Pre-computed W (build W)",     "#3A923A"),
]


def _parse_single_time(text):
    """Return the float in a one-number timing file, or None if the file is not
    a timing file (empty, or holds "(s2gxg,s2e)" result tuples)."""
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        # A result tuple like "(0.18,0.73)" is NOT a timing -- reject it.
        if "(" in line or "," in line:
            return None
        try:
            return float(line)
        except ValueError:
            return None
    return None


def read_series(folder):
    """Read {n -> time} from one folder's  *_n<N>m<M>.txt  timing files.

    Returns (ns, ts) sorted by n, or (None, None) if the folder has no valid
    single-number timing files (e.g. it holds result tuples instead)."""
    path = os.path.join(HERE, folder)
    n_to_t = {}
    saw_nontiming = False
    for f in glob.glob(os.path.join(path, "*_n*m*.txt")):
        fm = re.search(r"_n(\d+)m(\d+)\.txt$", os.path.basename(f))
        if not fm:
            continue
        with open(f) as fh:
            t = _parse_single_time(fh.read())
        if t is None:
            saw_nontiming = True
            continue
        n_to_t[int(fm.group(1))] = t

    if not n_to_t:
        if saw_nontiming:
            print(f"  [skip] {folder}: files are not timing numbers "
                  f"(look like estimation-result tuples).")
        else:
            print(f"  [skip] {folder}: no  *_n<N>m<M>.txt  files found.")
        return None, None

    ns = np.array(sorted(n_to_t))
    ts = np.array([n_to_t[n] for n in ns], dtype=float)
    return ns, ts


def _power_law(x, a, b):
    return a * np.power(x, b)


def fit_power_law(n, t):
    """Fit time = a * n^b via log-log seed + curve_fit (needs >= 2 points)."""
    if len(n) < 2:
        return None, None
    b0, loga0 = np.polyfit(np.log(n), np.log(t), 1)
    try:
        (a, b), _ = curve_fit(_power_law, n, t, p0=(np.exp(loga0), b0), maxfev=10000)
    except Exception:
        a, b = np.exp(loga0), b0
    return a, b


def main():
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 10,
        "axes.linewidth": 1,
        "figure.dpi": 150,
    })
    fig, ax = plt.subplots(figsize=(7, 5))

    plotted = []
    print("Reading series:")
    for folder, label, color in _SERIES:
        n, t = read_series(folder)
        if n is None:
            continue
        ax.scatter(n, t, marker="o", s=34, color=color, zorder=3)
        a, b = fit_power_law(n, t)
        if b is not None:
            n_fit = np.linspace(n.min(), n.max(), 200)
            ax.plot(n_fit, _power_law(n_fit, a, b), "--", color=color,
                    linewidth=1.6, label=f"{label}  (~n^{b:.2f})")
        else:
            ax.plot(n, t, "--", color=color, linewidth=1.6, label=label)
        plotted.append((label, n, t, b))

    if not plotted:
        raise SystemExit("No valid timing series found -- nothing to plot.")

    ax.set_xlabel("Individual Size (n)", fontsize=11)
    ax.set_ylabel("Time (s)", fontsize=11)
    ax.set_title("Time: matrix-free vs pre-computed W (estimation + build)",
                 fontsize=12, fontweight="bold")
    ax.set_ylim(bottom=0)
    ax.legend(title="Condition", frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", linestyle=":", linewidth=0.5, alpha=0.5)

    plt.tight_layout()
    pdf = os.path.join(HERE, "time_comparison_buildW.pdf")
    plt.savefig(pdf, bbox_inches="tight")
    plt.rcParams.update(plt.rcParamsDefault)

    print("\nPlotted:")
    for label, n, t, b in plotted:
        pts = ", ".join(f"n={nn}:{tt:.2f}s" for nn, tt in zip(n, t))
        exp = f"exponent~{b:.2f}" if b is not None else "single point"
        print(f"  {label:30s} {exp} | {pts}")
    print(f"\nSaved: {pdf}")


if __name__ == "__main__":
    main()
