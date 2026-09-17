# -*- coding: utf-8 -*-
"""
plot_time.py -- estimation-time comparison: matrix-free vs pre-computed W.

x-axis : individual size n
y-axis : average estimation wall-clock time (s)
legend : one line per CONDITION (matrix-free / pre-computed W), auto-discovered
         from the sub-folders next to this script.

Each condition lives in a folder named like
    RandomSNP_<cond>_s2gxg<..>_s2e<..>/
holding one-number files
    timing_<foldername>_n<N>m<M>.txt
(the single average produced by the pipeline's time step).  n (and m) are read
from the file names, so adding more n values -- or another condition folder --
needs no code change.

The y-axis is on a normal (linear) scale.  Matrix-free is ~10^2-10^3x slower
than pre-computed W, so the pre-computed-W line sits close to the axis; a
power-law trend (time = a * n^b) is fitted per condition and the exponent b
(shown in the legend) still recovers the complexity in n: ~1 for matrix-free
O(n), ~2 for the dense O(n^2) W @ b mat-vec.

    python plot_time.py
writes time_comparison.pdf next to this script.
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

# condition token (in the folder name) -> pretty legend label + colour
_COND_STYLE = {
    "matfree": ("Matrix-free", "#C03D3E"),          # slow: rebuilds W every CG pass
    "preW":    ("Pre-computed W", "#3274A1"),        # fast: cached dense W @ b
}
_FALLBACK_COLORS = ["#3A923A", "#9372B2", "#E07850", "#8C564B"]


def discover_conditions(root):
    """Find condition folders and read their (n -> time) series.

    Returns a list of dicts: {token, label, color, n, t} sorted by folder name.
    """
    conditions = []
    fallback_i = 0
    for path in sorted(glob.glob(os.path.join(root, "*_s2gxg*_s2e*"))):
        if not os.path.isdir(path):
            continue
        folder = os.path.basename(path)
        # condition token = the piece between the mode and "_s2gxg"
        mobj = re.match(r"[^_]+_(.+?)_s2gxg", folder)
        token = mobj.group(1) if mobj else folder

        # gather (n, time) from every timing_..._n<N>m<M>.txt in the folder
        n_to_t = {}
        for f in glob.glob(os.path.join(path, "timing_*_n*m*.txt")):
            fm = re.search(r"_n(\d+)m(\d+)\.txt$", os.path.basename(f))
            if not fm:
                continue
            n_val = int(fm.group(1))
            with open(f) as fh:
                s = fh.read().strip()
            if s:
                n_to_t[n_val] = float(s)
        if not n_to_t:
            continue

        if token in _COND_STYLE:
            label, color = _COND_STYLE[token]
        else:
            label, color = token, _FALLBACK_COLORS[fallback_i % len(_FALLBACK_COLORS)]
            fallback_i += 1

        ns = np.array(sorted(n_to_t))
        ts = np.array([n_to_t[n] for n in ns], dtype=float)
        conditions.append(dict(token=token, label=label, color=color, n=ns, t=ts))
    return conditions


def _power_law(x, a, b):
    return a * np.power(x, b)


def fit_power_law(n, t):
    """Fit time = a * n^b via a log-log linear least squares (robust seed)."""
    b0, loga0 = np.polyfit(np.log(n), np.log(t), 1)
    try:
        (a, b), _ = curve_fit(_power_law, n, t, p0=(np.exp(loga0), b0), maxfev=10000)
    except Exception:
        a, b = np.exp(loga0), b0
    return a, b


def main():
    conditions = discover_conditions(HERE)
    if not conditions:
        raise SystemExit(f"No condition folders (*_s2gxg*_s2e*) with timing files found in {HERE}")

    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 10,
        "axes.linewidth": 1,
        "figure.dpi": 150,
    })

    fig, ax = plt.subplots(figsize=(7, 5))

    for cond in conditions:
        n, t = cond["n"], cond["t"]
        ax.scatter(n, t, marker="o", s=34, color=cond["color"], zorder=3)

        a, b = fit_power_law(n, t)
        n_fit = np.linspace(n.min(), n.max(), 200)
        ax.plot(n_fit, _power_law(n_fit, a, b), "--", color=cond["color"],
                linewidth=1.6, label=f"{cond['label']}  (~n^{b:.2f})")

    ax.set_xlabel("Individual Size (n)", fontsize=11)
    ax.set_ylabel("Estimation Time (s)", fontsize=11)
    ax.set_title("Estimation time: matrix-free vs pre-computed W",
                 fontsize=12, fontweight="bold")
    ax.set_ylim(bottom=0)
    ax.legend(title="Condition", frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", linestyle=":", linewidth=0.5, alpha=0.5)

    plt.tight_layout()
    pdf = os.path.join(HERE, "time_comparison.pdf")
    plt.savefig(pdf, bbox_inches="tight")
    plt.rcParams.update(plt.rcParamsDefault)

    # console summary
    print("Conditions plotted:")
    for cond in conditions:
        a, b = fit_power_law(cond["n"], cond["t"])
        pts = ", ".join(f"n={n}:{t:.2f}s" for n, t in zip(cond["n"], cond["t"]))
        print(f"  {cond['label']:16s} exponent~{b:.2f} | {pts}")
    print(f"\nSaved: {pdf}")


if __name__ == "__main__":
    main()
