# -*- coding: utf-8 -*-
"""
null_experiment_plot.py -- does the epistasis estimator report nothing when
there is nothing?  Two panels from the raw s2gxg estimates (column 3):

  (a) NULL vs ALTERNATIVE.  One pair of boxes per sample size: the null run
      (true sigma^2_gxg = 0) beside a matched run with epistasis present, with
      dashed lines at both truths.  A working estimator puts the null boxes on
      0 and the alternative boxes on their truth, and separates the two more
      cleanly as n grows.  A null run alone cannot show the second half.

  (b) SPREAD vs n, log-log.  The SD of s2gxg_hat over replicates for both
      runs, with slope -1/2 and -1 guide lines through the first null point.
      A consistent estimator's SD falls with n; the slope says how fast.

Usage
-----
    python null_experiment_plot.py <null_dir> <alt_dir> [--outdir DIR]

Both folders hold row-per-replicate files "(s2a, s2d, s2gxg, s2e [, ...])";
only column 3 is read, so a 4-column and a 9-column layout can be compared.
The truths are parsed from the folder names, and only sample sizes present in
BOTH folders are plotted.  The PDF goes to <null_dir> unless --outdir is given.

What this does NOT show is calibration: whether a TEST of sigma^2_gxg = 0
rejects at its nominal rate.  That needs a per-replicate standard error, which
the result rows do not carry yet.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from realized_var_four_components_plot import _parse_truth, _discover, _tag
from gxg_col3_error_plot import _read, COL_GXG

NULL_COLOR = "#3274A1"
ALT_COLOR = "#E1812C"


def _load(dir_path):
    truth = _parse_truth(os.path.basename(dir_path))["s2gxg"]
    stem, ns, m_val, g_val = _discover(dir_path)
    est = {n: _read(os.path.join(dir_path, _tag(stem, n, m_val, g_val) + ".txt"))[:, COL_GXG]
           for n in ns}
    return truth, est, m_val, g_val


def _box(ax, data, pos, color):
    ax.boxplot(data, positions=[pos], widths=0.32, patch_artist=True,
               showfliers=False,
               boxprops=dict(linewidth=1.3, edgecolor=color, facecolor="white"),
               whiskerprops=dict(linewidth=1.1, color=color),
               capprops=dict(linewidth=1.1, color=color),
               medianprops=dict(linewidth=1.8, color=color))
    jitter = np.random.default_rng(0).uniform(-0.1, 0.1, data.size)
    ax.scatter(pos + jitter, data, s=5, color=color, alpha=0.35,
               linewidths=0, zorder=3)


def main():
    _usage = ("Usage: python null_experiment_plot.py <null_dir> <alt_dir> "
              "[--outdir DIR]")
    args = sys.argv[1:]
    out_dir = None
    if "--outdir" in args:
        i = args.index("--outdir")
        if i + 1 >= len(args):
            print("Error: --outdir needs a directory.")
            sys.exit(1)
        out_dir = os.path.abspath(args[i + 1])
        args = args[:i] + args[i + 2:]
    if len(args) != 2:
        print(_usage)
        sys.exit(1)

    null_dir, alt_dir = (os.path.abspath(a) for a in args)
    for d in (null_dir, alt_dir):
        if not os.path.isdir(d):
            print(f"Error: '{d}' is not a directory.")
            sys.exit(1)
    out_dir = out_dir or null_dir

    t0, est0, m0, g0 = _load(null_dir)
    t1, est1, m1, g1 = _load(alt_dir)
    if (m0, g0) != (m1, g1):
        raise SystemExit(f"m/G differ: null (m={m0}, G={g0}) vs alt (m={m1}, G={g1}).")
    if t0 != 0.0:
        print(f"Warning: the null folder's s2gxg is {t0}, not 0.")
    ns = sorted(set(est0) & set(est1))
    if not ns:
        raise SystemExit("The two folders share no sample size.")

    print(f"Null : {null_dir}  (s2gxg={t0})")
    print(f"Alt  : {alt_dir}  (s2gxg={t1})")
    print(f"m = {m0}" + ("" if g0 is None else f", G = {g0}"))
    sd0 = np.array([est0[n].std(ddof=1) for n in ns])
    sd1 = np.array([est1[n].std(ddof=1) for n in ns])
    for k, n in enumerate(ns):
        gap = (est1[n].mean() - est0[n].mean()) / np.hypot(sd0[k], sd1[k])
        print(f"  n={n:<6d} null mean={est0[n].mean():+.5f} SD={sd0[k]:.5f} | "
              f"alt mean={est1[n].mean():+.5f} SD={sd1[k]:.5f} | "
              f"separation (mean gap / pooled SD) = {gap:.2f}")
    slope0 = np.polyfit(np.log(ns), np.log(sd0), 1)[0]
    slope1 = np.polyfit(np.log(ns), np.log(sd1), 1)[0]
    print(f"  log-log SD slope: null {slope0:.2f}, alt {slope1:.2f}")

    plt.rcParams.update({"font.family": "Arial", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1.6, 1]})

    # --- (a) null vs alternative ------------------------------------------
    for k, n in enumerate(ns):
        if k % 2:
            ax_a.axvspan(k - 0.5, k + 0.5, facecolor="#E8E8E8", alpha=0.8, zorder=0)
        _box(ax_a, est0[n], k - 0.19, NULL_COLOR)
        _box(ax_a, est1[n], k + 0.19, ALT_COLOR)
    ax_a.axhline(t0, color=NULL_COLOR, linestyle="--", linewidth=0.9, zorder=1)
    ax_a.axhline(t1, color=ALT_COLOR, linestyle="--", linewidth=0.9, zorder=1)
    ax_a.set_xticks(range(len(ns)))
    ax_a.set_xticklabels([f"n = {n:,}" for n in ns], fontsize=9)
    ax_a.set_xlim(-0.6, len(ns) - 0.4)
    ax_a.set_xlabel("Sample size (n)")
    ax_a.set_ylabel(r"$\hat{\sigma}^2_{g \times g}$")
    ax_a.plot([], [], color=NULL_COLOR, lw=6, label=rf"null, $\sigma^2_{{g\times g}} = {t0:g}$")
    ax_a.plot([], [], color=ALT_COLOR, lw=6, label=rf"alternative, $\sigma^2_{{g\times g}} = {t1:g}$")
    ax_a.legend(frameon=False, loc="upper right", fontsize=9)
    ax_a.set_title("(a) Estimates under the null and the alternative",
                   fontsize=10, loc="left")

    # --- (b) SD vs n, log-log ----------------------------------------------
    ns_arr = np.array(ns, dtype=float)
    ax_b.plot(ns_arr, sd0, "o-", color=NULL_COLOR, lw=1.5, ms=5,
              label=f"null (slope {slope0:.2f})")
    ax_b.plot(ns_arr, sd1, "s-", color=ALT_COLOR, lw=1.5, ms=5,
              label=f"alternative (slope {slope1:.2f})")
    for p, ls in ((-0.5, ":"), (-1.0, "--")):
        ax_b.plot(ns_arr, sd0[0] * (ns_arr / ns_arr[0]) ** p, color="#888888",
                  linestyle=ls, lw=0.9, label=rf"$\propto n^{{{p:g}}}$")
    ax_b.set_xscale("log")
    ax_b.set_yscale("log")
    ax_b.set_xticks(ns_arr)
    ax_b.set_xticklabels([f"{n:,}" for n in ns], fontsize=9)
    ax_b.minorticks_off()
    # Log axes label only powers of ten, which leaves this narrow range with
    # one labelled tick; label a 1-2-5 ladder spanning the data instead.
    lo, hi = min(sd0.min(), sd1.min()), max(sd0.max(), sd1.max())
    yt = [b * 10.0 ** e for e in range(int(np.floor(np.log10(lo))),
                                       int(np.ceil(np.log10(hi))) + 1)
          for b in (1, 2, 5) if lo / 1.5 <= b * 10.0 ** e <= hi * 1.5]
    ax_b.set_yticks(yt)
    ax_b.set_yticklabels([f"{t:g}" for t in yt], fontsize=9)
    ax_b.set_xlabel("Sample size (n)")
    ax_b.set_ylabel(r"SD of $\hat{\sigma}^2_{g \times g}$ over replicates")
    ax_b.legend(frameon=False, fontsize=9)
    ax_b.set_title("(b) Spread shrinks with n", fontsize=10, loc="left")

    r0 = min(est0[n].size for n in ns)
    r1 = min(est1[n].size for n in ns)
    print(f"  replicates per box: null {r0}, alt {r1}")
    fig.suptitle(f"m = {m0:,}" + ("" if g0 is None else f",  G = {g0}")
                 + f",  R = {r0} (null) / {r1} (alternative) replicates",
                 fontsize=11, fontweight="bold")
    fig.tight_layout()

    gsuf = "" if g0 is None else f"_G{g0}"
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"null_experiment_m{m0}{gsuf}.pdf")
    fig.savefig(path, bbox_inches="tight", facecolor="white", format="pdf",
                pad_inches=0.1)
    print(f"Saved: {path}")


if __name__ == "__main__":
    main()
