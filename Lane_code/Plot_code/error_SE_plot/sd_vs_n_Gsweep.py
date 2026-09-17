# -*- coding: utf-8 -*-
"""
sd_vs_n_Gsweep.py -- SD of the s2gxg relative error vs sample size, one line per
number of pooled genes G (G = 1, 5, 10, 20).

x-axis : sample size n
y-axis : SD of the relative error of sigma^2_gxg  ( std(s2gxg_hat - truth), ddof=1
         -- shift-invariant, so identical to std of the estimates, and matches the
         "SD=" annotation on the fixed_m box plots ).
legend : G = 1 / 5 / 10 / 20.

All four conditions are the SAME dataset (chr1_10ksnp, m = 10000, 300 reps); only
the number of pooled genes differs.  Note the G1 / G5 / G20 result files live in
folders named "RandomSNP_..." but the data inside is chr1_10ksnp, and the G10 data
lives in the chr1_10ksnp_pooled_preW folder -- so the source folder is given
explicitly per G.

    python sd_vs_n_Gsweep.py
writes SE_G1_G5_G10_G20.pdf
"""
import os
import re
import glob

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")            # headless: save PDF without a display
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
TRUTH_S2GXG = 0.2
N_VALUES = [1000, 2000, 4000, 8000, 16000]

# G -> (source folder, colour).  Data type / m are auto-detected from the files.
SOURCES = {
    1:  ("RandomSNP_pooled_preW_G1_s2gxg0.2_s2e0.8",   "#C03D3E"),
    5:  ("RandomSNP_pooled_preW_G5_s2gxg0.2_s2e0.8",   "#3274A1"),
    10: ("chr1_10ksnp_pooled_preW_s2gxg0.2_s2e0.8",    "#E1812C"),
    20: ("RandomSNP_pooled_preW_G20_s2gxg0.2_s2e0.8",  "#3A923A"),
}


def file_basename_and_m(folder_path, gtoken):
    """Real result-file basename (matching _G<g>_) and its fixed m in folder."""
    for f in glob.glob(os.path.join(folder_path, f"*_{gtoken}_n*m*.txt")):
        mobj = re.search(r"^(.*)_n(\d+)m(\d+)\.txt$", os.path.basename(f))
        if mobj:
            return mobj.group(1), int(mobj.group(3))
    raise FileNotFoundError(f"No *_{gtoken}_n*m*.txt files in {folder_path}")


def sd_series(folder, g):
    """SD (ddof=1) of the s2gxg relative error at each n for one G."""
    folder_path = os.path.join(HERE, folder)
    basename, m_val = file_basename_and_m(folder_path, f"G{g}")
    sds = []
    for n in N_VALUES:
        fp = os.path.join(folder_path, f"{basename}_n{n}m{m_val}.txt")
        df = pd.read_csv(fp, header=None)
        s2gxg_hat = df[0].astype(str).str.replace("(", "", regex=False).astype(float)
        rel_err = s2gxg_hat - TRUTH_S2GXG          # = estimate - truth
        sds.append(rel_err.std(ddof=1))            # sample SD, matches box plots
    return np.array(sds), m_val


def main():
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 11,
        "axes.linewidth": 1,
        "figure.dpi": 150,
    })
    fig, ax = plt.subplots(figsize=(7, 5))

    m_seen = None
    print("SD of s2gxg relative error (ddof=1):")
    for g, (folder, color) in SOURCES.items():
        sds, m_val = sd_series(folder, g)
        m_seen = m_val
        ax.plot(N_VALUES, sds, marker="o", markersize=6, linewidth=1.8,
                color=color, label=f"G = {g}")
        print(f"  G={g:<2d} | " + ", ".join(f"n={n}:{s:.4f}" for n, s in zip(N_VALUES, sds)))

    ax.set_xscale("log")                            # geometric n -> even spacing
    ax.set_xticks(N_VALUES)
    ax.get_xaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
    ax.minorticks_off()

    ax.set_xlabel("Sample size (n)", fontsize=12)
    ax.set_ylabel(r"SD of relative error ($\sigma^2_{g\times g}$)", fontsize=12)
    ax.set_title(f"Estimation SD vs sample size  (pooled preW, m = {m_seen:,})",
                 fontsize=12, fontweight="bold")
    ax.set_ylim(bottom=0)
    ax.legend(title="Pooled genes", frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", linestyle=":", linewidth=0.5, alpha=0.5)

    plt.tight_layout()
    out = os.path.join(HERE, "SE_G1_G5_G10_G20.pdf")
    plt.savefig(out, bbox_inches="tight", format="pdf")
    plt.rcParams.update(plt.rcParamsDefault)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
