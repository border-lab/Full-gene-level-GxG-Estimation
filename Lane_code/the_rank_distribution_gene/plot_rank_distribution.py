# -*- coding: utf-8 -*-
"""Histogram of the per-gene rank r (one number per line)."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

path = sys.argv[1] if len(sys.argv) > 1 else "ContiguousSNP_n64000m32000_G320_thr0.99.txt"
r = np.loadtxt(path, dtype=int)
stem = path.rsplit(".txt", 1)[0]
q = {p: np.percentile(r, p) for p in (50, 90, 95, 99)}

ink, ink2, blue, grid = "#0b0b0b", "#52514e", "#2a78d6", "#e6e5e1"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": ink2, "axes.labelcolor": ink,
                     "xtick.color": ink2, "ytick.color": ink2})
fig, a1 = plt.subplots(figsize=(8, 4.5), facecolor="#fcfcfb")
for a in (a1,):
    a.set_facecolor("#fcfcfb"); a.grid(axis="y", color=grid, lw=0.8); a.set_axisbelow(True)
    for s in ("top", "right"): a.spines[s].set_visible(False)

# optional 3rd arg: bin width in ranks (default 1); widen it when few genes span a wide range
w = int(sys.argv[3]) if len(sys.argv) > 3 else 1
lo = r.min() - r.min() % w
bins = np.arange(lo, r.max() + w + 1, w) - 0.5
a1.hist(r, bins=bins, color=blue, edgecolor="#fcfcfb", linewidth=1)
a1.axvline(q[50], color=ink2, ls="--", lw=1, label=f"median {q[50]:.0f}")
a1.set_ylim(0, a1.get_ylim()[1] * 1.1)   # headroom so the legend clears the bars
a1.legend(loc="upper right", frameon=False, labelcolor=ink)
a1.yaxis.set_major_locator(MaxNLocator(integer=True))
a1.set_xlabel("rank R explaining 99% of the variance of $Z_g$"); a1.set_ylabel("number of genes")
a1.set_title("The distribution of r", loc="left", color=ink)
fig.tight_layout()
# optional 2nd arg: comma-separated output formats, e.g. "pdf"
exts = sys.argv[2].split(",") if len(sys.argv) > 2 else ["png", "pdf"]
for ext in exts:
    fig.savefig(f"{stem}.{ext}", dpi=160, facecolor=fig.get_facecolor())
print(f"n_genes={len(r)} min={r.min()} max={r.max()} mean={r.mean():.2f} sd={r.std(ddof=1):.2f} "
      + " ".join(f"p{p}={v:.1f}" for p, v in q.items()))
print("saved", ", ".join(f"{stem}.{ext}" for ext in exts))
