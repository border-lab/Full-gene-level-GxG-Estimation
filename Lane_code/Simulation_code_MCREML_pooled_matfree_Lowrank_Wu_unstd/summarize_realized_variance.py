# -*- coding: utf-8 -*-
"""Summarise the per-replicate realized variances Var-hat(H gamma).

    python3 summarize_realized_variance.py <rep_dir> <out_file>

<rep_dir> holds one rep*.txt per replicate (a single number each, written by
Simulate_Phenotype.py).  The realized variance is a random quantity -- a
quadratic form in the gamma draw -- with E[V_ell] = c * s2gxg (see
realized_variance.pdf and compute_c_pooled).  Rather than keeping the 300
draws, <out_file> records their mean and standard deviation:

    mean <m>       estimate of E[V_ell] = c * s2gxg
    std  <s>       draw-to-draw scatter of V_ell (ddof=0, as in calc_stats.py)
    se   <s/sqrt(n)>
    n    <n>       number of replicates found

Called by combine_code.sh; the pipeline's cleanup then deletes <rep_dir>.
"""
import glob
import math
import os
import sys

if len(sys.argv) != 3:
    print(__doc__)
    sys.exit(1)

rep_dir, out_file = sys.argv[1], sys.argv[2]

vals = []
for path in sorted(glob.glob(os.path.join(rep_dir, "rep*.txt"))):
    with open(path) as f:
        txt = f.read().strip()
    if txt:
        vals.append(float(txt))

if not vals:
    sys.exit(f"no rep*.txt files with a value in {rep_dir}")

n = len(vals)
mean = sum(vals) / n
std = math.sqrt(sum((v - mean) ** 2 for v in vals) / n)

os.makedirs(os.path.dirname(out_file) or ".", exist_ok=True)
with open(out_file, "w") as f:
    f.write(f"mean {mean:.10f}\n")
    f.write(f"std {std:.10f}\n")
    f.write(f"se {std / math.sqrt(n):.10f}\n")
    f.write(f"n {n}\n")
print(f"realized variance: mean={mean:.6f} std={std:.6f} (n={n}) -> {out_file}")
