# -*- coding: utf-8 -*-
"""Summarise the per-replicate realized variances of all four components.

    python3 summarize_realized_variance.py <pheno_dir> <out_file>

<pheno_dir> is a Phenotype/y_<TAG> directory holding one vell_rep*.txt per
replicate, each written by Simulate_Phenotype.py as labelled lines:

    vell_gxg <Var-hat(H gamma)>
    vell_a   <Var-hat(Z_a beta)>
    vell_d   <Var-hat(Z_d delta)>
    vell_e   <Var-hat(e)>

All four are random quantities -- quadratic forms in that replicate's effect
draw -- with

    E[vell_gxg] = c * s2gxg   (c != 1; see compute_c_pooled)
    E[vell_a]   = s2a         (c_a = 1 exactly, standardized design)
    E[vell_d]   = s2d         (c_d = 1 exactly, standardized design)
    E[vell_e]   = s2e

Rather than keeping the 200 draws, <out_file> records mean / std / se / n for
each component:

    <name>_mean <m>       estimate of E[.]
    <name>_std  <s>       draw-to-draw scatter (ddof=0)
    <name>_se   <s/sqrt(n)>
    <name>_n    <n>       number of replicates found

OPTIONAL -- the pipeline does NOT call this.  Every replicate's realized
variances are already carried into its estimate row by Simulate_MCREML.py
(columns 6-9 of result/<FILENAME>.txt), which keeps them paired with the fit
they belong to; this script is for when the reduction alone is wanted, and it
must be run BEFORE the combine job deletes the phenotype directory.
"""
import glob
import math
import os
import sys

NAMES = ("vell_gxg", "vell_a", "vell_d", "vell_e")

if len(sys.argv) != 3:
    print(__doc__)
    sys.exit(1)

pheno_dir, out_file = sys.argv[1], sys.argv[2]

vals = {k: [] for k in NAMES}
for path in sorted(glob.glob(os.path.join(pheno_dir, "vell_rep*.txt"))):
    with open(path) as f:
        rec = {k: float(v) for k, v in
               (line.split() for line in f if line.strip())}
    for k in NAMES:
        if k in rec:
            vals[k].append(rec[k])

if not vals["vell_gxg"]:
    sys.exit(f"no vell_rep*.txt files with a value in {pheno_dir}")

os.makedirs(os.path.dirname(out_file) or ".", exist_ok=True)
with open(out_file, "w") as f:
    for k in NAMES:
        v = vals[k]
        if not v:
            continue
        n = len(v)
        mean = sum(v) / n
        std = math.sqrt(sum((x - mean) ** 2 for x in v) / n)
        f.write(f"{k}_mean {mean:.10f}\n")
        f.write(f"{k}_std {std:.10f}\n")
        f.write(f"{k}_se {std / math.sqrt(n):.10f}\n")
        f.write(f"{k}_n {n}\n")
        print(f"{k}: mean={mean:.6f} std={std:.6f} (n={n})")
print(f"-> {out_file}")
