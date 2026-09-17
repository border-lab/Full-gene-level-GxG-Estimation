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
draw -- and in THIS pipeline all four are centred on their NOMINAL targets,
because the epistasis kernel is c-normalized (W = W_raw / c-hat):

    E[vell_gxg] = (c / c-hat) s2gxg
                              (c-hat divided out of the kernel; = c * s2gxg in
                               the _unstd_4VC sibling, where c != 1)
    E[vell_a]   = s2a         (c_a = 1 exactly, standardized design)
    E[vell_d]   = s2d         (c_d = 1 exactly, standardized design)
    E[vell_e]   = s2e

So a run at s2a = s2d = s2gxg = 0.1, s2e = 0.7 should show the four means at
0.1 / 0.1 / 0.1 / 0.7, each within a couple of its own reported se.  That is
the cheapest end-to-end check that the normalization is in force.

The epistasis mean is the one to read with its factor in hand: c-hat is the
O(nm) third-moment plug-in (Function_MCREML.C_METHOD = 'moment'), so vell_gxg
centres on (c_exact / c-hat) * s2gxg, and that ratio is written out per run as
c_gxg_after_normalization in result/c_<FILENAME>.txt.  A vell_gxg_mean that
sits several se from s2gxg but ON that number is the plug-in's residual, not a
broken normalization; one that sits off BOTH is worth chasing.

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
