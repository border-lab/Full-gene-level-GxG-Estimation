# -*- coding: utf-8 -*-
"""Summarise a row-per-replicate result file, whatever its column count.

    python3 calc_stats.py [--drop-boundary] <filename> [more files...]

Result rows written by this pipeline are 4 columns,

    (s2a, s2d, s2gxg, s2e)
     a    d    gxg    e

the four variance components of the additive + dominance + pooled within-gene
pairwise-epistasis model, and nothing else.

THIS PIPELINE'S EPISTASIS KERNEL IS C-NORMALIZED (W = W_raw / c-hat), so all
four fitted components are already on the realized-variance scale and all four
are directly comparable to the nominal targets -- with s2a = s2d = s2gxg = 0.1,
s2e = 0.7 the four column means should sit at 0.1 / 0.1 / 0.1 / 0.7.  No
post-fit correction is applied anywhere.

WHY THERE ARE NO REALIZED COLUMNS ANY MORE.  Rows used to carry five more:
Vl_hat (identically column 3 here) and the four realized variances Var-hat(.)
of that replicate's effect draws, so an estimate could be judged against its own
draw rather than the nominal target.  simulate_phenotype now rescales each drawn
component to hit its target exactly (force_realized=True), which makes those
four columns constants equal to the targets -- so the paired comparison and the
nominal one became the same comparison, and the columns went.  Measured over 30
replicates, the old paired spread and the new unpaired spread agree to within
Monte Carlo error.

It also removes the caveat that used to attach to column 3: c-hat is the O(nm)
third-moment plug-in (Function_MCREML.C_METHOD), not the exact c, so the
epistasis target used to be (c_exact / c-hat) * s2gxg rather than s2gxg.
Forcing the realized variance absorbs that factor at the draw, so column 3's
target is now the nominal number.  The ratio is still written out per run as
c_gxg_after_normalization in result/c_<FILENAME>.txt, as a property of the
genotype; no column depends on it.

WIDER FILES STILL PARSE, and that is the trap: an OLD 9-column file from this
same pipeline, or a sibling's, is labelled by the same first-four names and then
Vl_hat, Vl_real, ... -- the paired block below reappears and the summary is
correct for that file.  But the column count alone cannot tell a two- or
three-component sibling's layout apart from this one, so read a sibling's result
file with that pipeline's own calc_stats.py; this script would silently
mislabel it.

Two things beyond per-column moments, both of which exist because the columns
come in ESTIMATE / REALIZED pairs:

  PAIRED DIFFERENCES.  Each estimate has the draw-level quantity it should
  recover in the SAME row, so the mean of (estimate - realized) over replicates
  is the bias against what was actually drawn -- a sharper target than the
  nominal s2a / s2d / s2gxg / s2e, which the draws themselves scatter around by
  O(1/sqrt n).  Reported for every pair present.  The gxg pair is listed as
  Vl_hat - Vl_real; since Vl_hat = s2gxg_hat here, that IS the s2gxg pair, and
  it is not printed twice.

  BOUNDARY REPLICATES.  mc_reml clamps components to [1e-9, 5 var(y)], and a
  replicate whose likelihood peaks at a component = 0 sticks there -- after
  which the OTHER components converge to the wrong values, because the
  AI-Newton step is never re-projected onto the free subspace (the known defect
  documented in mc_reml).  Such rows are counted and named; --drop-boundary
  recomputes everything without them.  They are NOT dropped by default: they
  are a real property of the estimator at these sample sizes, and silently
  discarding them would flatter the summary.  Expect the DOMINANCE component to
  be the most frequent offender -- it carries the weakest signal of the four --
  so read the per-component counts, not just the total.

(realized_variance/<FILENAME>.txt is NOT a row file: it already holds the
mean/std of the realized variances over the reps, written by
summarize_realized_variance.py.)
"""
import sys
import math

ALL_LABELS = ["a", "d", "gxg", "e", "Vl_hat",
              "Vl_real", "Va_real", "Vd_real", "Ve_real"]

# (estimate label, realized label) -- the pairs the row layout puts side by side.
PAIRS = [("Vl_hat", "Vl_real"), ("a", "Va_real"), ("d", "Vd_real"),
         ("e", "Ve_real")]

# Components mc_reml can clamp to its lower bound (1e-9).  The threshold sits a
# decade above it, so a value that merely converged very small is not mistaken
# for a clamped one.
FITTED = ["a", "d", "gxg", "e"]
CLAMP_LO = 1e-8


def read_rows(filename):
    """Parse "(a,b,c)" rows into a list of float lists.

    Blank lines and #-comments are skipped; every row must have the same
    number of columns, since a short or long row means a truncated or
    mis-merged replicate rather than something worth averaging over.
    """
    rows = []
    # utf-8-sig, not the locale default: the rows are ASCII either way, but a
    # file that has been through a Windows editor carries a BOM, and a machine
    # whose default encoding is not UTF-8 (gbk, say) then dies on byte 0 with a
    # UnicodeDecodeError that says nothing about the real problem.
    with open(filename, 'r', encoding='utf-8-sig') as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            line = line.strip('()[] ')
            try:
                row = [float(v) for v in line.split(',')]
            except ValueError:
                raise SystemExit(f"{filename}:{lineno}: cannot parse '{line}' as numbers.")
            if rows and len(row) != len(rows[0]):
                raise SystemExit(
                    f"{filename}:{lineno}: {len(row)} columns, but the first row has "
                    f"{len(rows[0])}. Ragged file -- check for a truncated replicate.")
            rows.append(row)
    if not rows:
        raise SystemExit(f"{filename}: no data rows.")
    return rows


def labels_for(ncol):
    """Column names for a file of this width: the known ones, then col<i>."""
    if ncol == 1:
        return ["Vl_realized"]
    names = ALL_LABELS[:ncol]
    names += [f"col{i}" for i in range(len(names) + 1, ncol + 1)]
    return names


def _median(c):
    s = sorted(c)
    k = len(s)
    return s[k // 2] if k % 2 else 0.5 * (s[k // 2 - 1] + s[k // 2])


def _moments(c):
    """(mean, std, se) with std at ddof=0, as everywhere else in the pipeline."""
    n = len(c)
    mean = sum(c) / n
    std = math.sqrt(sum((x - mean) ** 2 for x in c) / n)
    return mean, std, std / math.sqrt(n)


def boundary_rows(rows, labels):
    """Indices of replicates with a fitted component at the lower clamp.

    Returns (indices, per-label counts).  Only the components mc_reml actually
    optimizes are checked: Vl_hat is a positive multiple of s2gxg_hat and would
    just re-report the same replicates, and the realized columns are draws, not
    fits.
    """
    idx = {}
    for name in FITTED:
        if name in labels:
            idx[name] = labels.index(name)
    counts = {name: 0 for name in idx}
    hits = []
    for j, row in enumerate(rows):
        clamped = [name for name, i in idx.items() if row[i] <= CLAMP_LO]
        if clamped:
            hits.append(j)
            for name in clamped:
                counts[name] += 1
    return hits, counts


def summarise(filename, show_name=False, drop_boundary=False):
    rows = read_rows(filename)
    labels = labels_for(len(rows[0]))
    width = max(len(x) for x in labels)

    hits, counts = boundary_rows(rows, labels)
    n_all = len(rows)
    if drop_boundary and hits:
        keep = set(range(n_all)) - set(hits)
        rows = [rows[j] for j in sorted(keep)]
        if not rows:
            raise SystemExit(f"{filename}: every replicate is boundary-clamped; "
                             f"nothing left to summarise.")
    n = len(rows)

    if show_name:
        print(f"== {filename}")
    print(f"n = {n}" + ("  (boundary replicates dropped)"
                        if drop_boundary and hits else ""))

    for i, label in enumerate(labels):
        c = [r[i] for r in rows]
        mean, std, se = _moments(c)
        print(f"Column {i + 1} ({label:>{width}}) - Mean: {mean:.6f}, "
              f"Median: {_median(c):.6f}, Std: {std:.6f}, "
              f"95% CI: [{mean - 1.96 * se:.6f}, {mean + 1.96 * se:.6f}]")

    # --- estimate vs. what this replicate actually realized -----------------
    pairs = [(e, r) for e, r in PAIRS if e in labels and r in labels]
    if pairs:
        print("Paired (estimate - realized), the bias against each replicate's own draw:")
        for est, real in pairs:
            ie, ir = labels.index(est), labels.index(real)
            d = [r[ie] - r[ir] for r in rows]
            mean, std, se = _moments(d)
            flag = "" if abs(mean) <= 1.96 * se else "   *"
            print(f"  {est:>{width}} - {real:<{width}} - Mean: {mean:+.6f}, "
                  f"Std: {std:.6f}, 95% CI: [{mean - 1.96 * se:+.6f}, "
                  f"{mean + 1.96 * se:+.6f}]{flag}")
        print("  (* = CI excludes 0)")

    # --- the optimizer's boundary cases -------------------------------------
    if hits:
        detail = ", ".join(f"{name}={counts[name]}" for name in FITTED
                           if counts.get(name))
        print(f"Boundary: {len(hits)}/{n_all} replicates have a component at the "
              f"lower clamp ({detail}).")
        if not drop_boundary:
            print("  Their OTHER components are biased too (see mc_reml's known "
                  "defect); rerun with --drop-boundary to exclude them.")


def main():
    args = [a for a in sys.argv[1:] if a != "--drop-boundary"]
    drop = "--drop-boundary" in sys.argv[1:]

    if not args:
        print("Usage: python3 calc_stats.py [--drop-boundary] <filename> [more files...]")
        sys.exit(1)

    for k, filename in enumerate(args):
        if k:
            print()
        summarise(filename, show_name=len(args) > 1, drop_boundary=drop)


if __name__ == "__main__":
    main()
