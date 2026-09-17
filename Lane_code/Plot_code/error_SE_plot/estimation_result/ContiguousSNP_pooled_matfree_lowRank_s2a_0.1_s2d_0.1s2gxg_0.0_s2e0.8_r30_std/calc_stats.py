# -*- coding: utf-8 -*-
"""Summarise a row-per-replicate result file, whatever its column count.

    python3 calc_stats.py [--drop-boundary] <filename> [more files...]

Result rows written by this pipeline are 9 columns,

    (s2a, s2d, s2gxg, s2e, Vl_hat, vell_gxg, vell_a, vell_d, vell_e)
     a    d    gxg    e    Vl_hat  Vl_real   Va_real  Vd_real  Ve_real

the four variance components of the additive + dominance + pooled within-gene
pairwise-epistasis model, then Vl_hat, then the realized variances Var-hat(.)
that this replicate's four effect draws actually produced.

THIS PIPELINE'S EPISTASIS KERNEL IS C-NORMALIZED (W = W_raw / c), so ALL FOUR
fitted components are already on the realized-variance scale and columns 1-4
are directly comparable to the nominal targets -- with s2a = s2d = s2gxg = 0.1,
s2e = 0.7 the four column means should sit at 0.1 / 0.1 / 0.1 / 0.7.  No
post-fit correction is applied anywhere, and column 5 (Vl_hat) is IDENTICALLY
column 3: it is kept only so the row layout matches the sibling pipelines'.
In the _unstd_4VC sibling column 3 lived on the raw-H component scale and only
column 5 = c_hat * s2gxg_hat was comparable to a variance -- do not carry that
reading over to a file from this pipeline.

Narrower files still parse: a 5-column file is the same layout without the
realized columns, and a file with MORE columns than the names below parses too
-- the extras are reported as col10, col11, ...  The labels assume the ADDITIVE
and DOMINANCE columns come first, so a two- or three-component sibling's result
file should be read with that pipeline's own calc_stats.py, not this one: the
column count alone cannot tell those layouts apart, and this script would
silently mislabel them.

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
