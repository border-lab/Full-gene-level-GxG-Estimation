# -*- coding: utf-8 -*-
"""Summarise a row-per-replicate result file.

    python3 calc_stats.py [--drop-boundary] <filename> [more files...]

Result rows written by this pipeline are 2 columns,

    (s2a, s2e)
     a    e

the two variance components of the additive-only model, and nothing else.
The additive design is column-standardized, so both fitted components are on
the realized-variance scale and directly comparable to the nominal targets --
with s2a = 0.1, s2e = 0.9 the two column means should sit at 0.1 / 0.9.  No
post-fit correction is applied anywhere.

The column count alone cannot tell a sibling pipeline's layout apart from
this one's, so read a sibling's result file with that pipeline's own
calc_stats.py; this script would silently mislabel it.

BOUNDARY REPLICATES.  mc_reml constrains s2a >= 0 and s2e >= 1e-9 var(y)
(BOLT-REML's parameter domain), so a replicate whose likelihood peaks at
s2a = 0 ends ON that bound.  A component <= 1e-8 counts as at its bound, and
the fraction of replicates with each component there is reported.
--drop-boundary recomputes the moments without those rows, but they are NOT
dropped by default: at these sample sizes a zero estimate is a real outcome,
and discarding it would bias the means up.  Whether a replicate converged is
in the MC-AI-REML job's stdout log, not in the result file.
"""
import sys
import math

COMPONENTS = ["a", "e"]

# A component at or below this counts as at its bound (s2a sits at exactly 0
# there).  A decade above s2e's 1e-9 var(y) floor, so a small converged value
# is not mistaken for a bound hit.
CLAMP_LO = 1e-8


def read_rows(filename):
    """Parse "(a,b)" rows into a list of float lists.

    Blank lines and #-comments are skipped; every row must have the same
    number of columns, since a short or long row means a truncated or
    mis-merged replicate rather than something worth averaging over.
    """
    rows = []
    # utf-8-sig, not the locale default: a file that has been through a
    # Windows editor carries a BOM, and a machine whose default encoding is not
    # UTF-8 (gbk, say) then dies on byte 0 with an unhelpful UnicodeDecodeError.
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
    names = COMPONENTS[:ncol]
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


def bound_flags(rows, labels):
    """Per replicate, the components at their lower bound: {name: [bool]}."""
    out = {}
    for name in COMPONENTS:
        if name in labels:
            i = labels.index(name)
            out[name] = [row[i] <= CLAMP_LO for row in rows]
    return out


def summarise(filename, show_name=False, drop_boundary=False):
    rows = read_rows(filename)
    labels = labels_for(len(rows[0]))
    width = max(len(x) for x in labels)
    n_all = len(rows)

    flags = bound_flags(rows, labels)
    hits = [j for j in range(n_all) if any(f[j] for f in flags.values())]

    # --- components at their bound, always over ALL replicates -------------
    detail = ", ".join(f"{name} {100.0 * sum(f) / n_all:.1f}% ({sum(f)})"
                       for name, f in flags.items())
    bound_line = (f"At lower bound (value <= {CLAMP_LO:g}): {detail}; "
                  f"any component {len(hits)}/{n_all}")

    if drop_boundary and hits:
        keep = set(range(n_all)) - set(hits)
        rows = [rows[j] for j in sorted(keep)]
        if not rows:
            raise SystemExit(f"{filename}: every replicate has a component at "
                             f"its bound; nothing left to summarise.")
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

    print(bound_line)


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
