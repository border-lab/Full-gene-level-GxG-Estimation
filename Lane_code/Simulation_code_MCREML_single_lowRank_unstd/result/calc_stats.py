# -*- coding: utf-8 -*-
"""Summarise a row-per-replicate result file, whatever its column count.

    python3 calc_stats.py <filename> [more files...]

Result rows are "(V_gamma,V_e,V_l,realized_variance)": the two variance
components of the pooled within-gene pairwise-epistasis model, Vl_hat =
c_hat * s2gxg_hat (the c-corrected ESTIMATE of the realized variance
Var(H gamma), realized_variance.pdf), and that replicate's REALIZED
Var-hat(H gamma) itself -- so columns 3 and 4 are paired per replicate and
comparing their means is the check the c-correction is meant to pass.

Older files still parse: three-column "(s2gxg,s2e,Vl_hat)" and two-column
"(s2gxg,s2e)" runs predate the fourth column, a one-column file of raw values
is summarised under the label "Vl_realized", and a file with MORE columns than
the names below parses too -- the extra ones are reported as col5, col6, ...
So a future estimator that writes another column needs no edit here; only its
name would be added to ALL_LABELS.
"""
import sys
import math

ALL_LABELS = ["gxg", "e", "Vl_hat", "Vl_realized"]


def read_rows(filename):
    """Parse "(a,b,c)" rows into a list of float lists.

    Blank lines and #-comments are skipped; every row must have the same
    number of columns, since a short or long row means a truncated or
    mis-merged replicate rather than something worth averaging over.
    """
    rows = []
    with open(filename, 'r') as f:
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


def summarise(filename, show_name=False):
    rows = read_rows(filename)
    ncol = len(rows[0])
    n = len(rows)
    labels = labels_for(ncol)
    width = max(len(x) for x in labels)

    if show_name:
        print(f"== {filename}")
    print(f"n = {n}")
    for i, label in enumerate(labels):
        c = [r[i] for r in rows]
        mean = sum(c) / n
        std = math.sqrt(sum((x - mean) ** 2 for x in c) / n)
        se = std / math.sqrt(n)
        print(f"Column {i + 1} ({label:>{width}}) - Mean: {mean:.6f}, "
              f"Median: {_median(c):.6f}, Std: {std:.6f}, "
              f"95% CI: [{mean - 1.96 * se:.6f}, {mean + 1.96 * se:.6f}]")


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 calc_stats.py <filename> [more files...]")
        sys.exit(1)

    files = sys.argv[1:]
    for k, filename in enumerate(files):
        if k:
            print()
        summarise(filename, show_name=len(files) > 1)


if __name__ == "__main__":
    main()
