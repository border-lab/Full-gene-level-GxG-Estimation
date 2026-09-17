# -*- coding: utf-8 -*-
import sys
import math

if len(sys.argv) < 2:
    print("Usage: python3 script.py <filename>")
    sys.exit(1)

filename = sys.argv[1]

# Result rows are "(s2gxg,s2e)" -- two variance components
# (pairwise-epistasis-only model).
labels = ["gxg", "e"]
cols = [[] for _ in labels]

with open(filename, 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        line = line.strip('()')
        parts = line.split(',')
        for i in range(len(labels)):
            cols[i].append(float(parts[i]))

n = len(cols[0])


def _median(c):
    s = sorted(c)
    k = len(s)
    return s[k // 2] if k % 2 else 0.5 * (s[k // 2 - 1] + s[k // 2])


for i, label in enumerate(labels):
    c = cols[i]
    mean = sum(c) / n
    std = math.sqrt(sum((x - mean) ** 2 for x in c) / n)
    se = std / math.sqrt(n)
    ci_low = mean - 1.96 * se
    ci_high = mean + 1.96 * se
    print(f"n = {n}") if i == 0 else None
    print(f"Column {i + 1} ({label:>3}) - Mean: {mean:.6f}, Median: {_median(c):.6f}, "
          f"Std: {std:.6f}, 95% CI: [{ci_low:.6f}, {ci_high:.6f}]")
