# -*- coding: utf-8 -*-
import sys
import math

if len(sys.argv) < 2:
    print("Usage: python3 script.py <filename>")
    sys.exit(1)

filename = sys.argv[1]
col1 = []
col2 = []

with open(filename, 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        line = line.strip('()')
        parts = line.split(',')
        col1.append(float(parts[0]))
        col2.append(float(parts[1]))

n = len(col1)
mean1 = sum(col1) / n
mean2 = sum(col2) / n

std1 = math.sqrt(sum((x - mean1) ** 2 for x in col1) / n)
std2 = math.sqrt(sum((x - mean2) ** 2 for x in col2) / n)

se1 = std1 / math.sqrt(n)
se2 = std2 / math.sqrt(n)

# 95% CI = mean +/- 1.96 * SE
ci1_low = mean1 - 1.96 * se1
ci1_high = mean1 + 1.96 * se1
ci2_low = mean2 - 1.96 * se2
ci2_high = mean2 + 1.96 * se2

print(f"n = {n}")
print(f"Column 1 - Mean: {mean1:.6f}, Std: {std1:.6f}, 95% CI: [{ci1_low:.6f}, {ci1_high:.6f}]")
print(f"Column 2 - Mean: {mean2:.6f}, Std: {std2:.6f}, 95% CI: [{ci2_low:.6f}, {ci2_high:.6f}]")