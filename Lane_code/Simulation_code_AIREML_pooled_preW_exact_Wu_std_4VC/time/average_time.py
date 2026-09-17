# -*- coding: utf-8 -*-
"""Average the per-rep estimation times and write a single number.

    python3 average_time.py <rep_times_dir> <out_file>

<rep_times_dir> holds one-line rep*.txt files, each a wall-clock estimation
time in seconds written by Simulate_AIREML.py.  This writes their mean (the
average time over all estimations) to <out_file>.

Pure standard library, so it runs under any python3 -- no conda needed in the
pipeline's combine job.
"""
import glob
import os
import sys


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 average_time.py <rep_times_dir> <out_file>")
        sys.exit(1)
    rep_dir, out_file = sys.argv[1], sys.argv[2]

    times = []
    for f in sorted(glob.glob(os.path.join(rep_dir, "rep*.txt"))):
        with open(f) as fh:
            s = fh.read().strip()
        if s:
            times.append(float(s))

    if not times:
        print(f"No rep*.txt timing files found in {rep_dir}")
        sys.exit(1)

    avg = sum(times) / len(times)
    out_dir = os.path.dirname(out_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(out_file, "w") as fh:
        fh.write(f"{avg:.4f}\n")
    print(f"Averaged {len(times)} estimations -> {avg:.4f} s  ({out_file})")


if __name__ == "__main__":
    main()
