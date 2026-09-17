# -*- coding: utf-8 -*-
"""
se_plot_G5_G20.py -- fixed-m relative-error (SE) box plots for the pooled-preW
G-sweep folders G5 and G20, both saved into the G5 folder.

Same figure as fixed_m.py (gxg component), but robust to the folder-name vs
file-name mismatch in these folders: the folder is named RandomSNP_pooled_preW_
G<..>_s2gxg0.2_s2e0.8 while the result files inside carry a DIFFERENT basename
(chr1_10ksnp_..._G<..>_..._m10000).  fixed_m.py globs "<foldername>_n*m*.txt" and
would find nothing, so here the real file basename is auto-detected from the
files themselves.  Truth (s2gxg, s2e) is still read from the folder name.

    python se_plot_G5_G20.py
"""
import os
import re
import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")

from plot_script import (read_MoM_results,
                         plot_relative_error_across_groups_combined)
from fixed_m import _parse_truth, _auto_ylimits

HERE = os.path.dirname(os.path.abspath(__file__))

# Source folders to plot, and the folder the PDFs are written into.
# The data in these folders is chr1_10ksnp (m=10000); the folders are named
# "RandomSNP_..." but the files inside are "chr1_10ksnp_...".  To match the
# existing G10 figure (chr1_10ksnp_pooled_preW_G10_..._fixed_m10000_gxg.pdf),
# each output is named with the folder's tokens but "RandomSNP" -> "chr1_10ksnp".
SRC_FOLDERS = [
    "RandomSNP_pooled_preW_G5_s2gxg0.2_s2e0.8",
    "RandomSNP_pooled_preW_G20_s2gxg0.2_s2e0.8",
]
OUT_FOLDER = "result_figure"


def discover_file_basename_and_grid(folder_path):
    """Return (file_basename, sorted n list, single m) from the *actual* result
    files in folder_path, ignoring the folder's own name."""
    n_by_m = {}
    basenames = set()
    for f in glob.glob(os.path.join(folder_path, "*_n*m*.txt")):
        b = os.path.basename(f)
        mobj = re.search(r"^(.*)_n(\d+)m(\d+)\.txt$", b)
        if not mobj:
            continue
        basenames.add(mobj.group(1))
        n_val, m_val = int(mobj.group(2)), int(mobj.group(3))
        n_by_m.setdefault(m_val, []).append(n_val)

    if not n_by_m:
        raise FileNotFoundError(f"No *_n<N>m<M>.txt result files in {folder_path}")
    if len(basenames) != 1:
        raise ValueError(f"Expected one file basename in {folder_path}, got {basenames}")
    if len(n_by_m) != 1:
        raise ValueError(f"Expected a single (fixed) m in {folder_path}, got {sorted(n_by_m)}")

    m_val = next(iter(n_by_m))
    return basenames.pop(), sorted(n_by_m[m_val]), m_val


def main():
    # plot_* saves to os.path.join(os.getcwd(), save_name + ".pdf"); run from
    # HERE so a save_name that starts with OUT_FOLDER lands inside that folder.
    os.chdir(HERE)

    for folder in SRC_FOLDERS:
        folder_path = os.path.join(HERE, folder)
        real_gxg, real_e = _parse_truth(folder)              # from folder name
        file_basename, ns, m_val = discover_file_basename_and_grid(folder_path)

        print(f"\nFolder     : {folder}")
        print(f"File base  : {file_basename}  (folder name differs)")
        print(f"Truth      : s2gxg={real_gxg}, s2e={real_e}")
        print(f"Fixed m    : {m_val}")
        print(f"n values   : {ns}")

        data_dict = read_MoM_results(ns, folder_path, file_basename, m_val)
        real_values = [real_gxg, real_e]

        # gxg component (col 0) -- matches the existing *_fixed_m*_gxg.pdf figures.
        col_num = 0
        ymin, ymax = _auto_ylimits(data_dict, ns, real_values, cols=(col_num,))

        # Name to match the existing G10 figure: folder tokens with the real
        # data type (RandomSNP -> chr1_10ksnp), written into result_figure/.
        fig_stem = folder.replace("RandomSNP", "chr1_10ksnp")
        save_name = os.path.join(OUT_FOLDER, f"{fig_stem}_fixed_m{m_val}_gxg")
        plot_relative_error_across_groups_combined(
            data_dict,
            x_labels=[f"m = {m_val:,}"],
            individual_sizes=ns,
            col_num=col_num,
            real_value=real_values[col_num],
            ymin=ymin,
            ymax=ymax,
            x_axis_name="Sample size (n)",
            save_name=save_name,
        )
        print(f"Saved      : {save_name}.pdf  (y in [{ymin:.3f}, {ymax:.3f}])")


if __name__ == "__main__":
    main()
