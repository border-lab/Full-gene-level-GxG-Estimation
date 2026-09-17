===============================================================================
 Plot_code  --  figures for the epistasis-only (GxG) MCREML simulations
===============================================================================

Files
-----
  fixed_m.py       Command-line tool: fixed-m, increasing-n box-plot figure.
  plot_script.py   Library of plotting / result-reading functions used by it.
  <result_dir>/    One folder per simulation setting (see "Input" below).


-------------------------------------------------------------------------------
 fixed_m.py
-------------------------------------------------------------------------------
Makes a box-plot of the REML relative error (estimate - truth) with the SNP
count m held fixed and the sample size n increasing left-to-right.

  Usage:
      python fixed_m.py <result_dir> [component] [ymin ymax]

  NOTE: use "python", not "python3".  On Windows "python3" is often a
  Microsoft Store app-execution-alias stub that prints nothing and exits
  (code 49) without running the script -- so the command appears to "return
  nothing".  This project's real interpreter (Anaconda) is on "python".

Arguments
  <result_dir>   (required) Folder of results, e.g. Random_s2gxg0.05_s2e0.95.
                 The true s2gxg / s2e are read from the folder name, and the
                 fixed m and the list of n are discovered from the files in it.

  component      (optional) Which variance to plot. Default: both.
                   both   two stacked panels: s2gxg on top, s2e below
                   gxg    single panel, s2gxg only
                   e      single panel, s2e only

  ymin ymax      (optional) Two numbers that fix the y-axis range,
                 e.g.  -0.5 0.5.  If omitted, the range is auto-fitted to the
                 data.  ymin must be less than ymax.

  component and "ymin ymax" may be given in any order after the folder.


Examples
  # both components, y-axis auto
  python fixed_m.py Random_s2gxg0.05_s2e0.95

  # s2gxg only, y-axis auto
  python fixed_m.py Random_s2gxg0.05_s2e0.95 gxg

  # s2e only, y-axis auto
  python fixed_m.py Random_s2gxg0.05_s2e0.95 e

  # s2gxg only, y-axis fixed to [-0.5, 0.5]
  python fixed_m.py Random_s2gxg0.05_s2e0.95 gxg -0.5 0.5

  # both components, y-axis fixed to [-0.4, 0.4]
  python fixed_m.py Random_s2gxg0.05_s2e0.95 -0.4 0.4


Input: the result folder
  Folder name encodes the truth:   <basename>_s2gxg<VALUE>_s2e<VALUE>
  Files inside are named:          <basename>_n<N>m<M>.txt
  Each file has one row per replicate, formatted  "(s2gxg_hat, s2e_hat)".
  All files in a folder must share the same m (one fixed m per figure).

  Example:
      Random_s2gxg0.05_s2e0.95/
          Random_s2gxg0.05_s2e0.95_n1000m1000.txt
          Random_s2gxg0.05_s2e0.95_n2000m1000.txt
          ...
          Random_s2gxg0.05_s2e0.95_n32000m1000.txt


Output
  A PDF written to the current directory:
      both :  <basename>_fixed_m<M>.pdf
      gxg  :  <basename>_fixed_m<M>_gxg.pdf
      e    :  <basename>_fixed_m<M>_e.pdf
  Each box shows the mean, SD, and a one-sample t-test vs 0 (ns / * / ** / ***).
  The script prints the truth, m, n values, chosen component, and y-axis range
  so you can confirm what was plotted.
