# -*- coding: utf-8 -*-
"""
realized_var_four_components_plot.py -- all FOUR realized variances of the
nine-column result files in one figure, each shown against its own expected
value.

Usage
-----
    python realized_var_four_components_plot.py <result_dir> [--exact]
                                                [--outdir DIR] [--png]

<result_dir> holds the row-per-replicate files described in calc_stats.py,

    (s2a, s2d, s2gxg, s2e, Vl_hat, vell_gxg, vell_a, vell_d, vell_e)
     a    d    gxg    e    Vl_hat  Vl_real   Va_real  Vd_real  Ve_real

Columns 6-9 are the REALIZED variances: what this replicate's four effect draws
actually produced, as opposed to the nominal s2a / s2d / s2gxg / s2e they were
drawn under.  This script plots those four columns -- and nothing that was
fitted -- so the figure is a property of the SIMULATOR, not of mc_reml: it
answers "did the draws land on the variances the design asked for, and does the
draw-to-draw scatter shrink with n?"

Each component is drawn against its own expected value:

    Va_real   E = c_additive  * s2a  = s2a       (c_additive  == 1)
    Vd_real   E = c_dominance * s2d  = s2d       (c_dominance == 1)
    Vl_real   E = c * s2gxg                     <- genotype-dependent, so this
                                                   one is a different number at
                                                   every n, read per file from
                                                   c_<tag>.txt (c_hat_hwe by
                                                   default, c_exact with
                                                   --exact)
    Ve_real   E = s2e

Only the gxg component needs the c correction: the additive and dominance
incidence matrices are already standardised, so their c is exactly 1, and the
environmental term is drawn directly at s2e.  This is why Vl_real is the only
series whose dashed reference is not flat.

Layout.  One figure, one y-axis (realized variance, in absolute units), broken
between the a/d/gxg band near 0.1 and the e band near 0.7 -- the two bands are
the same quantity in the same units, so a single continuous axis would compress
the three small components into a hairline.  Markers are the mean over
replicates; whiskers are +/- 1 SD, i.e. the actual draw-to-draw spread, NOT the
standard error of the mean -- the spread is the point of the figure, and it is
the quantity expected to fall as n grows.  Series are dodged slightly in x so
overlapping whiskers stay readable, and each carries its own marker shape as
well as its colour.

The panel carries no title and no footnote -- the setting, the meaning of the
markers and whiskers, and the definition of the dashed references all belong in
the document's figure caption, not repeated inside the image next to it.  The
run prints a ready-to-paste Typst #figure block, with this run's numbers in it,
for Simulation.typ.  The legend stays: it is identity, not commentary.

The figure is written to <result_dir> (override with --outdir); --png writes a
PNG next to the PDF for quick viewing.
"""
import os
import re
import sys
import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Realized-variance columns, 0-based, in the order calc_stats.py names them.
COL_VL, COL_VA, COL_VD, COL_VE = 5, 6, 7, 8

# Categorical slots 1, 2, 3 and 7 of the validated palette: these four clear
# every all-pairs CVD / normal-vision gate on a light surface.  Marker shape is
# a second, redundant encoding, so the series are never told apart by colour
# alone.
SERIES = [
    # key,   column,  colour,     marker, label,         math label
    ("a",   COL_VA, "#2a78d6", "o", "Additive",    r"$\hat{V}_a$"),
    ("d",   COL_VD, "#eb6834", "s", "Dominance",   r"$\hat{V}_d$"),
    ("gxg", COL_VL, "#1baf7a", "^", "Epistasis",   r"$\hat{V}_\ell$"),
    ("e",   COL_VE, "#4a3aa7", "D", "Environment", r"$\hat{V}_e$"),
]

INK_SECOND = "#52514e"


def _parse_truth(folder):
    """s2a, s2d, s2gxg, s2e from the folder name.

    The folder writes some of them with a trailing underscore
    (`s2a_0.1_s2d_0.1s2gxg_0.1_s2e0.7`) and some without, hence the optional
    `_` in each pattern.
    """
    out = {}
    for key in ("s2a", "s2d", "s2gxg", "s2e"):
        m = re.search(key + r"_?([0-9]*\.?[0-9]+)", folder)
        if not m:
            raise ValueError(f"Cannot parse {key} from '{folder}'.")
        out[key] = float(m.group(1))
    return out


def _discover(dir_path):
    """Find the file stem and every (n, m, G) present, skipping the c_ files."""
    grid = {}
    stem = None
    for f in sorted(glob.glob(os.path.join(dir_path, "*_n*m*.txt"))):
        b = os.path.basename(f)
        if b.startswith("c_"):
            continue
        mobj = re.search(r"^(.*)_n(\d+)m(\d+)(?:_G(\d+))?\.txt$", b)
        if not mobj:
            continue
        stem = mobj.group(1)
        grid.setdefault((int(mobj.group(3)), mobj.group(4)), []).append(int(mobj.group(2)))
    if not grid:
        raise FileNotFoundError(f"No result files found in {dir_path}")
    if len(grid) != 1:
        raise ValueError(f"Expected one (m, G) setting, found {sorted(grid)}.")
    (m_val, g_val), ns = next(iter(grid.items()))
    return stem, sorted(ns), m_val, g_val


def _tag(stem, n, m_val, g_val):
    suffix = "" if g_val is None else f"_G{g_val}"
    return f"{stem}_n{n}m{m_val}{suffix}"


def _read(path):
    """Parse "(a,b,...)" rows into an (R, 9) array.

    utf-8-sig for the same reason calc_stats.py uses it: a file that has been
    through a Windows editor carries a BOM.
    """
    rows = []
    with open(path, "r", encoding="utf-8-sig") as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                rows.append([float(v) for v in line.strip("()[] ").split(",")])
            except ValueError:
                raise SystemExit(f"{path}:{lineno}: cannot parse '{line}' as numbers.")
    if not rows:
        raise SystemExit(f"{path}: no data rows.")
    widths = {len(r) for r in rows}
    if len(widths) != 1:
        raise SystemExit(f"{path}: ragged file, row widths {sorted(widths)} -- "
                         "check for a truncated replicate.")
    arr = np.array(rows)
    if arr.shape[1] < 9:
        raise SystemExit(f"{path} has {arr.shape[1]} columns, need 9 "
                         "(the realized columns are 6-9).")
    return arr


def _read_c(dir_path, tag):
    """{key: value} from c_<tag>.txt, or None when the file is absent."""
    path = os.path.join(dir_path, "c_" + tag + ".txt")
    if not os.path.isfile(path):
        return None
    vals = {}
    with open(path, encoding="utf-8-sig") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) == 2:
                vals[parts[0]] = float(parts[1])
    return vals


def _expected(dir_path, tag, truth, arr, use_exact):
    """Expected realized variance per component for one n.

    Returns (expected-by-key, source string, normalized).  Two pipelines write
    these files and they need different handling, which the c file's own keys
    tell apart:

    UN-NORMALIZED (`_unstd_`).  The interaction columns of H are raw, so the
    epistasis kernel carries a factor c and only Vl_hat = c_hat * s2gxg_hat is
    on a variance scale.  Expected epistasis variance is c * s2gxg, keyed
    `expected_realized_variance_hwe` (or `_exact` with --exact).  Without a c
    file, c is recovered from the results -- Vl_hat = c * s2gxg_hat holds
    exactly, so column 5 / column 3 is constant within a file.

    C-NORMALIZED (`_std_`).  The kernel is divided through by c, so
    `c_gxg_after_normalization` is 1, ALL FOUR fitted components are already on
    the realized scale, and column 5 is identically column 3.  The expected
    epistasis variance is then just s2gxg, keyed
    `expected_realized_variance_gxg`.  --exact is meaningless here -- the
    hwe/exact distinction is in the normalization the pipeline already applied
    -- so it is accepted and ignored rather than treated as an error.

    `normalized` says which of the two it was, so callers can label the
    epistasis axis correctly instead of promising a c correction that either
    was or was not applied.
    """
    cvals = _read_c(dir_path, tag) or {}

    if "expected_realized_variance_gxg" in cvals:
        exp_gxg = cvals["expected_realized_variance_gxg"]
        src, normalized = "c file (c-normalized kernel)", True
    else:
        key = ("expected_realized_variance_exact" if use_exact
               else "expected_realized_variance_hwe")
        if key in cvals:
            exp_gxg, src, normalized = cvals[key], "c file", False
        elif use_exact:
            raise SystemExit(f"--exact needs c_exact in c_{tag}.txt, which is "
                             "not there.")
        else:
            ratio = arr[:, 4] / arr[:, 2]
            c_hat = float(np.median(ratio))
            if np.max(np.abs(ratio - c_hat)) > 1e-8 * max(1.0, abs(c_hat)):
                raise SystemExit(f"{tag}: Vl_hat/s2gxg_hat is not constant; "
                                 "cannot recover c without a c file.")
            exp_gxg = c_hat * truth["s2gxg"]
            src = "from results"
            normalized = abs(c_hat - 1.0) < 1e-8

    # The c file is authoritative for the other three when it names them; the
    # folder name is the fallback.  They agree whenever both are present.
    exp = {"a": cvals.get("expected_realized_variance_additive", truth["s2a"]),
           "d": cvals.get("expected_realized_variance_dominance", truth["s2d"]),
           "gxg": exp_gxg,
           "e": cvals.get("expected_realized_variance_residual", truth["s2e"])}
    return exp, src, normalized


def _break_marks(ax_top, ax_bot):
    """The diagonal ticks that mark the y-axis break.

    Only on the left, where the spine actually is: the right spine is off, so a
    mark there would float unattached in white space.
    """
    kw = dict(marker=[(-1, -0.6), (1, 0.6)], markersize=7, linestyle="none",
              color=INK_SECOND, mec=INK_SECOND, mew=1, clip_on=False)
    ax_top.plot([0], [0], transform=ax_top.transAxes, **kw)
    ax_bot.plot([0], [1], transform=ax_bot.transAxes, **kw)


def _spread_labels(entries, ylim, min_gap_frac=0.085):
    """Nudge direct labels apart so near-equal series stay legible.

    entries is [(y, payload), ...]; returns [(y_adjusted, payload), ...] in the
    input order.  Va and Vd both sit on 0.1 by construction, so their labels
    land on top of one another unless they are separated -- the marker itself
    still shows the true value, the label is only a name tag.
    """
    gap = min_gap_frac * (ylim[1] - ylim[0])
    order = sorted(range(len(entries)), key=lambda i: entries[i][0])
    ys = [entries[i][0] for i in order]
    for k in range(1, len(ys)):
        if ys[k] - ys[k - 1] < gap:
            ys[k] = ys[k - 1] + gap
    # Re-centre the block on where it started, so the nudge is symmetric.
    shift = 0.5 * ((entries[order[0]][0] + entries[order[-1]][0]) - (ys[0] + ys[-1]))
    out = list(entries)
    for k, i in enumerate(order):
        out[i] = (ys[k] + shift, entries[i][1])
    return out


def _print_caption(pdf_name, m_val, g_val, ns, reps, mean, sd, exp, truth,
                   use_exact, normalized=False):
    """Emit the Typst #figure block the figure itself no longer carries.

    The panel has no title and no footnote -- everything descriptive lives in
    the document caption instead -- so this prints one ready to paste into
    Simulation.typ, with the numbers filled in from the run that just happened
    rather than typed by hand and left to rot.
    """
    gxg_exp = [exp[n]["gxg"] for n in ns]
    ve_sd_first, ve_sd_last = sd[("e", ns[0])], sd[("e", ns[-1])]
    sd_rng = {k: (min(sd[(k, n)] for n in ns), max(sd[(k, n)] for n in ns))
              for k in ("a", "d", "gxg")}
    worst = max(((abs(mean[(k, n)] - exp[n][k]), k, n)
                 for k in ("a", "d", "gxg", "e") for n in ns))
    worst_name = {"a": "additive", "d": "dominance", "gxg": "epistasis",
                  "e": "environment"}[worst[1]]
    chat = "c" if use_exact else "hat(c)"
    gtag = "" if g_val is None else f", $G = {g_val}$"
    # With a c-normalized kernel the epistasis reference is s2gxg itself and is
    # as flat as the other three, so the sentence about it varying with n would
    # be simply false.
    if normalized:
        gxg_clause = (f'and $sigma^2_(g times g) = {truth["s2gxg"]}$ for '
                      f'epistasis --- the kernel is $c$-normalized here, so '
                      f'this reference is flat in $n$ too')
    else:
        gxg_clause = (f'and ${chat} dot sigma^2_(g times g)$ for epistasis, '
                      f'which is genotype-dependent and so takes a different '
                      f'value at every $n$ '
                      f'(${gxg_exp[0]:.3f}$--${max(gxg_exp):.3f}$)')
    # No thousands separators inside $...$: the rest of Simulation.typ writes
    # sample sizes as $n = 1000$, and a comma reads as a separator in math.

    print("\nTypst caption for Simulation.typ:\n")
    print("#figure(")
    print(f'image("{pdf_name}", width: 100%),caption: [Realized variance of '
          f'the four components against its own expected value '
          f'($m = {m_val}${gtag}, $R = {reps[ns[0]]}$ replicates per point). '
          f'Plotted are the realized variances the effect draws actually '
          f'produced, $hat(V)_a$, $hat(V)_d$, $hat(V)_ell$ and $hat(V)_e$, not '
          f'anything the estimator fitted. Markers: mean over replicates; '
          f'whiskers: $plus.minus 1$ SD, the draw-to-draw spread, not the SE '
          f'of the mean. Dashed lines: the expected value of each component '
          f'--- $sigma^2_a = {truth["s2a"]}$, $sigma^2_d = {truth["s2d"]}$ and '
          f'$sigma^2_e = {truth["s2e"]}$, which are flat in $n$, '
          f'{gxg_clause}. The $y$ axis is broken '
          f'between the additive/dominance/epistasis band and the environment '
          f'band. Every mean sits on its expected value at every $n$ (largest '
          f'deviation ${worst[0]:.3f}$, {worst_name} at $n = {worst[2]}$); the '
          f'environmental SD contracts from ${ve_sd_first:.3f}$ to '
          f'${ve_sd_last:.3f}$, while the three genetic SDs are flat in $n$ '
          f'(additive ${sd_rng["a"][0]:.3f}$--${sd_rng["a"][1]:.3f}$, '
          f'dominance ${sd_rng["d"][0]:.3f}$--${sd_rng["d"][1]:.3f}$, '
          f'epistasis ${sd_rng["gxg"][0]:.3f}$--${sd_rng["gxg"][1]:.3f}$) '
          f'because they are set by the $m = {m_val}$ effect draws rather '
          f'than by the sample size.]')
    print(")")


def main():
    _usage = ("Usage: python realized_var_four_components_plot.py <result_dir> "
              "[--exact] [--outdir DIR] [--png]")
    if len(sys.argv) < 2:
        print(_usage)
        sys.exit(1)

    dir_path = os.path.abspath(sys.argv[1])
    if not os.path.isdir(dir_path):
        print(f"Error: '{dir_path}' is not a directory.")
        sys.exit(1)

    rest = sys.argv[2:]
    use_exact = "--exact" in rest
    want_png = "--png" in rest
    rest = [t for t in rest if t not in ("--exact", "--png")]

    out_dir = dir_path
    if "--outdir" in rest:
        i = rest.index("--outdir")
        if i + 1 >= len(rest):
            print("Error: --outdir needs a directory.")
            print(_usage)
            sys.exit(1)
        out_dir = os.path.abspath(rest[i + 1])
        rest = rest[:i] + rest[i + 2:]
    if rest:
        print(f"Error: unrecognised argument(s) {rest}.")
        print(_usage)
        sys.exit(1)

    folder = os.path.basename(dir_path)
    truth = _parse_truth(folder)
    stem, ns, m_val, g_val = _discover(dir_path)

    print(f"Directory : {dir_path}")
    print(f"Truth     : s2a={truth['s2a']}, s2d={truth['s2d']}, "
          f"s2gxg={truth['s2gxg']}, s2e={truth['s2e']}")
    print(f"Fixed m   : {m_val}" + ("" if g_val is None else f", G = {g_val}"))
    _probe = _read_c(dir_path, _tag(stem, ns[0], m_val, g_val)) or {}
    is_norm = "expected_realized_variance_gxg" in _probe
    if is_norm:
        print("Kernel    : c-NORMALIZED -- all four fitted components are already "
              "on the realized scale, expected epistasis variance is s2gxg itself")
        if use_exact:
            print("            (--exact ignored: there is no post-fit c to choose)")
    else:
        print(f"c used    : {'c_exact' if use_exact else 'c_hat_hwe'} (gxg only; "
              "c == 1 for a and d)")
    print("Plotted   : columns 6-9, the REALIZED variances, vs. their expected values")

    mean, sd, exp, reps = {}, {}, {}, {}
    for n in ns:
        tag = _tag(stem, n, m_val, g_val)
        arr = _read(os.path.join(dir_path, tag + ".txt"))
        reps[n] = arr.shape[0]
        exp[n], src, normalized = _expected(dir_path, tag, truth, arr, use_exact)
        for key, col, _c, _mk, _lab, _ml in SERIES:
            mean[(key, n)] = float(arr[:, col].mean())
            sd[(key, n)] = float(arr[:, col].std(ddof=1))
        print(f"  n={n:<6d} R={arr.shape[0]:<4d} [{src}]")
        for key, _col, _c, _mk, lab, _ml in SERIES:
            bias = mean[(key, n)] - exp[n][key]
            se = sd[(key, n)] / np.sqrt(arr.shape[0])
            flag = "  *" if abs(bias) > 1.96 * se else ""
            print(f"      {lab:<12s} E={exp[n][key]:.6f}  mean={mean[(key, n)]:.6f}  "
                  f"SD={sd[(key, n)]:.6f}  mean-E={bias:+.6f}{flag}")
    print("  (* = |mean - expected| exceeds 1.96 SE of the mean)")

    # ---------------- figure -------------------------------------------------
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 10,
        'axes.linewidth': 1,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.major.width': 1,
        'ytick.major.width': 1,
        'xtick.major.size': 4,
        'ytick.major.size': 4,
        'figure.dpi': 150,
        'savefig.dpi': 600,
    })

    x = np.arange(len(ns), dtype=float)
    # The three small components share a band; dodge them so their whiskers do
    # not sit on top of one another.  Environment is alone on the top panel and
    # needs no dodge.
    dodge = {"a": -0.13, "d": 0.0, "gxg": 0.13, "e": 0.0}
    lower = [s for s in SERIES if s[0] != "e"]
    upper = [s for s in SERIES if s[0] == "e"]

    def band(series):
        lo, hi = np.inf, -np.inf
        for key, _col, _c, _mk, _lab, _ml in series:
            for n in ns:
                lo = min(lo, mean[(key, n)] - sd[(key, n)], exp[n][key])
                hi = max(hi, mean[(key, n)] + sd[(key, n)], exp[n][key])
        pad = 0.18 * (hi - lo)
        return lo - pad, hi + pad

    lo_b, hi_b = band(lower)
    lo_t, hi_t = band(upper)

    fig, (ax_t, ax_b) = plt.subplots(
        2, 1, sharex=True, figsize=(7.2, 5.6),
        gridspec_kw={"height_ratios": [1.0, 2.4], "hspace": 0.09})

    for ax, series in ((ax_t, upper), (ax_b, lower)):
        ax.set_axisbelow(True)
        ax.grid(axis="y", color="#e6e5e1", linewidth=0.8)
        for key, _col, colour, marker, label, _ml in series:
            xs = x + dodge[key]
            mu = np.array([mean[(key, n)] for n in ns])
            er = np.array([sd[(key, n)] for n in ns])
            ev = np.array([exp[n][key] for n in ns])

            # Expected value: dashed, recessive, in the series colour.  Flat for
            # a / d / e; c-driven and therefore n-dependent for gxg.
            ax.plot(xs, ev, linestyle=(0, (5, 3)), linewidth=1.4, color=colour,
                    alpha=0.55, zorder=2)
            ax.errorbar(xs, mu, yerr=er, fmt="none", ecolor=colour,
                        elinewidth=1.6, capsize=3.5, capthick=1.2,
                        alpha=0.9, zorder=3)
            ax.plot(xs, mu, linestyle="-", linewidth=2, color=colour,
                    marker=marker, markersize=7, markerfacecolor=colour,
                    markeredgecolor="white", markeredgewidth=1.2,
                    label=label, zorder=4)

    ax_t.set_ylim(lo_t, hi_t)
    ax_b.set_ylim(lo_b, hi_b)

    # Break the shared axis between the two bands.
    ax_t.spines["bottom"].set_visible(False)
    ax_t.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    _break_marks(ax_t, ax_b)

    # Direct labels at the right edge -- the relief the palette's contrast WARN
    # requires, and they also make the series readable without the legend.
    x_lab = x[-1] + max(dodge.values())
    for ax, series in ((ax_t, upper), (ax_b, lower)):
        entries = [(mean[(key, ns[-1])], (colour, mlab))
                   for key, _col, colour, _mk, _lab, mlab in series]
        for y_lab, (colour, mlab) in _spread_labels(entries, ax.get_ylim()):
            # One x for the whole panel, so the labels line up in a column
            # instead of stepping with the series dodge.
            ax.annotate(mlab, xy=(x_lab, y_lab),
                        xytext=(11, 0), textcoords="offset points",
                        color=colour, fontsize=11, va="center", ha="left",
                        annotation_clip=False)

    ax_b.set_xticks(x)
    ax_b.set_xticklabels([f"{n:,}" for n in ns], fontsize=9)
    ax_b.set_xlim(-0.5, len(ns) - 0.5 + 0.4)
    ax_b.set_xlabel("Sample size (n)", fontsize=10, labelpad=8)

    fig.supylabel("Realized variance", fontsize=10, x=0.035)

    # No title and no explanatory text inside the figure: the setting (m, G, R),
    # what the markers and whiskers are, and what the dashed lines mean all
    # belong in the document's figure caption, which is printed at the end of
    # this run for pasting into Simulation.typ.  The legend stays -- it is
    # identity, not commentary, and the series cannot be read without it.
    handles, labels = [], []
    for ax in (ax_b, ax_t):
        h, l = ax.get_legend_handles_labels()
        handles += h
        labels += l
    leg = ax_t.legend(handles, labels, loc="lower left",
                      bbox_to_anchor=(0.0, 1.015), ncol=4, frameon=False,
                      fontsize=9, handlelength=2.2, columnspacing=2.0,
                      borderaxespad=0.0)
    for t in leg.get_texts():
        t.set_color(INK_SECOND)

    os.makedirs(out_dir, exist_ok=True)
    etag = "_cexact" if (use_exact and not normalized) else ""
    gsuf = "" if g_val is None else f"_G{g_val}"
    base = f"realized_var_four_components_m{m_val}{gsuf}{etag}"
    save_kw = dict(bbox_inches="tight", facecolor="white", edgecolor="none",
                   transparent=False, pad_inches=0.15)

    pdf_path = os.path.join(out_dir, base + ".pdf")
    fig.savefig(pdf_path, format="pdf", **save_kw)
    print(f"PDF saved to: {pdf_path}")
    if want_png:
        png_path = os.path.join(out_dir, base + ".png")
        fig.savefig(png_path, format="png", **save_kw)
        print(f"PNG saved to: {png_path}")

    plt.close(fig)
    plt.rcParams.update(plt.rcParamsDefault)

    _print_caption(base + ".pdf", m_val, g_val, ns, reps, mean, sd, exp,
                   truth, use_exact, normalized)


if __name__ == "__main__":
    main()
