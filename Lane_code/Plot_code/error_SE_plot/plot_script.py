import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

def read_MoM_results(individual_sizes, path, file_name, num_snp):
    """
    read_MoM_results from cluster output files.

    parameters:
    individual_sizes: list of individual sizes (e.g., [1000, 2000, 4000])
    path: directory path where the result files are stored
    num_snp: number of SNPs used in the analysis
    file_name: base name of the result files (e.g., "MoM_results")

    returns:
    A dictionary with individual sizes as keys and corresponding DataFrames as values.
    """ 

    data_dict = {}
    for n in individual_sizes:
        full_file_name = f"{file_name}_n{n}m{num_snp}.txt"  # different variable
        file_path = os.path.join(path, full_file_name)
        
        df = pd.read_csv(file_path, header=None)
        df[0] = df[0].astype(str).str.replace('(', '', regex=False).astype(float)
        df[1] = df[1].astype(str).str.replace(')', '', regex=False).astype(float)
        
        data_dict[n] = df
        
    return data_dict

def read_MoM_results_add(individual_sizes, path, file_name, num_snp):
    """
    read_MoM_results from cluster output files with three columns (a, gxg, e).
    parameters:
    individual_sizes: list of individual sizes (e.g., [1000, 2000, 4000])
    path: directory path where the result files are stored
    num_snp: number of SNPs used in the analysis
    file_name: base name of the result files (e.g., "MoM_results")
    returns:
    A dictionary with individual sizes as keys and corresponding DataFrames as values.
    Each DataFrame has three columns: [0]=a (additive), [1]=gxg (epistatic), [2]=e (error).
    """
    data_dict = {}
    for n in individual_sizes:
        full_file_name = f"{file_name}_n{n}m{num_snp}.txt"
        file_path = os.path.join(path, full_file_name)

        df = pd.read_csv(file_path, header=None)
        df[0] = df[0].astype(str).str.replace('(', '', regex=False).astype(float)
        df[1] = df[1].astype(str).astype(float)
        df[2] = df[2].astype(str).str.replace(')', '', regex=False).astype(float)

        data_dict[n] = df

    return data_dict
    
def _nice_tick_step(span, target=6):
    """A round tick step (1/2/2.5/5 x 10^k) giving roughly `target` ticks."""
    if span <= 0 or not np.isfinite(span):
        return 0.25
    raw = span / target
    mag = 10.0 ** np.floor(np.log10(raw))
    for mult in (1, 2, 2.5, 5, 10):
        if mult * mag >= raw:
            return mult * mag
    return 10 * mag


def plot_relative_error_across_groups_combined(*data_dicts, x_labels, individual_sizes, col_num, real_value, ymin, ymax, x_axis_name="Group", title=None, save_name=None, custom_bottom_labels=None, y_axis_name=None, ytick_step=None, p_values=None, show_significance=True, ref_line=0.0):
    """
    Plot relative errors with box plots in journal style.

    Parameters:
    -----------
    col_num : int
        Column index to plot.
        0 = additive (a), 1 = epistatic (gxg), 2 = error (e)
        For two-column data: 0 = gxg, 1 = e
    save_name : str
        Filename without extension (e.g., 'my_plot').
        PDF will be saved to the same directory as the script.
        If None, plot is only displayed.
    custom_bottom_labels : list of str, optional
        Custom labels for the x-axis (one per box). If None, defaults to
        "n = {value}" using individual_sizes. Length must equal
        len(x_labels) * len(individual_sizes).
    x_axis_name : str
        Label for the x-axis. Defaults to "Group" but you can pass e.g.
        "Block size", "Number of SNPs (m)", etc.
    y_axis_name : str, optional
        Full label for the y-axis. If provided, it overrides the default
        "Relative error (<parameter symbol>)" label entirely. If None, the
        label is auto-generated from col_num as before.
    """

    # Set high-quality rendering
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

    # Gather data for all combinations
    data_list = []
    combined_labels = []
    bottom_labels = []

    for label, data_dict in zip(x_labels, data_dicts):
        for n in individual_sizes:
            df = data_dict[n]
            col_values = df.iloc[:, col_num].values - real_value
            combined_labels.append(f"{label}\nn={n:,}")
            bottom_labels.append(f"n = {n:,}")
            data_list.append(pd.DataFrame({"value": col_values, "group": combined_labels[-1]}))

    # Override bottom labels if custom ones are provided
    if custom_bottom_labels is not None:
        if len(custom_bottom_labels) != len(bottom_labels):
            raise ValueError(
                f"custom_bottom_labels has {len(custom_bottom_labels)} entries "
                f"but {len(bottom_labels)} boxes are being plotted "
                f"(len(x_labels)={len(x_labels)} * len(individual_sizes)={len(individual_sizes)})."
            )
        bottom_labels = list(custom_bottom_labels)

    data = pd.concat(data_list, ignore_index=True)
    data["group"] = pd.Categorical(data["group"], categories=combined_labels, ordered=True)

    # Compute summary statistics
    summary = (
        data.groupby("group", observed=True)["value"]
        .agg(["mean", "std", "count"])
        .loc[combined_labels]
    )

    # Significance per box.  By default a one-sample t-test of the plotted
    # values against 0.  A caller whose reference point is itself estimated
    # from the same replicates must supply the correctly-paired p-values
    # instead, one per box in plotting order -- the default test would use the
    # wrong standard error and over-report significance.
    if not show_significance:
        pvals = {label: np.nan for label in combined_labels}
    elif p_values is None:
        pvals = {}
        for combined_label in combined_labels:
            group_data = data[data["group"] == combined_label]["value"].values
            t_stat, p_val = stats.ttest_1samp(group_data, 0)
            pvals[combined_label] = p_val
    else:
        if len(p_values) != len(combined_labels):
            raise ValueError(
                f"p_values has {len(p_values)} entries but {len(combined_labels)} "
                "boxes are being plotted.")
        pvals = dict(zip(combined_labels, p_values))

    summary["p_value"] = summary.index.map(pvals)

    # x-axis positions
    x_positions = np.arange(len(summary))

    # Create figure
    fig, ax = plt.subplots(figsize=(3 + len(combined_labels) * 0.8, 5))

    # Extend y-axis for labels above plot
    ymax_extended = ymax + 0.25 * (ymax - ymin)
    ax.set_ylim(ymin, ymax_extended)
    ax.set_xlim(-0.6, len(combined_labels) - 0.4)

    # Add alternating white/gray background for EACH BOX
    for i in range(len(combined_labels)):
        x_start = i - 0.5
        x_end = i + 0.5
        if i % 2 == 0:
            ax.axvspan(x_start, x_end, facecolor='white', alpha=1.0, zorder=0)
        else:
            ax.axvspan(x_start, x_end, facecolor='#E8E8E8', alpha=0.8, zorder=0)

    # Prepare data for box plots
    box_data = [data[data["group"] == label]["value"].values for label in combined_labels]

    # Define colors
    box_color = '#3274A1'
    median_color = '#CC0000'
    label_color = '#000000'

    # Create box plots - NO OUTLIERS
    bp = ax.boxplot(
        box_data,
        positions=x_positions,
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        boxprops=dict(linewidth=1.5, edgecolor=box_color, facecolor='white'),
        whiskerprops=dict(linewidth=1.2, color=box_color),
        capprops=dict(linewidth=1.2, color=box_color),
        medianprops=dict(linewidth=2, color=median_color)
    )

    # Add m = label at top
    for i, label in enumerate(x_labels):
        group_center = (i + 0.5) * len(individual_sizes) - 0.5
        ax.text(
            group_center, ymax_extended - 0.01 * (ymax_extended - ymin),
            label,
            ha='center', va='top', fontsize=11, fontweight='bold', color=label_color
        )

    # Significance stars, Mean, and SD labels.  A figure with no reference point
    # has no null hypothesis, so show_significance=False drops the stars and
    # moves Mean / SD up into the space they vacate.
    y_mean, y_sd = (0.15, 0.22) if show_significance else (0.08, 0.15)
    for i, (combined_label, row) in enumerate(summary.iterrows()):
        if show_significance:
            p_val = row['p_value']
            if p_val < 0.001:
                sig_stars = '***'
            elif p_val < 0.01:
                sig_stars = '**'
            elif p_val < 0.05:
                sig_stars = '*'
            else:
                sig_stars = 'ns'

            ax.text(
                i, ymax_extended - 0.08 * (ymax_extended - ymin),
                sig_stars,
                ha='center', va='top', fontsize=11, fontweight='bold',
                color=label_color
            )

        ax.text(
            i, ymax_extended - y_mean * (ymax_extended - ymin),
            f"Mean={row['mean']:.3f}",
            ha='center', va='top', fontsize=9, fontweight='bold',
            color=label_color
        )

        ax.text(
            i, ymax_extended - y_sd * (ymax_extended - ymin),
            f"SD={row['std']:.3f}",
            ha='center', va='top', fontsize=9, fontweight='bold',
            color=label_color
        )

    # Determine y-axis label. If y_axis_name is provided, use it verbatim.
    # Otherwise fall back to the auto-generated label based on col_num.
    if y_axis_name is not None:
        ylabel = y_axis_name
    else:
        first_df = list(data_dicts[0].values())[0]
        num_cols = first_df.shape[1]

        if num_cols == 3:
            if col_num == 0:
                theta_simple = r"$\sigma^2_{a}$"
            elif col_num == 1:
                theta_simple = r"$\sigma^2_{g \times g}$"
            elif col_num == 2:
                theta_simple = r"$\sigma^2_{e}$"
            else:
                theta_simple = "Parameter"
        else:
            if col_num == 0:
                theta_simple = r"$\sigma^2_{g \times g}$"
            elif col_num == 1:
                theta_simple = r"$\sigma^2_{e}$"
            else:
                theta_simple = "Parameter"

        ylabel = f"Relative error ({theta_simple})"

    # Reference line.  ref_line=None draws none, for raw-value figures where
    # zero is not a meaningful reference.
    if ref_line is not None:
        ax.axhline(ref_line, color='#666666', linestyle='--', linewidth=0.8, zorder=1)

    # Axis labels (x-axis label uses x_axis_name, y-axis label uses ylabel)
    ax.set_xlabel(x_axis_name, fontsize=10, labelpad=8)
    ax.set_ylabel(ylabel, fontsize=10, labelpad=8)

    # Y-axis ticks.  A fixed 0.25 step leaves narrow-range figures with a
    # single labelled tick, so the step adapts to the range unless one is given.
    step = ytick_step if ytick_step else _nice_tick_step(ymax - ymin)
    k0, k1 = int(np.ceil(ymin / step)), int(np.floor(ymax / step))
    yticks = [round(k * step, 10) for k in range(k0, k1 + 1)]
    ax.set_yticks(yticks)

    # X-tick labels (uses bottom_labels, which may be custom)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(bottom_labels, fontsize=9)

    # Title
    if title:
        ax.set_title(title, fontsize=11, pad=10, fontweight='normal')

    # Clean up spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#333333')
    ax.spines['bottom'].set_color('#333333')

    ax.tick_params(axis='both', which='major', labelsize=9, colors='#333333')

    plt.tight_layout()

    # Save PDF to script directory
    if save_name:
        script_dir = os.getcwd()
        save_path = os.path.join(script_dir, f"{save_name}.pdf")

        plt.savefig(
            save_path,
            bbox_inches='tight',
            facecolor='white',
            edgecolor='none',
            format='pdf',
            transparent=False,
            pad_inches=0.1
        )
        print(f"PDF saved to: {save_path}")

    plt.show()

    plt.rcParams.update(plt.rcParamsDefault)

def plot_relative_error_combined_vertical(data_dicts_list, x_labels_list, individual_sizes, col_num, real_value, ymin, ymax, save_name=None):
    """
    Plot multiple relative error box plots as vertical subplots in one figure.
    
    Parameters:
    -----------
    data_dicts_list : list of data_dicts
        e.g., [Contiguous_raw_1k, Contiguous_raw_2k, Contiguous_raw_4k]
    x_labels_list : list of x_labels
        e.g., [["m = 1000"], ["m = 2000"], ["m = 4000"]]
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import os
    from scipy import stats
    
    # Set high-quality rendering
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
    
    num_plots = len(data_dicts_list)
    # Flattened figure: wider and shorter per subplot
    fig, axes = plt.subplots(num_plots, 1, figsize=(10, 3 * num_plots), sharex=True)
    
    if num_plots == 1:
        axes = [axes]
    
    for ax_idx, (data_dict, x_labels) in enumerate(zip(data_dicts_list, x_labels_list)):
        ax = axes[ax_idx]
        
        # Gather data for all combinations
        data_list = []
        combined_labels = []
        bottom_labels = []
        
        for label in x_labels:
            for n in individual_sizes:
                df = data_dict[n]
                col_values = df.iloc[:, col_num].values - real_value
                combined_labels.append(f"{label}\nn={n:,}")
                bottom_labels.append(f"n = {n:,}")
                data_list.append(pd.DataFrame({"value": col_values, "group": combined_labels[-1]}))
        
        data = pd.concat(data_list, ignore_index=True)
        data["group"] = pd.Categorical(data["group"], categories=combined_labels, ordered=True)
        
        # Compute summary statistics
        summary = (
            data.groupby("group", observed=True)["value"]
            .agg(["mean", "std", "count"])
            .loc[combined_labels]
        )
        
        # Perform one-sample t-test for each group (mean != 0)
        p_values = {}
        for combined_label in combined_labels:
            group_data = data[data["group"] == combined_label]["value"].values
            t_stat, p_val = stats.ttest_1samp(group_data, 0)
            p_values[combined_label] = p_val
        
        summary["p_value"] = summary.index.map(p_values)
        
        # x-axis positions
        x_positions = np.arange(len(summary))
        
        # Extend y-axis for labels above plot
        ymax_extended = ymax + 0.25 * (ymax - ymin)
        ax.set_ylim(ymin, ymax_extended)
        ax.set_xlim(-0.6, len(combined_labels) - 0.4)
        
        # Add alternating white/gray background for EACH BOX
        for i in range(len(combined_labels)):
            x_start = i - 0.5
            x_end = i + 0.5
            if i % 2 == 0:
                ax.axvspan(x_start, x_end, facecolor='white', alpha=1.0, zorder=0)
            else:
                ax.axvspan(x_start, x_end, facecolor='#E8E8E8', alpha=0.8, zorder=0)
        
        # Prepare data for box plots
        box_data = [data[data["group"] == label]["value"].values for label in combined_labels]
        
        # Define colors
        box_color = '#3274A1'
        median_color = '#CC0000'
        label_color = '#000000'
        
        # Create box plots - NO OUTLIERS
        bp = ax.boxplot(
            box_data,
            positions=x_positions,
            widths=0.5,
            patch_artist=True,
            showfliers=False,
            boxprops=dict(linewidth=1.5, edgecolor=box_color, facecolor='white'),
            whiskerprops=dict(linewidth=1.2, color=box_color),
            capprops=dict(linewidth=1.2, color=box_color),
            medianprops=dict(linewidth=2, color=median_color)
        )
        
        # Add m = label at top (as title)
        for i, label in enumerate(x_labels):
            group_center = (i + 0.5) * len(individual_sizes) - 0.5
            ax.text(
                group_center, ymax_extended - 0.01 * (ymax_extended - ymin),
                label,
                ha='center', va='top', fontsize=11, fontweight='bold', color=label_color
            )
        
        # Significance stars, Mean, and SD labels
        for i, (combined_label, row) in enumerate(summary.iterrows()):
            p_val = row['p_value']
            if p_val < 0.001:
                sig_stars = '***'
            elif p_val < 0.01:
                sig_stars = '**'
            elif p_val < 0.05:
                sig_stars = '*'
            else:
                sig_stars = 'ns'
            
            ax.text(
                i, ymax_extended - 0.08 * (ymax_extended - ymin),
                sig_stars,
                ha='center', va='top', fontsize=10, fontweight='bold',
                color=label_color
            )
            
            ax.text(
                i, ymax_extended - 0.15 * (ymax_extended - ymin),
                f"Mean={row['mean']:.3f}",
                ha='center', va='top', fontsize=8, fontweight='bold',
                color=label_color
            )
            
            ax.text(
                i, ymax_extended - 0.22 * (ymax_extended - ymin),
                f"SD={row['std']:.3f}",
                ha='center', va='top', fontsize=8, fontweight='bold',
                color=label_color
            )
        
        # Theta label for y-axis
        if col_num == 0:
            theta_simple = r"$\sigma^2_{g \times g}$"
        elif col_num == 1:
            theta_simple = r"$\sigma^2_{e}$"
        else:
            theta_simple = "Parameter"
        
        # Reference line at zero
        ax.axhline(0, color='#666666', linestyle='--', linewidth=0.8, zorder=1)
        
        # Y-axis label
        ax.set_ylabel(f"Relative error ({theta_simple})", fontsize=10, labelpad=8)
        
        # Y-axis ticks
        yticks = np.arange(np.ceil(ymin * 4) / 4, ymax + 0.01, 0.25)
        yticks = [y for y in yticks if y <= ymax]
        ax.set_yticks(yticks)
        
        # X-tick labels (only for bottom plot)
        ax.set_xticks(x_positions)
        if ax_idx == num_plots - 1:
            ax.set_xticklabels(bottom_labels, fontsize=9)
            ax.set_xlabel("Sample size (n)", fontsize=10, labelpad=8)
        else:
            ax.set_xticklabels([])
        
        # Panel label (A), (B), (C)
        panel_label = chr(65 + ax_idx)  # A, B, C, ...
        ax.text(-0.08, 1.02, f"({panel_label})", transform=ax.transAxes, 
                fontsize=12, fontweight='bold', va='bottom', ha='left')
        
        # Clean up spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#333333')
        ax.spines['bottom'].set_color('#333333')
        
        ax.tick_params(axis='both', which='major', labelsize=9, colors='#333333')
    
    plt.tight_layout()
    
    # Save PDF to script directory
    if save_name:
        script_dir = os.getcwd()
        save_path = os.path.join(script_dir, f"{save_name}.pdf")
        
        plt.savefig(
            save_path,
            bbox_inches='tight',
            facecolor='white',
            edgecolor='none',
            format='pdf',
            transparent=False,
            pad_inches=0.1
        )
        print(f"PDF saved to: {save_path}")
    
    plt.show()
    
    # Reset rcParams
    plt.rcParams.update(plt.rcParamsDefault)


def plot_three_effects_single_m(data_dict, individual_sizes, real_values, m_label, ymin, ymax, save_name=None):
    """
    Plot relative errors for three variance components (a, gxg, e) in stacked panels.
    Compact version with reduced vertical spacing.
    """
    from scipy import stats
    
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 9,
        'axes.linewidth': 1,
        'figure.dpi': 150,
        'savefig.dpi': 600,
    })
    
    components = [
        (0, r"$\sigma^2_{a}$"),
        (1, r"$\sigma^2_{g \times g}$"),
        (2, r"$\sigma^2_{e}$"),
    ]
    
    num_boxes = len(individual_sizes)
    
    # Compact figure: each panel is shorter
    fig, axes = plt.subplots(
        3, 1,
        figsize=(3 + num_boxes * 0.7, 7.5),  # Reduced from 12 to 7.5
        sharex=True,
        gridspec_kw={'hspace': 0.25}  # Tighter spacing between panels
    )
    
    box_color = '#3274A1'
    median_color = '#CC0000'
    label_color = '#000000'
    
    ymax_extended = ymax + 0.30 * (ymax - ymin)  # Slightly more room for labels
    
    for ax_idx, (col_num, theta_label) in enumerate(components):
        ax = axes[ax_idx]
        real_value = real_values[col_num]
        
        box_data = []
        means = []
        sds = []
        p_values = []
        bottom_labels = []
        
        for n in individual_sizes:
            df = data_dict[n]
            col_values = df.iloc[:, col_num].values - real_value
            box_data.append(col_values)
            means.append(np.mean(col_values))
            sds.append(np.std(col_values, ddof=1))
            _, p_val = stats.ttest_1samp(col_values, 0)
            p_values.append(p_val)
            bottom_labels.append(f"n = {n:,}")
        
        x_positions = np.arange(len(individual_sizes))
        
        ax.set_ylim(ymin, ymax_extended)
        ax.set_xlim(-0.6, len(individual_sizes) - 0.4)
        
        # Alternating background
        for i in range(len(individual_sizes)):
            color = 'white' if i % 2 == 0 else '#E8E8E8'
            ax.axvspan(i - 0.5, i + 0.5, facecolor=color, alpha=1.0 if i % 2 == 0 else 0.8, zorder=0)
        
        # Box plots
        ax.boxplot(
            box_data,
            positions=x_positions,
            widths=0.5,
            patch_artist=True,
            showfliers=False,
            boxprops=dict(linewidth=1.2, edgecolor=box_color, facecolor='white'),
            whiskerprops=dict(linewidth=1, color=box_color),
            capprops=dict(linewidth=1, color=box_color),
            medianprops=dict(linewidth=1.5, color=median_color)
        )
        
        # m = label at top of top panel only
        if ax_idx == 0:
            ax.text(
                (len(individual_sizes) - 1) / 2,
                ymax_extended - 0.01 * (ymax_extended - ymin),
                m_label,
                ha='center', va='top', fontsize=10, fontweight='bold', color=label_color
            )
        
        # Offset for labels
        top_offset = 0.09 if ax_idx == 0 else 0.03
        
        for i, (mean_val, sd_val, p_val) in enumerate(zip(means, sds, p_values)):
            if p_val < 0.001:
                sig_stars = '***'
            elif p_val < 0.01:
                sig_stars = '**'
            elif p_val < 0.05:
                sig_stars = '*'
            else:
                sig_stars = 'ns'
            
            ax.text(
                i, ymax_extended - top_offset * (ymax_extended - ymin),
                sig_stars,
                ha='center', va='top', fontsize=9, fontweight='bold', color=label_color
            )
            ax.text(
                i, ymax_extended - (top_offset + 0.08) * (ymax_extended - ymin),
                f"Mean={mean_val:.3f}",
                ha='center', va='top', fontsize=7.5, fontweight='bold', color=label_color
            )
            ax.text(
                i, ymax_extended - (top_offset + 0.16) * (ymax_extended - ymin),
                f"SD={sd_val:.3f}",
                ha='center', va='top', fontsize=7.5, fontweight='bold', color=label_color
            )
        
        # Reference line at zero
        ax.axhline(0, color='#666666', linestyle='--', linewidth=0.8, zorder=1)
        
        # Y-axis label
        ax.set_ylabel(f"Relative error ({theta_label})", fontsize=9, labelpad=6)
        
        # Y-axis ticks: fewer ticks for compactness
        yticks = np.arange(np.ceil(ymin * 2) / 2, ymax + 0.01, 0.5)
        yticks = [y for y in yticks if y <= ymax]
        ax.set_yticks(yticks)
        
        # X-axis (only bottom)
        ax.set_xticks(x_positions)
        if ax_idx == len(components) - 1:
            ax.set_xticklabels(bottom_labels, fontsize=8)
            ax.set_xlabel("Sample size (n)", fontsize=9, labelpad=6)
        else:
            ax.set_xticklabels([])
        
        # Spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#333333')
        ax.spines['bottom'].set_color('#333333')
        ax.tick_params(axis='both', which='major', labelsize=8, colors='#333333')
    

    
    if save_name:
        script_dir = os.getcwd()
        save_path = os.path.join(script_dir, f"{save_name}.pdf")
        plt.savefig(save_path, bbox_inches='tight', facecolor='white',
                    edgecolor='none', format='pdf', transparent=False, pad_inches=0.1)
        print(f"PDF saved to: {save_path}")
    
    plt.show()
    plt.rcParams.update(plt.rcParamsDefault)


def plot_two_effects_single_m(data_dict, individual_sizes, real_values, m_label, ymin, ymax, save_name=None):
    """
    Plot relative errors for the two variance components (gxg, e) of the
    epistasis-only model in stacked panels, fixed m and increasing n.

    Adapted from plot_three_effects_single_m for two-column result files whose
    rows are "(s2gxg, s2e)": column 0 = gxg, column 1 = e.

    Parameters
    ----------
    data_dict : dict
        Maps n -> DataFrame with col 0 = gxg, col 1 = e (as returned by
        read_MoM_results).
    individual_sizes : list of int
        Sample sizes plotted left-to-right (e.g. [1000, 2000, 4000, ...]).
    real_values : sequence of float
        True values [sigma^2_gxg, sigma^2_e]; relative error = estimate - truth.
    m_label : str
        Label shown at the top of the top panel (e.g. "m = 1000").
    ymin, ymax : float
        Shared y-axis limits for both panels.
    save_name : str, optional
        Filename without extension; PDF saved to the current working directory.
    """
    from scipy import stats

    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 9,
        'axes.linewidth': 1,
        'figure.dpi': 150,
        'savefig.dpi': 600,
    })

    components = [
        (0, r"$\sigma^2_{g \times g}$"),
        (1, r"$\sigma^2_{e}$"),
    ]

    num_boxes = len(individual_sizes)

    fig, axes = plt.subplots(
        2, 1,
        figsize=(3 + num_boxes * 0.7, 5.5),
        sharex=True,
        gridspec_kw={'hspace': 0.25}
    )

    box_color = '#3274A1'
    median_color = '#CC0000'
    label_color = '#000000'

    ymax_extended = ymax + 0.30 * (ymax - ymin)

    for ax_idx, (col_num, theta_label) in enumerate(components):
        ax = axes[ax_idx]
        real_value = real_values[col_num]

        box_data = []
        means = []
        sds = []
        p_values = []
        bottom_labels = []

        for n in individual_sizes:
            df = data_dict[n]
            col_values = df.iloc[:, col_num].values - real_value
            box_data.append(col_values)
            means.append(np.mean(col_values))
            sds.append(np.std(col_values, ddof=1))
            _, p_val = stats.ttest_1samp(col_values, 0)
            p_values.append(p_val)
            bottom_labels.append(f"n = {n:,}")

        x_positions = np.arange(len(individual_sizes))

        ax.set_ylim(ymin, ymax_extended)
        ax.set_xlim(-0.6, len(individual_sizes) - 0.4)

        # Alternating background
        for i in range(len(individual_sizes)):
            color = 'white' if i % 2 == 0 else '#E8E8E8'
            ax.axvspan(i - 0.5, i + 0.5, facecolor=color, alpha=1.0 if i % 2 == 0 else 0.8, zorder=0)

        # Box plots
        ax.boxplot(
            box_data,
            positions=x_positions,
            widths=0.5,
            patch_artist=True,
            showfliers=False,
            boxprops=dict(linewidth=1.2, edgecolor=box_color, facecolor='white'),
            whiskerprops=dict(linewidth=1, color=box_color),
            capprops=dict(linewidth=1, color=box_color),
            medianprops=dict(linewidth=1.5, color=median_color)
        )

        # m = label at top of top panel only
        if ax_idx == 0:
            ax.text(
                (len(individual_sizes) - 1) / 2,
                ymax_extended - 0.01 * (ymax_extended - ymin),
                m_label,
                ha='center', va='top', fontsize=10, fontweight='bold', color=label_color
            )

        top_offset = 0.09 if ax_idx == 0 else 0.03

        for i, (mean_val, sd_val, p_val) in enumerate(zip(means, sds, p_values)):
            if p_val < 0.001:
                sig_stars = '***'
            elif p_val < 0.01:
                sig_stars = '**'
            elif p_val < 0.05:
                sig_stars = '*'
            else:
                sig_stars = 'ns'

            ax.text(
                i, ymax_extended - top_offset * (ymax_extended - ymin),
                sig_stars,
                ha='center', va='top', fontsize=9, fontweight='bold', color=label_color
            )
            ax.text(
                i, ymax_extended - (top_offset + 0.08) * (ymax_extended - ymin),
                f"Mean={mean_val:.3f}",
                ha='center', va='top', fontsize=7.5, fontweight='bold', color=label_color
            )
            ax.text(
                i, ymax_extended - (top_offset + 0.16) * (ymax_extended - ymin),
                f"SD={sd_val:.3f}",
                ha='center', va='top', fontsize=7.5, fontweight='bold', color=label_color
            )

        # Reference line at zero
        ax.axhline(0, color='#666666', linestyle='--', linewidth=0.8, zorder=1)

        # Y-axis label
        ax.set_ylabel(f"Relative error ({theta_label})", fontsize=9, labelpad=6)

        # Y-axis ticks
        yticks = np.arange(np.ceil(ymin * 2) / 2, ymax + 0.01, 0.5)
        yticks = [y for y in yticks if y <= ymax]
        ax.set_yticks(yticks)

        # X-axis (only bottom)
        ax.set_xticks(x_positions)
        if ax_idx == len(components) - 1:
            ax.set_xticklabels(bottom_labels, fontsize=8)
            ax.set_xlabel("Sample size (n)", fontsize=9, labelpad=6)
        else:
            ax.set_xticklabels([])

        # Spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#333333')
        ax.spines['bottom'].set_color('#333333')
        ax.tick_params(axis='both', which='major', labelsize=8, colors='#333333')

    if save_name:
        script_dir = os.getcwd()
        save_path = os.path.join(script_dir, f"{save_name}.pdf")
        plt.savefig(save_path, bbox_inches='tight', facecolor='white',
                    edgecolor='none', format='pdf', transparent=False, pad_inches=0.1)
        print(f"PDF saved to: {save_path}")

    plt.show()
    plt.rcParams.update(plt.rcParamsDefault)


def plot_gxg_across_densities_xaxis(data_dicts, density_labels, individual_sizes, real_value, m_label, ymin, ymax, save_name=None):
    """
    Plot relative errors of gxg estimates with density on x-axis.
    Each panel = one sample size; within each panel, boxplots = different densities.
    
    Parameters:
    -----------
    data_dicts : list of dict
        List of dictionaries, one per density. Each dict maps n -> DataFrame with col 1 = gxg.
    density_labels : list of str
        Labels for each density (e.g., ["1.0", "0.5", "0.1", "0.02"])
    individual_sizes : list of int
        Sample sizes (e.g., [1000, 2000, 4000, 8000, 16000, 32000]).
    real_value : float
        True sigma^2_gxg value.
    m_label : str
        Label for SNP count shown at top (e.g., "m = 1000").
    ymin, ymax : float
        Y-axis limits for each panel.
    save_name : str, optional
        Filename without extension.
    """
    from scipy import stats
    
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 9,
        'axes.linewidth': 1,
        'figure.dpi': 150,
        'savefig.dpi': 600,
    })
    
    n_panels = len(individual_sizes)  # One panel per sample size
    num_boxes = len(data_dicts)  # Number of density values per panel
    
    # Create stacked panels (one per sample size)
    fig, axes = plt.subplots(
        n_panels, 1,
        figsize=(3 + num_boxes * 0.7, 2.5 * n_panels),
        sharex=True
    )
    
    if n_panels == 1:
        axes = [axes]
    
    box_color = '#3274A1'
    median_color = '#CC0000'
    label_color = '#000000'
    
    ymax_extended = ymax + 0.30 * (ymax - ymin)
    
    for ax_idx, n in enumerate(individual_sizes):
        ax = axes[ax_idx]
        
        # Gather data for this sample size across densities
        box_data = []
        means = []
        sds = []
        p_values = []
        bottom_labels = []
        
        for data_dict, density_label in zip(data_dicts, density_labels):
            df = data_dict[n]
            col_values = df.iloc[:, 1].values - real_value  # column 1 is gxg
            box_data.append(col_values)
            means.append(np.mean(col_values))
            sds.append(np.std(col_values, ddof=1))
            _, p_val = stats.ttest_1samp(col_values, 0)
            p_values.append(p_val)
            bottom_labels.append(f"density = {density_label}")
        
        x_positions = np.arange(num_boxes)
        
        ax.set_ylim(ymin, ymax_extended)
        ax.set_xlim(-0.6, num_boxes - 0.4)
        
        # Alternating background
        for i in range(num_boxes):
            color = 'white' if i % 2 == 0 else '#E8E8E8'
            ax.axvspan(i - 0.5, i + 0.5, facecolor=color, alpha=1.0 if i % 2 == 0 else 0.8, zorder=0)
        
        # Box plots
        ax.boxplot(
            box_data,
            positions=x_positions,
            widths=0.5,
            patch_artist=True,
            showfliers=False,
            boxprops=dict(linewidth=1.2, edgecolor=box_color, facecolor='white'),
            whiskerprops=dict(linewidth=1, color=box_color),
            capprops=dict(linewidth=1, color=box_color),
            medianprops=dict(linewidth=1.5, color=median_color)
        )
        
        # m = label at top of top panel only
        if ax_idx == 0:
            ax.text(
                (num_boxes - 1) / 2,
                ymax_extended - 0.01 * (ymax_extended - ymin),
                m_label,
                ha='center', va='top', fontsize=10, fontweight='bold', color=label_color
            )
        
        # Significance, Mean, SD labels at top
        top_offset = 0.09 if ax_idx == 0 else 0.03
        for i, (mean_val, sd_val, p_val) in enumerate(zip(means, sds, p_values)):
            if p_val < 0.001:
                sig_stars = '***'
            elif p_val < 0.01:
                sig_stars = '**'
            elif p_val < 0.05:
                sig_stars = '*'
            else:
                sig_stars = 'ns'
            
            ax.text(i, ymax_extended - top_offset * (ymax_extended - ymin),
                    sig_stars, ha='center', va='top', fontsize=9,
                    fontweight='bold', color=label_color)
            ax.text(i, ymax_extended - (top_offset + 0.07) * (ymax_extended - ymin),
                    f"Mean={mean_val:.3f}", ha='center', va='top',
                    fontsize=7.5, fontweight='bold', color=label_color)
            ax.text(i, ymax_extended - (top_offset + 0.14) * (ymax_extended - ymin),
                    f"SD={sd_val:.3f}", ha='center', va='top',
                    fontsize=7.5, fontweight='bold', color=label_color)
        
        # Reference line at zero
        ax.axhline(0, color='#666666', linestyle='--', linewidth=0.8, zorder=1)
        
        # Y-axis label uses sample size
        ax.set_ylabel(f"Relative error ($\\sigma^2_{{g \\times g}}$, n = {n:,})", 
                     fontsize=9, labelpad=6)
        
        # Y-axis ticks
        yticks = np.arange(np.ceil(ymin * 2) / 2, ymax + 0.01, 0.5)
        yticks = [y for y in yticks if y <= ymax]
        ax.set_yticks(yticks)
        
        # X-axis (only bottom)
        ax.set_xticks(x_positions)
        if ax_idx == n_panels - 1:
            ax.set_xticklabels(bottom_labels, fontsize=9)
            ax.set_xlabel("Density", fontsize=10, labelpad=6)
        else:
            ax.set_xticklabels([])
        
        # Spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#333333')
        ax.spines['bottom'].set_color('#333333')
        ax.tick_params(axis='both', which='major', labelsize=8, colors='#333333')
    
    plt.subplots_adjust(left=0.12, right=0.97, top=0.96, bottom=0.08, hspace=0.25)
    
    if save_name:
        script_dir = os.getcwd()
        save_path = os.path.join(script_dir, f"{save_name}.pdf")
        plt.savefig(save_path, bbox_inches='tight', facecolor='white',
                    edgecolor='none', format='pdf', transparent=False, pad_inches=0.1)
        print(f"PDF saved to: {save_path}")
    
    plt.show()
    plt.rcParams.update(plt.rcParamsDefault)


def plot_four_effects_single_m(panels, individual_sizes, m_label,
                               save_name=None, shared_ylim=None,
                               boundary_note=None):
    """Four stacked box panels, one per variance component, in the house style.

    The four-component sibling of plot_three_effects_single_m.  Each panel is
    one component's deviation from its reference, drawn exactly like the
    single-panel figures this family produces: one box per sample size,
    alternating background bands, a dashed line at zero, and a one-sample
    t-test against zero reported as ns / * / ** / *** above each box with its
    mean and SD.

    panels
        [(y_axis_label, {n: 1-D array of DEVIATIONS}), ...] in draw order.  The
        caller centres the data, not this function: the epistasis reference is
        c * s2gxg and therefore differs at every n, and with --paired every
        reference is per-replicate, neither of which a single `real_value`
        could express.
    shared_ylim
        (ymin, ymax) to force one scale on all four panels.  Left None, each
        panel is scaled to its own whiskers -- the components' errors differ by
        a factor of three at the same n, and a common scale spends most of the
        environment panel on white space.
    boundary_note
        Optional string printed under the x-axis, for reporting how many
        replicates sit at mc_reml's lower clamp.
    """
    from scipy import stats

    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 9,
        'axes.linewidth': 1,
        'figure.dpi': 150,
        'savefig.dpi': 600,
    })

    num_boxes = len(individual_sizes)
    fig, axes = plt.subplots(
        len(panels), 1,
        figsize=(3 + num_boxes * 0.7, 2.45 * len(panels)),
        sharex=True,
        gridspec_kw={'hspace': 0.22}
    )

    box_color = '#3274A1'
    median_color = '#CC0000'
    label_color = '#000000'

    x_positions = np.arange(num_boxes)
    bottom_labels = [f"n = {n:,}" for n in individual_sizes]

    for ax_idx, (y_label, dev) in enumerate(panels):
        ax = axes[ax_idx]

        box_data, means, sds, p_values = [], [], [], []
        for n in individual_sizes:
            v = np.asarray(dev[n], dtype=float)
            box_data.append(v)
            means.append(v.mean())
            sds.append(v.std(ddof=1))
            p_values.append(stats.ttest_1samp(v, 0)[1])

        # Panel scale: the whiskers, not the outliers -- showfliers is off, so
        # letting the extremes set the limits would shrink the boxes to slivers.
        if shared_ylim is not None:
            lo, hi = shared_ylim
        else:
            lo, hi = np.inf, -np.inf
            for v in box_data:
                q1, q3 = np.percentile(v, [25, 75])
                w = 1.5 * (q3 - q1)
                lo = min(lo, v[v >= q1 - w].min())
                hi = max(hi, v[v <= q3 + w].max())
            span = (hi - lo) or 1.0
            lo, hi = lo - 0.08 * span, hi + 0.08 * span

        # Headroom for the three annotation rows, plus the m label on panel 0.
        head = 0.50 if ax_idx == 0 else 0.40
        hi_ext = hi + head * (hi - lo)
        ax.set_ylim(lo, hi_ext)
        ax.set_xlim(-0.6, num_boxes - 0.4)

        for i in range(num_boxes):
            colour = 'white' if i % 2 == 0 else '#E8E8E8'
            ax.axvspan(i - 0.5, i + 0.5, facecolor=colour,
                       alpha=1.0 if i % 2 == 0 else 0.8, zorder=0)

        ax.boxplot(
            box_data,
            positions=x_positions,
            widths=0.5,
            patch_artist=True,
            showfliers=False,
            boxprops=dict(linewidth=1.2, edgecolor=box_color, facecolor='white'),
            whiskerprops=dict(linewidth=1, color=box_color),
            capprops=dict(linewidth=1, color=box_color),
            medianprops=dict(linewidth=1.5, color=median_color)
        )

        rng = hi_ext - lo
        if ax_idx == 0:
            ax.text((num_boxes - 1) / 2, hi_ext - 0.01 * rng, m_label,
                    ha='center', va='top', fontsize=10, fontweight='bold',
                    color=label_color)

        top_offset = 0.11 if ax_idx == 0 else 0.03
        for i, (mean_val, sd_val, p_val) in enumerate(zip(means, sds, p_values)):
            if p_val < 0.001:
                stars = '***'
            elif p_val < 0.01:
                stars = '**'
            elif p_val < 0.05:
                stars = '*'
            else:
                stars = 'ns'
            ax.text(i, hi_ext - top_offset * rng, stars,
                    ha='center', va='top', fontsize=9, fontweight='bold',
                    color=label_color)
            ax.text(i, hi_ext - (top_offset + 0.09) * rng, f"Mean={mean_val:.3f}",
                    ha='center', va='top', fontsize=7.5, fontweight='bold',
                    color=label_color)
            ax.text(i, hi_ext - (top_offset + 0.18) * rng, f"SD={sd_val:.3f}",
                    ha='center', va='top', fontsize=7.5, fontweight='bold',
                    color=label_color)

        ax.axhline(0, color='#666666', linestyle='--', linewidth=0.8, zorder=1)
        ax.set_ylabel(y_label, fontsize=9, labelpad=6)

        step = _nice_tick_step(hi - lo, target=5)
        ticks = [t for t in np.arange(np.ceil(lo / step) * step,
                                      hi + 0.5 * step, step) if lo <= t <= hi]
        ax.set_yticks(ticks)
        # One fixed precision per panel: matplotlib's default formatter prints
        # 0.1 next to 0.00 on the same axis, which reads as two scales.
        dec = next((d for d in range(7)
                    if all(abs(round(t, d) - t) < 1e-9 for t in ticks)), 6)
        ax.set_yticklabels([f"{t:.{dec}f}" for t in ticks])

        ax.set_xticks(x_positions)
        if ax_idx == len(panels) - 1:
            ax.set_xticklabels(bottom_labels, fontsize=8)
            ax.set_xlabel("Sample size (n)", fontsize=9, labelpad=6)
        else:
            ax.set_xticklabels([])

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#333333')
        ax.spines['bottom'].set_color('#333333')
        ax.tick_params(axis='both', which='major', labelsize=8, colors='#333333')

    if boundary_note:
        fig.text(0.5, 0.055, boundary_note, ha='center', va='top',
                 fontsize=7.5, color='#666666')

    if save_name:
        script_dir = os.getcwd()
        save_path = os.path.join(script_dir, f"{save_name}.pdf")
        plt.savefig(save_path, bbox_inches='tight', facecolor='white',
                    edgecolor='none', format='pdf', transparent=False,
                    pad_inches=0.1)
        print(f"PDF saved to: {save_path}")

    plt.show()
    plt.rcParams.update(plt.rcParamsDefault)
    return fig
