import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

def plot_relative_error_accross_sample_size(*dfs, basic_individual, col_num, real_value, ymin, ymax, save_path=None):
    """
    Plot relative errors with error bars (sorted by numeric sample size).
    
    Parameters:
    -----------
    dfs: list of DataFrames containing estimation results at different sample sizes
    basic_individual: base sample size to compute actual sample sizes
    col_num: column number in DataFrames to plot (0 for s2gxg, 1 for s2e)
    real_value: the true value of the parameter being estimated
    ymin, ymax: y-axis limits
    save_path: optional path to save the plot image

    Returns:
    --------
    A plot showing relative errors across different sample sizes with error bars.
    """
    # Gather data 
    sample_sizes = [basic_individual * (2 ** i) for i in range(len(dfs))]  
    data_list = []

    for N, df in zip(sample_sizes, dfs):
        col_values = df.iloc[:, col_num].values - real_value
        data_list.append(
            pd.DataFrame({"value": col_values, "N": N})
        )

    data = pd.concat(data_list, ignore_index=True)

    # Compute summary statistics
    summary = (
        data.groupby("N")["value"]
        .agg(["mean", "std", "count"])
        .sort_index()         
    )
    
    # Calculate SE = std / sqrt(count)
    summary["se"] = summary["std"] / np.sqrt(summary["count"])

    # Get baseline SE (first sample size)
    baseline_se = summary["se"].iloc[0]

    # x-axis positions
    x_positions = np.arange(len(summary))

    plt.figure(figsize=(10, 5))

    # Strip plot 
    sns.stripplot(
        x="N", y="value", data=data,
        color="black", size=3, jitter=True, alpha=0.6,
        order=summary.index  
    )

    # Error bars (using SE) + mean dots + text annotations
    for i, (N, row) in enumerate(summary.iterrows()):
        plt.errorbar(
            x=i, y=row["mean"], yerr=row["se"],  # Using SE for error bars
            fmt="none", ecolor="red", elinewidth=3,
            capsize=8, capthick=2.5, alpha=0.9, zorder=5
        )
        plt.plot(
            i, row["mean"], marker="o", color="blue",
            markersize=6, markeredgecolor="black",
            markeredgewidth=0.5, zorder=6
        )
        
        # Calculate fold change relative to baseline
        fold_change = row["se"] / baseline_se
        
        # Determine color: green if SE decreased (fold < 1), red if increased (fold > 1)
        if i == 0:
            # First sample size - no fold change, use black
            fold_text = ""
            box_color = 'white'
        else:
            fold_text = f"\n({fold_change:.2f}x)"
            box_color = 'lightgreen' if fold_change < 1 else 'lightcoral'
        
        # Add text annotation for mean and SE with fold change
        plt.text(
            i, ymax - 0.05 * (ymax - ymin),  # Position near top of plot
            f"Mean: {row['mean']:.4f}\nSE: {row['se']:.4f}{fold_text}",
            ha='center', va='top', fontsize=9,
            bbox=dict(boxstyle='round,pad=0.3', facecolor=box_color, edgecolor='gray', alpha=0.8)
        )

    # Theta label
    if col_num == 0:
        theta = r"$\sigma^2_{g \times g}$"
    elif col_num == 1:
        theta = r"$\sigma^2_{e}$"
    else:
        theta = "Parameter"

    # Labels and title
    plt.axhline(0, color="gray", linestyle="--", linewidth=1)
    plt.title(
        f"Change of relative error of {theta} by Sample Size\n(real value = {real_value}, error bars = SE)",
        fontsize=14, pad=10
    )
    plt.ylim(ymin, ymax)
    plt.xlabel("Sample Size (N)", fontsize=12)
    plt.ylabel(f"Relative error: ({theta} - {real_value}) / {real_value}", fontsize=12)

    # Correct tick labels
    plt.xticks(ticks=x_positions, labels=[f"N={N}" for N in summary.index], fontsize=11)

    plt.tight_layout()

    # Save or show plot
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=600, bbox_inches="tight")
        print(f"Plot saved to: {save_path}")
    else:
        plt.show()

def plot_var_pairwise_products_against_ld(var, ld, q=0.999, save_path=None):
    
    """
    Plot variance of pairwise products against LD (r).
    parameters:
    var: 1D array of variances of pairwise products
    ld: 1D array of LD (r) values corresponding to the variances
    q: quantile threshold to remove outliers (e.g., 0.999 removes top 0.1% of variances)
    save_path: optional path to save the plot image

    returns:
    A scatter plot of variance of pairwise products vs LD (r) with outliers removed
    """

    var = np.asarray(var)
    ld = np.asarray(ld)

    # remove top (1-q)% variance outliers
    cutoff = np.quantile(var, q)
    mask = var <= cutoff
     
    # PLOT
    plt.figure(figsize=(8, 5))
    plt.scatter(
        ld[mask],
        var[mask],
        s=0.5,           
        alpha=0.2,      
        edgecolors='none'
    )
    plt.xlabel("LD (r)", fontsize=13)
    plt.ylabel("variance of pairwise products", fontsize=13)
    plt.title("variance of pairwise products vs LD (r)")
    plt.grid(alpha=0.2)

    plt.xlim(-1, 1)
    plt.axvline(x=0, linestyle=':', color='gray', alpha=0.5, linewidth=1)
    plt.axhline(y=1, linestyle='--', linewidth=1.5, label='y = 1')
    plt.legend()
    plt.tight_layout()

    # SAVE 
    if save_path is not None:
        plt.savefig(save_path, dpi=300)
        print(f"Plot saved to: {save_path}")

    plt.show()



def plot_distribution_varzizj(varZiZj_array, start, end, num_bin):
    """
    Plot the distribution of Var(ZiZj) within a specified range.

    parameters:
    varZiZj_array: 1D array of Var(ZiZj) values
    start: start of the range to plot
    end: end of the range to plot
    num_bin: number of bins in the histogram"

    returns:
    A histogram plot showing the distribution of Var(ZiZj) within the specified range.
    """
    plt.figure(figsize=(10, 6))
    # Histogram
    plt.hist(varZiZj_array, 
             bins=num_bin, 
             range=(start, end), 
             color='blue', 
             edgecolor='white', 
             alpha=0.6)
  
    plt.xlim(start, end)
  
    plt.title(fr'Distribution of $Var(Z_i Z_j)$ in range [{start}, {end}]', fontsize=15)
    plt.xlabel(r'Variance of ($Var(Z_i Z_j)$)', fontsize=12)
    plt.ylabel('Count', fontsize=12) 
  
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.show()


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


def plot_relative_error_across_groups_combined(*data_dicts, x_labels, individual_sizes, col_num, real_value, ymin, ymax, x_axis_name="Group", save_path=None):
    """
    Plot relative errors with error bars across different groups and sample sizes combined.
    
    Parameters:
    -----------
    data_dicts: dictionaries with sample sizes as keys (e.g., highLD_results, middleLD_results, ...)
    x_labels: list of labels for groups (e.g., ["High LD", "Middle LD", "Low LD"])
    individual_sizes: list of sample sizes (e.g., [1000, 2000, 4000, 8000])
    col_num: column number in DataFrames to plot
    real_value: the true value of the parameter being estimated
    ymin, ymax: y-axis limits
    x_axis_name: name for the x-axis
    save_path: optional path to save the plot image
    """
    # Gather data for all combinations
    data_list = []
    combined_labels = []
    
    for label, data_dict in zip(x_labels, data_dicts):
        for n in individual_sizes:
            df = data_dict[n]
            col_values = df.iloc[:, col_num].values - real_value
            combined_label = f"{label}\nN={n}"
            combined_labels.append(combined_label)
            data_list.append(pd.DataFrame({"value": col_values, "group": combined_label}))
    
    data = pd.concat(data_list, ignore_index=True)
    data["group"] = pd.Categorical(data["group"], categories=combined_labels, ordered=True)
    
    # Compute summary statistics
    summary = (
        data.groupby("group", observed=True)["value"]
        .agg(["mean", "std", "count"])
        .loc[combined_labels]
    )
    summary["se"] = summary["std"] / np.sqrt(summary["count"])
    
    # Calculate baseline SE for each LD group (first sample size within each group)
    baseline_ses = {}
    for i, label in enumerate(x_labels):
        first_combined_label = f"{label}\nN={individual_sizes[0]}"
        baseline_ses[label] = summary.loc[first_combined_label, "se"]
    
    # x-axis positions
    x_positions = np.arange(len(summary))
    
    # Create figure with appropriate width
    fig_width = max(12, len(combined_labels) * 1.2)
    plt.figure(figsize=(fig_width, 6))
    
    # Define colors for each LD level
    colors = plt.cm.tab10(np.linspace(0, 1, len(x_labels)))
    color_map = {label: colors[i] for i, label in enumerate(x_labels)}
    
    # Strip plot with colors by LD level (smaller points)
    for i, (label, data_dict) in enumerate(zip(x_labels, data_dicts)):
        for j, n in enumerate(individual_sizes):
            combined_label = f"{label}\nN={n}"
            subset = data[data["group"] == combined_label]
            x_pos = combined_labels.index(combined_label)
            
            plt.scatter(
                x=np.random.normal(x_pos, 0.1, len(subset)),
                y=subset["value"],
                color=color_map[label], s=5, alpha=0.4  # Smaller points: s=5
            )
    
    # Error bars + mean dots + text annotations
    for i, (combined_label, row) in enumerate(summary.iterrows()):
        plt.errorbar(
            x=i, y=row["mean"], yerr=row["se"],
            fmt="none", ecolor="red", elinewidth=2,
            capsize=5, capthick=2, alpha=0.9, zorder=5
        )
        plt.plot(
            i, row["mean"], marker="o", color="blue",
            markersize=5, markeredgecolor="black",
            markeredgewidth=0.5, zorder=6
        )
        
        # Determine which LD group this belongs to
        current_ld_label = combined_label.split("\n")[0]
        baseline_se = baseline_ses[current_ld_label]
        
        # Calculate fold change relative to baseline within the same LD group
        fold_change = row["se"] / baseline_se
        
        # Check if this is the first sample size in its LD group
        is_first_in_group = combined_label.endswith(f"N={individual_sizes[0]}")
        
        if is_first_in_group:
            fold_text = ""
            box_color = 'white'
        else:
            fold_text = f"\n({fold_change:.2f}x)"
            box_color = 'lightgreen' if fold_change < 1 else 'lightcoral'
        
        # Add text annotation
        plt.text(
            i, ymax - 0.03 * (ymax - ymin),
            f"Mean:{row['mean']:.3f}\nSE:{row['se']:.4f}{fold_text}",
            ha='center', va='top', fontsize=7,
            bbox=dict(boxstyle='round,pad=0.2', facecolor=box_color, edgecolor='gray', alpha=0.8)
        )
    
    # Add vertical separators between LD groups
    for i in range(1, len(x_labels)):
        sep_pos = i * len(individual_sizes) - 0.5
        plt.axvline(sep_pos, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    
    # Theta label
    if col_num == 0:
        theta = r"$\sigma^2_{g \times g}$"
    elif col_num == 1:
        theta = r"$\sigma^2_{e}$"
    else:
        theta = "Parameter"
    
    # Labels and title
    plt.axhline(0, color="gray", linestyle="--", linewidth=1)
    plt.title(
        f"Change of relative error of {theta} by {x_axis_name} and Sample Size\n(real value = {real_value}, error bars = SE)",
        fontsize=14, pad=10
    )
    plt.ylim(ymin, ymax)
    plt.xlabel(f"{x_axis_name} / Sample Size", fontsize=12)
    plt.ylabel(f"Relative error: ({theta} - {real_value})", fontsize=12)
    
    plt.xticks(ticks=x_positions, labels=combined_labels, fontsize=9, rotation=0)
    
    # Add legend for LD levels
    legend_handles = [plt.Line2D([0], [0], marker='o', color='w', 
                                  markerfacecolor=color_map[label], markersize=8, label=label)
                      for label in x_labels]
    plt.legend(handles=legend_handles, loc='upper right', fontsize=10)
    
    plt.tight_layout()
    
    # Save or show plot
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=600, bbox_inches="tight")
        print(f"Plot saved to: {save_path}")
    else:
        plt.show()