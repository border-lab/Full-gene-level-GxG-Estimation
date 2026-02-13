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

    # Get baseline SD (first sample size)
    baseline_std = summary["std"].iloc[0]

    # x-axis positions
    x_positions = np.arange(len(summary))

    plt.figure(figsize=(10, 5))

    # Strip plot 
    sns.stripplot(
        x="N", y="value", data=data,
        color="black", size=3, jitter=True, alpha=0.6,
        order=summary.index  
    )

    # Error bars (using std) + mean dots + text annotations
    for i, (N, row) in enumerate(summary.iterrows()):
        plt.errorbar(
            x=i, y=row["mean"], yerr=row["std"],  # Using std for error bars
            fmt="none", ecolor="red", elinewidth=3,
            capsize=8, capthick=2.5, alpha=0.9, zorder=5
        )
        plt.plot(
            i, row["mean"], marker="o", color="blue",
            markersize=6, markeredgecolor="black",
            markeredgewidth=0.5, zorder=6
        )
        
        # Calculate fold change relative to baseline
        fold_change = row["std"] / baseline_std
        
        # Determine color: green if SD decreased (fold < 1), red if increased (fold > 1)
        if i == 0:
            # First sample size - no fold change, use black
            fold_text = ""
            box_color = 'white'
        else:
            fold_text = f"\n({fold_change:.2f}x)"
            box_color = 'lightgreen' if fold_change < 1 else 'lightcoral'
        
        # Add text annotation for mean and std with fold change
        plt.text(
            i, ymax - 0.05 * (ymax - ymin),  # Position near top of plot
            f"Mean: {row['mean']:.4f}\nSD: {row['std']:.4f}{fold_text}",
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
        f"Change of relative error of {theta} by Sample Size\n(real value = {real_value}, error bars = SD)",
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


def read_MoM_results(individual_sizes, path, num_snp):
    """
    read_MoM_results from cluster output files.

    parameters:
    individual_sizes: list of individual sizes (e.g., [1000, 2000, 4000])
    path: directory path where the result files are stored
    num_snp: number of SNPs used in the analysis

    returns:
    A dictionary with individual sizes as keys and corresponding DataFrames as values.
    """ 

    data_dict = {}
    for n in individual_sizes:
        file_name = f"{n}m{num_snp}.txt"
        file_path = os.path.join(path, file_name)
        
        df = pd.read_csv(file_path, header=None)
        df[0] = df[0].astype(str).str.replace('(', '', regex=False).astype(float)
        df[1] = df[1].astype(str).str.replace(')', '', regex=False).astype(float)
        
        data_dict[n] = df
        
    return data_dict