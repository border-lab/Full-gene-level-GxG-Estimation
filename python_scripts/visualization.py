import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os


def plot_relative_error_accross_sample_size(*dfs, basic_individual, col_num, real_value, save_path=None):

    """Plot relative errors with error bars (sorted by numeric sample size)."""
    """Accept multiple dataframes corresponding to different sample sizes"""
    """basic_individual: the base individual size (e.g., 1000)"""
    """col_num: column index to plot (generally 0 for s2gxg, 1 for s2e)"""
    """real_value: the true value of the parameter being estimated"""

    # Gather data 
    sample_sizes = [basic_individual * (2 ** i) for i in range(len(dfs))]  
    data_list = []

    for N, df in zip(sample_sizes, dfs):
        col_values = df.iloc[:, col_num].values - real_value
        data_list.append(
            pd.DataFrame({"value": col_values, "N": N})
        )

    data = pd.concat(data_list, ignore_index=True)

    #  Compute summary
    summary = (
        data.groupby("N")["value"]
        .agg(["mean", "std", "count"])
        .sort_index()         
    )
    summary["se"] = summary["std"]

    # x-axis positions
    x_positions = np.arange(len(summary))

    plt.figure(figsize=(10, 5))

    #  Strip plot 
    sns.stripplot(
        x="N", y="value", data=data,
        color="black", size=3, jitter=True, alpha=0.6,
        order=summary.index  
    )

    # Error bars + mean dots
    for i, (N, row) in enumerate(summary.iterrows()):
        plt.errorbar(
            x=i, y=row["mean"], yerr=row["se"],
            fmt="none", ecolor="red", elinewidth=3,
            capsize=8, capthick=2.5, alpha=0.9, zorder=5
        )
        plt.plot(
            i, row["mean"], marker="o", color="blue",
            markersize=6, markeredgecolor="black",
            markeredgewidth=0.5, zorder=6
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
        f"Change of relative error of {theta} by Sample Size\n(real value = {real_value}) ",
        fontsize=14, pad=10
    )
    plt.xlabel("Sample Size (N)", fontsize=12)
    plt.ylabel(f"Relative error: ({theta} - {real_value}) / {real_value}", fontsize=12)

    # Correct tick labels
    plt.xticks(ticks=x_positions, labels=[f"N={N}" for N in summary.index], fontsize=11)

    plt.tight_layout()

    #  Save or show plot
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=600, bbox_inches="tight")
        print(f" Plot saved to: {save_path}")
    else:
        plt.show()


def plot_var_pairwise_products_against_ld(var, ld, q=0.999, save_path=None):
    
    """Plot variance of pairwise products against LD (r). """
    """var: 1D array of variances of pairwise products"""
    """ld: 1D array of LD (r) values corresponding to each variance"""
    """q: quantile threshold to remove outliers (e.g., 0.999 removes top 0.1% of variances)"""

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
        alpha=0.3,      
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
    """Plot the distribution of Var(ZiZj) within a specified range."""

    """varZiZj_array： 1D array of Var(ZiZj) values"""
    """start, end 1D: range for x-axis"""
    """num_bin: number of bins in the histogram"""
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
    plt.ylabel('Percentage (%)', fontsize=12) 
  
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.show()