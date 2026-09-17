import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define sample sizes and SNP counts
n_sizes = [1000, 2000, 4000, 8000, 16000, 32000]
m_sizes = [1000, 2000, 4000]

# Read all files and compute means
time_data = {}

for m in m_sizes:
    for n in n_sizes:
        df = pd.read_csv(f"data/time/n{n}m{m}time.txt", header=None)
        time_data[(n, m)] = df[0].mean()

# Set style
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.linewidth': 1,
    'figure.dpi': 150,
})

# Define fitting functions
def linear(x, a, b):
    return a * x + b

def quadratic_simple(x, a):
    return a * x**2

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))


colors_m = ['#E07850', '#3274A1', '#3A923A']

for i, m in enumerate(m_sizes):
    means = [time_data[(n, m)] for n in n_sizes]
    n_array = np.array(n_sizes)
    
    ax1.scatter(n_array, means, marker='o', s=30, color=colors_m[i], zorder=3)
    
    popt, _ = curve_fit(linear, n_array, means)
    n_fit = np.linspace(min(n_sizes), max(n_sizes), 100)
    ax1.plot(n_fit, linear(n_fit, *popt), '--', color=colors_m[i], linewidth=1.5,
             label=f'm = {m}')

ax1.set_xlabel('Individual Size (n)', fontsize=11)
ax1.set_ylabel('Estimation Time (s)', fontsize=11)
ax1.legend(title='Number of SNPs', frameon=False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.set_title('(A) Linear scaling with n: O(n)', fontsize=12, fontweight='bold', loc='center')

colors_n = ['#E07850', '#3274A1', '#3A923A', '#C03D3E', '#9372B2', '#8C564B']

for i, n in enumerate(n_sizes):
    means = [time_data[(n, m)] for m in m_sizes]
    m_array = np.array(m_sizes)
    
    ax2.scatter(m_array, means, marker='o', s=30, color=colors_n[i], zorder=3)
    
    popt, _ = curve_fit(quadratic_simple, m_array, means)
    m_fit = np.linspace(min(m_sizes), max(m_sizes), 100)
    ax2.plot(m_fit, quadratic_simple(m_fit, *popt), '--', color=colors_n[i], linewidth=1.5,
             label=f'n = {n//1000}K')

ax2.set_xlabel('Number of SNPs (m)', fontsize=11)
ax2.set_ylabel('Estimation Time (s)', fontsize=11)
ax2.legend(title='Individual Size', frameon=False, loc='upper left', fontsize=8)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.set_title('(B) Quadratic scaling with m: O(m²)', fontsize=12, fontweight='bold', loc='center')

plt.tight_layout()
plt.savefig('time_scaling.pdf', bbox_inches='tight', format='pdf')
plt.show()

plt.rcParams.update(plt.rcParamsDefault)