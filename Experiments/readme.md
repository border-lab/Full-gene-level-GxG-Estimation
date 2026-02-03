# Experiments

This file contains Jupyter notebooks for validating and exploring the gene-level GxG (gene-gene interaction) estimation method using Method of Moments (MoM).

## Notebooks

### 1. simulate_MoM.ipynb

Validates the MoM estimation method through Monte Carlo simulations.

**Experiments:**
- **Experiment 1**: Performance validation using randomly selected SNPs
- **Experiment 2**: Investigation of estimation challenges with contiguous SNPs


### 2. solve_ContiguousSNP.ipynb

Investigates the challenges of GxG estimation when using contiguous (neighboring) SNPs, particularly the effects of linkage disequilibrium (LD).

**Experiments:**
- **Experiment 1**: Validate the closed-form expression for Var$(Z_i Z_j)$ by comparing formula-based variance against empirical variance using simulated SNPs across increasing sample sizes

- **Experiment 2**: Compare the overall distribution of pairwise product variance Var$(Z_i Z_j)$ between contiguous and random SNPs to understand why the MoM estimator produces biased results with contiguous SNPs

- **Experiment 3**: Examine the relationship between pairwise product variance Var$(Z_i Z_j)$ and LD ($r_{ij}$) to identify high LD as the main factor behind inflated variance and biased estimation
- **Experiment 4**: Systematic investigation of three factors affecting $h^2_{g \times g}$ estimation — LD level (high/medium/low), MAF threshold (all variants vs. common variants), and SNP number (small vs. large dataset) — across 12 groups. Key findings:
  - Low LD outperforms medium LD, which outperforms high LD
  - MAF filtering shows no clear benefit and can increase bias in medium/high LD settings
  - SNP number does not affect bias, but its effect on standard error depends on LD level
---

### 3. simulation_weight.ipynb

Explores weight optimization strategies to reduce variance in GxG estimates.

To be continue

---

##  Data and Code Dependencies

The notebooks use test data from `../test_data/`.
All Python modules are located in `../python_scripts/`.