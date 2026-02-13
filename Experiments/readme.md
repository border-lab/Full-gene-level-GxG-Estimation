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
- **Experiment 4**: Systematic investigation of three factors affecting $h^2_{g \times g}$ estimation — LD level, MAF, and number of SNPs — across 12 groups (3 LD levels × 2 MAF thresholds × 2 SNP sizes):
  - **4.1 LD level**: Bias increases with LD level; high-LD shows greatest bias, low-LD is nearly unbiased. Greater independence among SNPs leads to more consistent estimators.
  - **4.2 MAF**: Common variants reduce standard error and improve convergence speed, but unexpectedly increase bias under high LD conditions.
  - **4.3 Number of SNPs**: SNP count does not affect bias but influences standard error — larger datasets provide lower SE and faster convergence, particularly under low LD.
  - **Summary table**: Relative error of heritability estimates across all 12 groups at N = 1K, 2K, 4K, 8K.

---

### 3. simulation_weight.ipynb

Explores weight optimization strategies to reduce variance in GxG estimates.

To be continue

---

##  Data and Code Dependencies

The notebooks use test data from `../test_data/`.
All Python modules are located in `../python_scripts/`.


##  Update detail
1. Sample the dataset based on one big dataset
2. Filter MAF first, then LD
3. the num_MC increase to 150
4. Write the SE(fold-change) and mean values in the figures