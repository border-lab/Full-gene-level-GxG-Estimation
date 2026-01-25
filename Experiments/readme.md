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
- **Experiment 1**: Analysis of pairwise product variance relationships
- **Experiment 2**: Distribution analysis of Var(ZiZj) for standardized genotypes
- **Experiment 3**: Correlation between variance and LD strength
- **Experiment 4**: Evaluation of LD-pruning strategies
- **Experiment 5**: Impact of removing highly correlated SNP pairs

---

### 3. simulation_weight.ipynb

Explores weight optimization strategies to reduce variance in GxG estimates.

---

##  Data and Code Dependencies

The notebooks use test data from `../test_data/`.
All Python modules are located in `../python_scripts/`.

