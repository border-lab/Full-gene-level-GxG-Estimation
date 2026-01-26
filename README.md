# Full-gene-level-GxG-Estimation

A Stochastic Method of Moments (MoM) approach for estimating gene-level epistatic (GxG) variance components.

## Overview

This project implements statistical methods to estimate the contribution of epistatic (gene-gene) interactions to phenotypic variance. The approach uses a computationally efficient Method of Moments estimator that avoids explicit formation of large kernel matrices.

## Methodology

### Variance Decomposition

Phenotypic variance is decomposed into: $Var(y) = \sigma^2_{g \times g} \times W + \sigma^2_e \times I$

Where:
- **σ²_gxg**: GxG epistatic variance component (target parameter)
- **σ²_e**: Residual/environmental variance
- **W**: GxG interaction kernel matrix
- **I**: Identity matrix

### GxG Interaction Kernel

The interaction kernel W captures pairwise epistatic effects:
W = 0.5 × (K⊙K - D@Dᵀ) / p


Where K = ZZᵀ is the linear kernel, D = Z⊙Z (element-wise square), and p is the number of interaction pairs.

### Method of Moments Estimation

The estimator solves moment equations using quadratic forms with stochastic trace estimation for computational efficiency.

## Project Structure
This project has three part....

See the `readme.md` in each subdirectory for detailed documentation.


## Installation

```bash
# Clone the repository
git clone https://github.com/[username]/Full-gene-level-GxG-Estimation.git
cd Full-gene-level-GxG-Estimation

# Install dependencies
pip install numpy pandas matplotlib jupyter