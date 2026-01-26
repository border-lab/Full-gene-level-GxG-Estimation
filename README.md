# Full-gene-level-GxG-Estimation

A Stochastic Method of Moments (MoM) approach for estimating gene-level epistatic (GxG) variance components.

## Overview

This project implements statistical methods to estimate the contribution of epistatic (gene-gene) interactions to phenotypic variance. The approach uses a computationally efficient Method of Moments estimator that avoids explicit formation of large kernel matrices.

## Methodology

### Variance Decomposition

Phenotypic variance is decomposed into: $Var(y) = \sigma^2_{g \times g} \times W + \sigma^2_e \times I$

Where:
- **$\sigma^2_{g \times g}$**: GxG epistatic variance component (target parameter)
- **$\sigma^2_e$**: Residual/environmental variance
- **$W$**: GxG interaction kernel matrix
- **$I$**: Identity matrix

### GxG Interaction Kernel

The interaction kernel W captures pairwise epistatic effects:
$W = 0.5 × (K \circ K - D \circ D^T) / p$


Where $K = ZZ^T$ is the linear kernel, $D = Z \circ Z $(element-wise square), and p is the number of interaction pairs. $\circ$ is the haardmard product.


## Project Structure

This project has three parts:

### python_scripts
Core Python modules implementing the MoM estimator and supporting utilities.

### test_data
Genomic datasets used for validation experiments.

### Experiments
Jupyter notebooks for validating the MoM estimator and perform other analysis.

See the `readme.md` in each subdirectory for detailed documentation.

## Installation

```bash
# Clone the repository
git clone https://github.com/ZiyanZhang7315/Full-gene-level-GxG-Estimation.git

```
## How to use 
```bash
# cd to the dir of Experiments
cd Full-gene-level-GxG-Estimation/Experiments
```

Go to the direction of Experiments, there will be jupyter notebook about how to run simulation and introduction