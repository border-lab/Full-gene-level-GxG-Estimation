# Scripts

This file contains the core Python modules for gene-level GxG (gene-gene interaction) estimation.

### gene_level_MoM.py

Core estimation module implementing the MoM.

| Function | Description |
|----------|-------------|
| `simulate_from_raw_only_W` | Simulate phenotypes from raw genotype data with specified GxG and e variance components |
| `MoM_only_M` | Method of Moments estimator for variance components |

### weight_function.py

Weight function for the estimation procedure.

| Function | Description |
|----------|-------------|
| `weight_vector` | Compute weight vectors from genotype matrix |
| `check_weighted_products` | Validate weighted pairwise products |

### pairwise_product_var_calculation.py

Variance calculation for pairwise SNP products.

| Function | Description |
|----------|-------------|
| `calculate_varXiXj` | Calculate variance of pairwise products XᵢXⱼ (raw data) |
| `calculate_varZiZj` | Calculate variance of pairwise products ZᵢZⱼ (standardized data) |

### visualization.py

visualization part.

| Function | Description |
|----------|-------------|
| `plot_relative_error_accross_sample_size` | Plot relative estimation error across different sample sizes |
| `plot_var_pairwise_products_against_ld` | Scatter plot of pairwise product variance against LD (r) values |

## How to use

```python
from Scripts.gene_level_MoM import *
from Scripts.weight_function import *
from Scripts.pairwise_product_var_calculation import *
from Scripts.visualization import *
