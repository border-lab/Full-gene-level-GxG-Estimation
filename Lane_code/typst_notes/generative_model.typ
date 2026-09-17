= Pooled Model
\
Here is the simplified pooled epistatic model with no fixed effect:

$ y_k = sum_(g=1)^(G) sum^(m_g)_(i < j) Z_(k i) Z_(k j) gamma_(i j) + e_k $

$ gamma_(i j) ~ N(0, sigma^2_(g times g) / P), quad P = sum_(g=1)^G binom(m_g, 2) $

- $y_k$: the (mean-centred) phenotype of individual $k$.\
- $e_k ~ N(0, sigma_e^2)$: residual noise.

- $g = 1, dots, G$: the index of gene.

- $m_g$: the number of SNPs in gene $g$.

- $i < j$: a within-gene SNP pair; the inner sum runs over the $binom(m_g, 2)$ pairs of gene $g$, so cross-gene pairs are excluded.

- $Z in RR^(N times M)$: the standardized genotype matrix.

- $gamma_(i j)$: the epistatic effect of pair $(i, j)$, i.i.d. Gaussian with per-pair variance $sigma^2_(g times g) \/ P$.

- $sigma^2_(g times g)$: the total epistatic variance.

- $P$: the total number of within-gene pairs summed over all genes.


\
The same model can be written in covariance (variance-component) form, collecting
the per-pair effects into a single epistatic genetic value:

$ y = g_(g times g) + e $

$ g_(g times g) ~ N(0, sigma^2_(g times g) W), quad e ~ N(0, sigma_e^2 I) $

$ V = "Var"(y) = sigma^2_(g times g) W + sigma_e^2 I, quad W = 1/P  sum_(g=1)^(G) K_g $

- $g_(g times g) = H gamma$: the epistatic genetic value vector, i.e. the per-pair effects collected together; marginalizing $gamma$ gives $g_(g times g) ~ N(0, sigma^2_(g times g) W)$.

- $H in RR^(N times P)$: the interaction design matrix, stacked gene-by-gene as $H = [H_1 | dots | H_G]$; the column for pair $(i, j)$ is the element-wise product $Z_(dot.c i) circle.small Z_(dot.c j)$, and $H_g in RR^(N times P_g)$ holds the columns of gene $g$.

- $K_g = H_g H_g^T$: the within-gene kernel of gene $g$. Because the columns of $H$ split by gene, $H H^T = sum_(g=1)^G K_g$, so $W = 1/P sum_g K_g = 1/P H H^T$.

- $W$: the pooled epistatic relationship matrix.


This model assumes that each SNP pair belongs to one gene: the genes are
disjoint sets of SNPs and only within-gene (cis) pairs enter the sum, with no
between-gene (trans) interaction terms. This is precisely what makes the columns of
$H$ partition gene-by-gene, so the cross-gene blocks $H_g H_h^T$ ($g eq.not h$)
never appear and $H H^T = sum_(g=1)^G K_g$ holds.

- *infinitesimal architecture:* every within-gene pair is causal, and the effects share a single per-pair variance $sigma^2_(g times g) \/ P$ --- one pooled component $sigma^2_(g times g)$ rather than a separate variance per gene.

- *Independent genes:* between-gene LD is ignored, so the per-gene contributions add without double-counting and the pooled kernel $W = 1/P sum_g K_g$ is unbiased for the total cis epistatic variance.

= The Model in Notes
\
The notes use the same two-component structure, but normalize the kernel *per gene*:
the pooled kernel is the equal-weight average of $G$ self-normalized within-gene
kernels. It is written in covariance (variance-component) form,

$ y = g_(g times g) + e, quad
  g_(g times g) ~ N(0, sigma^2_(g times g) W), quad e ~ N(0, sigma_e^2 I) $

$ V = "Var"(y) = sigma^2_(g times g) W + sigma_e^2 I $

The construction is stated first in full generality --- genes of *arbitrary,
possibly unequal, size* --- and then specialized to the equal-size case that the
simulations use.

== General case: unequal-size genes
\
Let gene $g$ contain $m_g$ SNPs, with $sum_(g=1)^G m_g = m$; the gene sizes
$m_1, dots, m_G$ need not be equal. Each gene is self-normalized and the genes are
then pooled with equal weight:

$ K_g = 1/p_g sum_(a < b, thick a\,b in cal(G)_g) H_(a b) H_(a b)^T,
  quad p_g = binom(m_g, 2) , $

$ W = 1/G sum_(g=1)^G K_g . $

- $cal(G)_g$: the set of SNPs in gene $g$ --- a contiguous block of $m_g$ SNPs.

- $H_(a b) in RR^N$: the *column-standardized* within-gene product of SNPs $a, b$ (mean $0$, variance $1$ across individuals).

- $p_g = binom(m_g, 2)$: the number of within-gene pairs in gene $g$; it *varies across genes* whenever the $m_g$ differ.

- $K_g = 1/p_g sum_(a<b in cal(G)_g) H_(a b) H_(a b)^T$: the *self-normalized* within-gene kernel of gene $g$, with $"tr"(K_g) = N$ regardless of $m_g$.

- $W = 1/G sum_g K_g$: the pooled kernel with *equal weight per gene*, $"tr"(W) = N$; every gene contributes $sigma^2_(g times g) \/ G$ to the epistatic variance, *independent of its size* $m_g$.

The per-gene self-normalization ($1\/p_g$) followed by the uniform $1\/G$ average
is what makes the modelling assumption explicit: all genes share one epistatic
variance $sigma^2_(g times g)$, split *evenly* across genes, so a large gene and a
small gene each carry $sigma^2_(g times g) \/ G$ even though the large gene holds
many more pairs. This is a genuinely different weighting from the pooled model of
the previous section, which weights each *pair* equally and hence lets a larger
gene carry proportionally more heritability ($prop m_g^2$).

== Special case: equal-size genes
\
When every gene has the same size $m_g = m\/G$, all the pair counts coincide,

$ p_g = binom(m\/G, 2) equiv p quad "for all" g , $

so the self-normalization constant $1\/p_g$ is the same for every gene and the
per-gene average collapses onto the per-pair pooling:

$ W = 1/G sum_(g=1)^G K_g = 1/P H H^T, quad P = G p = G binom(m\/G, 2) . $

Equivalently, in effect-level form with i.i.d. per-pair effects,

$ y_k = sum_(g=1)^(G) sum_(a < b, thick a\,b in cal(G)_g) (H_(a b))_k thin gamma_(a b) + e_k,
  quad gamma_(a b) ~ N(0, sigma^2_(g times g) \/ P), quad P = G p . $

In this equal-size regime the notes' per-gene weighting and the pooled model's
per-pair weighting *coincide exactly*, because equal $m_g$ makes "equal per gene"
and "equal per pair" the same thing. Its distinctive consequence is the
effective-marker reduction: restricting the kernel to within-gene pairs lowers its
effective number of markers by $sqrt(G)$,

$ M_e^"within" = sqrt(G) thin M_e^"gene" = M_e^"genome" \/ sqrt(G) . $

*Relationship to the pooled model.* The two models coincide exactly when genes are
equal-sized and diverge otherwise. They differ only in the normalization
convention: the pooled model weights each *pair* equally ($W = 1/P sum_g K_g$ with
un-normalized $K_g = H_g H_g^T$), so a larger gene carries proportionally more
heritability ($prop m_g^2$); the notes weight each *gene* equally
($W = 1/G sum_g K_g$ with self-normalized $K_g$), giving every gene the same
$sigma^2_(g times g) \/ G$. For unequal gene sizes the two diverge; for equal gene
sizes the divergence vanishes.
