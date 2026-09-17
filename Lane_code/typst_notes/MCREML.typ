#set page(margin: 1in)
#set text(size: 11pt)
#set par(justify: true)
#set heading(numbering: "1.1")

#align(center)[
  #text(size: 16pt, weight: "bold")[
    Monte-Carlo AI-REML for a Pairwise-Epistasis Model
  ]
  #v(4pt)
  #text(size: 10pt)[Simulation and estimation pipelines: precomputed-$W$ (`Simulation_code_MCREML_gxg`) and matrix-free (`Simulation_code_MCREML_gxg_matfree`)]
]

= Identifiability <sec-identifiability>

- With only two components, the estimator must separate the epistasis kernel $W$ from the residual identity $I$. These are weakly identified: interactions of column-standardized SNPs are close to independent, so $W$ is often nearly collinear with $I$ ($W approx I$ in the bulk of its spectrum).

- As a consequence the average-information matrix $cal(A)$ can be near-singular and an undamped Newton step is unstable. The update is stabilized with a Levenberg--Marquardt ridge, a trust region on the step, and a box clamp $sigma_i^2 gt.eq 0$, which leaves well-identified cases unperturbed.


= Model

With $n$ individuals, $m$ SNPs, and no fixed effects ($y$ mean-centred):

$ y = g_"gxg" + e, $
$ V = "Var"(y) = sigma_"gxg"^2 W + sigma_e^2 I, $

with parameter vector $theta = (sigma_"gxg"^2, sigma_e^2)$
so $k = 2$.

= Relationship matrix

Let $Z$ be the column-standardized dosage matrix (each SNP column mean 0,
variance 1), and write $Z_a$ for its $a$-th column.

*Pairwise epistasis.* For each SNP pair $(a, b)$, $a < b$, form the element-wise
product $Z_a circle.small Z_b$ and standardize it to $H_(a b)$ (mean 0, variance
1). With $p = binom(m, 2) = m(m-1)\/2$ pairs,

$ W = 1/p sum_(a<b) H_(a b) H_(a b)^T . $


= AI-REML update

Every
$V^(-1)$ and every trace is formed directly from a dense factorization, with no
conjugate gradient and no stochastic probes.

Writing $x = V^(-1) y$ and $K_1 = W$, $K_2 = I$, each
iteration forms $V = sum_i sigma_i^2 K_i$ explicitly, factorizes the $n times n$
matrix (Cholesky / direct inverse, $O(n^3)$), and takes the average-information
Newton step $theta arrow.l theta + (cal(A) + lambda I)^(-1) s$, clamped to
$sigma_i^2 gt.eq 0$, with

$ s_i = 1/2 (x^T K_i x - "tr"(V^(-1) K_i)), quad
  cal(A)_(i j) = 1/2 (K_i x)^T V^(-1) (K_j x) . $

The difference from the MC update is purely computational: here every piece is
exact. $x = V^(-1) y$ and $V^(-1)(K_j x)$ are direct solves against the
factorization, and each $"tr"(V^(-1) K_i)$ is the full contraction
$sum_(r s) (V^(-1))_(r s) (K_i)_(r s)$ of two dense $n times n$ matrices
($O(n^2)$) -- so the estimate is deterministic given $y$, with no Hutchinson
noise and no CG tolerance.

*Cost.* Per iteration, assembling $V$ from the prebuilt $K_i$ is $O(k n^2)$, its
factorization $O(n^3)$, and each trace / $cal(A)_(i j)$ a further $O(n^2)$
contraction; the $O(n^3)$ factorization dominates. Exact AI-REML is thus
$O("iters" dot n^3)$ time and $O(n^2)$ memory, holding the dense
$W, V, V^(-1)$ (with $W$ built in $O(n^2 p)$ once). This is exactly the $O(n^3)$ /
$O(n^2)$ scaling that the matrix-free Monte-Carlo update below removes, at the
price of CG iterations and stochastic traces.

= MC AI-REML update

Writing $x = V^(-1) y$ and $K_1 = W$, $K_2 = I$, each
iteration takes the average-information Newton step
  $ theta arrow.l theta + (cal(A) + lambda I)^(-1) s $ with

$ s_i = 1/2 (x^T K_i x - "tr"(V^(-1) K_i)), quad
  cal(A)_(i j) = 1/2 (K_i x)^T V^(-1) (K_j x) . $

Data quadratics are exact given $x$:
$x^T W x = x^T (W x)$, $x^T I x = norm(x)^2$; the traces use $"Nmc"$ Rademacher
probes $U = [u_1, ..., u_"Nmc"] in RR^(n times "Nmc")$,
$"tr"(V^(-1) K_i) approx "Nmc"^(-1) sum_(b=1)^"Nmc" (V^(-1) u_b)^T K_i u_b$. Every
$V^(-1)$ is applied by batched CG whose only primitive is the $V$ mat-vec

$ V U = sigma_"gxg"^2 (W U) + sigma_e^2 U . $


$W U$ is the standardized linear operator used in MoM.
= Cost


Let $n$ individuals, $m$ SNPs, $"Nmc"$ probes, $k = 2$ components, $"iters"$ REML
iterations, and $t_"cg"$ CG iterations per solve. `mc_reml` is fully matrix-free:
every epistasis apply $W v$ goes through `compute_WU` at $c_W = O(n m^2)$ -- a
dense $W$ is never formed or passed in. The single CG primitive is the $V$
mat-vec on $c$ columns,

$ V U : quad O((c_W + n) thin c) quad
  (W " " c_W thin c ", residual " O(n c)) . $

*Setup* (once per replicate):
+ Design $Z$ (column-standardize) -- $O(n m)$.
+ Weight matrices $S, R, T$ (`compute_weight_matrices`) -- $O(n m^2)$, stored $m times m$.
+ Draw Rademacher probes $U in {plus.minus 1}^(n times "Nmc")$ -- $O(n "Nmc")$.
+ Fixed-probe product $W U$ -- $c_W "Nmc"$; built once, reused every iteration.

*Per REML iteration* ($times "iters"$):
+ Solve $x = V^(-1) y$, 1 RHS -- $O(t_"cg" (c_W + n))$.
+ Data quadratics $x^T K_i x$: $W x$ at $c_W$, $norm(x)^2$ at $O(n)$.
+ Trace probes $P = V^(-1) U$, $"Nmc"$ RHS -- $O(t_"cg" (c_W + n) "Nmc")$; then $"tr"(V^(-1) K_i) approx "Nmc"^(-1) sum_b P_b^T (K_i U)_b$ at $O(n "Nmc")$.
+ AI matrix: form $K_i x$ ($n times k$) ($W x$ reused, $I x = x$); solve $G = V^(-1)(K_i x)$, $k$ RHS -- $O(t_"cg" (c_W + n) k)$; then $cal(A) = 1/2 (K x)^T G$ at $O(n k^2)$.
+ AI-Newton step: solve the $k times k$ system -- $O(k^3)$, negligible.

The three solve groups carry $1 + "Nmc" + k = "Nmc" + 3$ RHS columns, so one
iteration costs $O(t_"cg" ("Nmc" + 3)(c_W + n))$ and the whole estimation is

$ O("iters" dot t_"cg" dot ("Nmc" + 3) dot (c_W + n))
  = O("iters" dot t_"cg" dot ("Nmc" + 3) dot n m^2) , $

to leading order, since the epistasis apply $c_W = O(n m^2)$ dominates the
residual term $O(n)$.

*A precomputed-$W$ speed-up is possible but not implemented here.* `compute_WU`
recomputes the $m$-by-$m$ contraction on every CG iteration, so at $m tilde n$ it
is $tilde m^2 \/ n$ times slower than applying a prebuilt dense $W$ would be
($O(n^2)$ per apply instead of $O(n m^2)$; $~2500 times$ at $n = m = 1000$). In
the current code `mc_reml` / `_v_matvec` take no $W$ argument and always rebuild
$(S, R, T)$ from $Z_a$, so estimation stays matrix-free. For $n gt.tilde m$ one
could form $W$ once (it is deterministic per genotype) and replace the
`compute_WU` in `_v_matvec` with a dense $W U$; the matrix-free apply remains
preferable in the large-$m$, memory-bound regime.


= Simulation


== Relative error and SE for random SNP for pre-computed and mat-free

We compare the two estimators of the epistasis variance -- precomputed-$W$ and
matrix-free -- on the same random-SNP phenotypes. @fig-matfree and @fig-prew show
box plots of the relative error (estimate minus truth) with $m = 1000$ held fixed
and $n$ increasing left-to-right. As $W u$ for precomputed-$W$ = $W u$ for linear operator in matrix-free method. These two polt should be similar.



#figure(
  image("RandomSNP_matfree_s2gxg0.2_s2e0.8_fixed_m1000_gxg.pdf", width: 100%),
  caption: [Matrix-free MC-AI-REML on random SNPs: relative error of $sigma_"gxg"^2$, fixed $m = 1000$, increasing $n$.],
) <fig-matfree>

#figure(
  image("RandomSNP_preW_s2gxg0.2_s2e0.8_fixed_m1000_gxg.pdf", width: 100%),
  caption: [Precomputed-$W$ MC-AI-REML on random SNPs: relative error of $sigma_"gxg"^2$, fixed $m = 1000$, increasing $n$.],
) <fig-prew>


- *Near-unbiased and consistent.* The mean relative error is small at every $n$ and
  converges to $0$, while the SE falls from $tilde 0.02$ at $n = 1000$ to $tilde 0.002$
  at $n = 16000$ -- the $tilde 1\/sqrt(n)$ decay expected of a consistent estimator.

- *Small-$n$ inflation.* At $n = 1000$, $sigma_"gxg"^2$ carries a mild upward bias
  ($+0.05$ to $+0.09$) with a wide spread ($"SD" approx 0.33$). This is the
  $W approx I$ weak identifiability of @sec-identifiability: with few individuals the
  epistasis and residual kernels are hard to separate, inflating the variance and
  pulling $sigma_"gxg"^2$ up (and $sigma_e^2$ down). It resolves as $n$ grows.

- *The two pipelines agree.* Precomputed-$W$ and matrix-free match within Monte-Carlo
  error at every $n$ (identical to the reported precision for $n gt.eq 2000$),
  confirming that applying $W$ matrix-free reproduces the pre-computed dense-$W$
  estimate: the pipelines differ only in how $W b$ is computed, not in the statistical
  behavior of the estimate.

== Running time check

We time the estimation step alone -- the `MC_REML` call, averaged over the
10 replicates (  ` --nodes=1, --ntasks=1  NMC=100 at mzhang partition`)  -- with $m = 1000$ fixed and $n$ increasing.

#figure(
  image("time_comparison.pdf", width: 100%),
  caption: [Average estimation (MC-AI-REML) wall-clock time versus sample size $n$ at fixed $m = 1000$, for the precomputed-$W$ and matrix-free $W$-apply strategies.],
) <fig-time>

- *Scaling matches the per-apply cost.* The fitted exponents are $tilde n^(1.9)$ for
  precomputed-$W$ and $tilde n^(1.0)$ for matrix-free, recovering the $O(n^2)$ dense
  $W b$ and the (linear-in-$n$) $O(n m^2)$ matrix-free apply from the Cost section --
  the probe / REML / CG iteration counts do not grow with $n$, so each apply's cost sets
  the slope.

- *the gap closes as $n$ grows.* Being quadratic against the matrix-free linear
  cost, the precomputed-$W$ speed-up shrinks from $247 times$ ($n = 1000$) to $29 times$
  ($n = 16000$); extrapolating the fits, the two cross near $n tilde 5 times 10^5$,
  beyond which the matrix-free apply is faster -- and it already stores only the
  $O(n m)$ genotype rather than the $O(n^2)$ dense $W$. Precomputed-$W$ is thus preferable
  in the small-$n$ / large-$m$ regime, matrix-free in the large-$n$, memory-bound regime.



= Pooled kernel

The pooled kernel encodes exactly this:
split the SNPs into genes, build one epistasis GRM per gene, and average them into
a single kernel that still fits with the *same* two-component model
$V = sigma_"gxg"^2 W + sigma_e^2 I$.

== Construction

Partition the $m$ SNP columns into $G$ contiguous, equal-size gene blocks
$cal(G)_1, ..., cal(G)_G$ of $m\/G$ SNPs each (e.g. $G = 10$, $m = 1000$ gives ten
100-SNP genes). For gene $g$, form its within-gene epistasis GRM from *only* its
own pairs, exactly as in the genome-wide case but restricted to $cal(G)_g$:

$ K_g = 1/p_g sum_(a < b, thick a\,b in cal(G)_g) H_(a b) H_(a b)^T,
  quad p_g = binom(m\/G, 2) , $

where $H_(a b)$ is again the standardized element-wise product of columns $a$ and
$b$. Cross-gene pairs $(a in cal(G)_g, b in cal(G)_(g'))$, $g eq.not g'$, are
*excluded*. The pooled kernel is the uniform average of the per-gene GRMs,

$ W = 1/G sum_(g=1)^G K_g . $

This uses $G dot binom(m\/G, 2)$ within-gene pairs instead of the genome-wide
$binom(m, 2)$ -- a fraction $tilde 1\/G$ of the pairs -- so it isolates the
within-gene interaction signal and discards the (here null) cross-gene one. The pooling is *uniform* ($1\/G$ weights): every gene contributes equally
regardless of its pair count, which is the modelling assumption that all genes
share one epistasis variance $sigma_"gxg"^2$. 
\


== Relative error and SE for 1000 Random SNP, G=10 for pre-computed pooled kernal and mat-free method

We repeat the random-SNP comparison of @sec-identifiability with the *pooled*
kernel of the previous section: the $m = 1000$ SNPs are split into $G = 10$
equal 100-SNP genes, each gene contributes its own within-gene epistasis GRM,
and the ten are averaged into the single $W = 1/G sum_g K_g$ that the same
two-component model fits. @fig-pooled-prew and @fig-pooled-matfree show box
plots of the relative error of $sigma_"gxg"^2$ (estimate minus truth $= 0.2$)
with $m = 1000$ held fixed and $n$ increasing left-to-right. Because the pooled
$W u$ is computed by the identical operator in both pipelines, the two figures
should coincide up to Monte-Carlo error.

#figure(
  image("RandomSNP_pooled_preW_s2gxg0.2_s2e0.8_fixed_m1000_gxg.pdf", width: 100%),
  caption: [Precomputed pooled-$W$ MC-AI-REML on random SNPs ($G = 10$ genes,
    $m\/G = 100$ SNPs each): relative error of $sigma_"gxg"^2$, fixed $m = 1000$,
    increasing $n$ left-to-right. Each box shows the mean, SD, and a one-sample
    $t$-test of the relative error against $0$.],
) <fig-pooled-prew>

#figure(
  image("RandomSNP_pooled_matfree_s2gxg0.2_s2e0.8_fixed_m1000_gxg.pdf", width: 100%),
  caption: [Matrix-free MC-AI-REML with the pooled kernel on the same random-SNP
    phenotypes ($G = 10$, $m = 1000$): relative error of $sigma_"gxg"^2$,
    increasing $n$. The pooled $W u$ is applied gene-by-gene without ever forming
    a dense $W$.],
) <fig-pooled-matfree>


- *The two methods agree.* 

== Running time for 8.2

We time the pooled estimation step alone -- the `MC_REML` call, averaged over the
replicates -- with $m = 1000$ (G = 10) fixed and $n$ increasing.

#figure(
  image("time_comparison_pooled.pdf", width: 100%),
  caption: [Average estimation (MC-AI-REML) wall-clock time versus sample size
    $n$ at fixed $m = 1000$, $G = 10$, for the precomputed pooled-$W$ and
    matrix-free $W$-apply strategies.],
) <fig-time-pooled>

- *Precomputed-$W$ is faster throughout, but the gap closes with $n$.* At every
  tested $n$ the dense pooled-$W$ apply is cheaper in absolute time (from
  $0.39 "s"$ at $n = 1000$ to $186 "s"$ at $n = 16000$, versus $55 "s"$ to
  $1587 "s"$ matrix-free), yet its steeper exponent ($tilde n^(2.3)$ against
  $tilde n^(1.4)$ matrix-free) means the $tilde 140 times$ speed-up at $n = 1000$
  shrinks to $tilde 8 times$ at $n = 16000$. The dense $O(n^2)$ $W b$ grows from a
  tiny base while the matrix-free apply stays near-linear in $n$, so precomputed-$W$
  remains preferable in the small-$n$ / large-$m$ regime and matrix-free in the
  large-$n$, memory-bound regime.




== Relative error and SE for Contiguous SNP, G=10 for pre-computed pooled kernal



#figure(
  image("ContiguousSNP_pooled_preW_s2gxg0.2_s2e0.8_fixed_m1000_gxg.pdf", width: 100%),
  caption: [Precomputed pooled-$W$ on $m = 1000$ contiguous SNPs ($G = 10$ genes,
    $100$ SNPs each): relative error of $sigma_"gxg"^2$, increasing $n$
    left-to-right.],
) <fig-contig-1k>


#figure(
  image("chr1_10ksnp_pooled_preW_s2gxg0.2_s2e0.8_fixed_m10000_gxg.pdf", width: 100%),
  caption: [Precomputed pooled-$W$ on $m = 10000$ contiguous chromosome-1 SNPs
    ($G = 10$ genes, $1000$ SNPs each): relative error of $sigma_"gxg"^2$,
    increasing $n$.],
) <fig-chr1-pooled>

- *Near-unbiased at every $n$.* On contiguous (LD-correlated) genotypes the
  pooled estimator stays centred on $0$ throughout ($|"mean rel. err."| lt.eq
  0.007$ in both panels), with none of the small-$n$ inflation seen for random
  SNPs.

- *Precision tracks the within-gene pair count.* At $m = 1000$ ($100$ SNPs\/gene)
  the SD is $tilde 0.04$ and roughly flat in $n$ -- LD among neighbouring SNPs
  caps precision, so the spread stops shrinking past $n approx 4000$
  (@fig-contig-1k). At $m = 10000$ ($1000$ SNPs\/gene) the $tilde 100 times$
  larger within-gene pair count enriches $W$ and restores the clean
  $tilde 1\/sqrt(n)$ decay, $"SD"$ falling $0.066 arrow.r 0.011$
  (@fig-chr1-pooled).

