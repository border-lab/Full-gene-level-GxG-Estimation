#set page(margin: 1in)
#set text(size: 11pt)
#set par(justify: true)
#set heading(numbering: "1.1")

#align(center)[
  #text(size: 16pt, weight: "bold")[
    Monte-Carlo AI-REML for an Additive + Dominance Variance Model
  ]
  #v(4pt)
  #text(size: 10pt)[Matrix-free simulation and estimation pipeline (`Simulation_code_MCREML_two_var`)]
]

= Model

With $n$ individuals and $m$ SNPs, and no fixed effects (the phenotype $y$ is
mean-centred), the phenotype is decomposed into additive genetic, dominance
genetic, and residual parts:

$ y = g_a + g_d + e, quad
  g_a ~ N(0, sigma_a^2 K_a), quad
  g_d ~ N(0, sigma_d^2 K_d), quad
  e ~ N(0, sigma_e^2 I). $

The covariance of $y$ is therefore

$ V = "Var"(y) = sigma_a^2 K_a + sigma_d^2 K_d + sigma_e^2 I, $

with the two genomic relationship matrices (GRMs) built from thin design
matrices,

$ K_a = (Z_a Z_a^T) / m, quad K_d = (Z_d Z_d^T) / m . $

This is the two-component ($K_a$, $K_d$) generalisation of the additive-only
model $V = sigma_a^2 K_a + sigma_e^2 I$; the parameter vector is
$ theta = (sigma_a^2, sigma_d^2, sigma_e^2) $, so $k = 3$.

= Design matrices

*Additive design $Z_a$.* Column-standardised allele dosages: for SNP $j$ with
dosage column $x_j in {0,1,2}^n$,
$ (Z_a)_(i j) = (x_(i j) - macron(x)_j) / "sd"(x_j) . $

*Dominance design $Z_d$.* The orthogonal (statistical) dominance coding used by
GCTA / Zhu et al. (2015). With allele frequency $p_j = macron(x)_j \/ 2$ and
$q_j = 1 - p_j$, genotypes are recoded

$ 0 |-> -p_j / q_j, quad 1 |-> 1, quad 2 |-> -q_j / p_j, $

which is uncorrelated with the additive coding under Hardy–Weinberg
equilibrium. The columns are then standardised to mean $0$, variance $1$.

Both designs are standardised so that
$ "tr"(K_a) / n = "tr"(K_d) / n = 1 . $
This makes the simulation's per-component variance rescaling
(Section 3) unbiased regardless of LD structure, exactly as in the additive
pipeline.

= Simulation with sampling-error removal

The Cholesky factors $L_a, L_d$ satisfy
$ L_a L_a^T = sigma_a^2 K_a, quad L_d L_d^T = sigma_d^2 K_d, $
and are built once per $(text("genotype"), sigma_a^2, sigma_d^2)$. Each
replicate draws $u_1, u_2, u_3 ~ N(0, I_n)$ and forms

$ a = L_a u_1, quad d = L_d u_2, quad e = sqrt(sigma_e^2) u_3, $

then rescales each vector to its *exact* target sample variance
($a arrow.l a dot sqrt(sigma_a^2 \/ "var"(a))$, and likewise $d, e$) before
setting $y = a + d + e$ and centring. Because $"tr"(K_.)\/n = 1$ and the
designs are column-centred, $EE["var"(a)] = sigma_a^2$ etc., so the rescaling
removes Monte-Carlo sampling error without introducing bias.

= Matrix-free estimation

$K_a$, $K_d$, and $V$ are *never* formed. The only primitive is the $V$
mat-vec,

$ V B = (sigma_a^2 / m) Z_a (Z_a^T B)
      + (sigma_d^2 / m) Z_d (Z_d^T B)
      + sigma_e^2 B, $

costing $O(n m c)$ for $B in RR^(n times c)$. Systems $V X = B$ are solved by
batched conjugate gradient (one $V$-pass advances all columns), and traces are
Hutchinson estimates over $B$ Rademacher probes $r_1, ..., r_B$.

== AI-REML update

Writing $u = V^(-1) y$ and $K_1 = K_a$, $K_2 = K_d$, $K_3 = I$, each iteration
takes the average-information Newton step

$ theta arrow.l theta + (cal(A) + lambda I)^(-1) s, $

with score and average-information matrix

$ s_i = 1/2 (u^T K_i u - "tr"(V^(-1) K_i)), quad
  cal(A)_(i j) = 1/2 (K_i u)^T V^(-1) (K_j u) . $

The data quadratics are exact given $u$:
$ u^T K_a u = norm(Z_a^T u)^2 \/ m, quad
  u^T K_d u = norm(Z_d^T u)^2 \/ m, quad
  u^T I u = norm(u)^2, $
while the traces use the probes,
$ "tr"(V^(-1) K_i) approx 1/B sum_(b=1)^B (V^(-1) r_b)^T K_i r_b . $

Per iteration this is $1$ solve for $u$, $B$ shared solves for the probes, and
$3$ solves for $V^(-1)(K_j u)$ — total $O("iters" dot (B + 4) dot t_"cg" dot n m)$,
with no $n^2$ storage or $n^3$ factorisation. Estimates are clamped to be
non-negative.

= Pipeline

A 4-step SLURM dependency chain mirroring the additive pipeline:

+ *Cholesky* (1 job) — build and save $L_a, L_d$.
+ *Phenotype* (array) — draw $y$ per replicate; save $Z_a, Z_d$ once.
+ *MC-AI-REML* (array) — estimate $(sigma_a^2, sigma_d^2, sigma_e^2)$ per replicate.
+ *Combine + clean-up* — concatenate results, remove intermediates.

Results are written as `(s2a,s2d,s2e)` per replicate; `calc_stats.py` reports
mean, median, std, and a 95% CI for each component.
