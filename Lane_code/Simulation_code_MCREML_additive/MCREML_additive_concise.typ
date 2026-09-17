#set page(margin: 1in)
#set text(size: 11pt)
#set par(justify: true)
#set heading(numbering: "1.1")

#align(center)[
  #text(size: 16pt, weight: "bold")[
    Monte-Carlo AI-REML for an Additive Model
  ]
  #v(4pt)
  #text(size: 10pt)[Matrix-free simulation and estimation pipeline (`Simulation_code_MCREML_additive`)]
]

= Model

With $n$ individuals, $m$ SNPs, and no fixed effects ($y$ mean-centred):

$ y = g_a + e, $
$ V = "Var"(y) = sigma_a^2 K + sigma_e^2 I, $

with parameter vector $theta = (sigma_a^2, sigma_e^2)$ so $k = 2$.

= Relationship matrix

*Additive* $K = Z Z^T \/ m$, with $Z$ the column-standardized dosages (each SNP
column centred and scaled to mean 0, variance 1). $K$ is never formed in
estimation; it enters only through the fast mat-vec $K b = Z (Z^T b) \/ m$ at
$O(n m)$.

= Exact AI-REML update

With no fixed effects ($y$ mean-centred, as in the Model section) the reference
estimator fits the model *exactly*: every $V^(-1)$ and every trace is formed
directly from a dense factorization, with no conjugate gradient and no stochastic
probes. It is the ground truth used to check the Monte-Carlo estimator
(unbiasedness, sampling variance vs. the Cramér--Rao bound) at small $n$.

Writing $x = V^(-1) y$ and $K_1 = K$, $K_2 = I$, each iteration forms
$V = sigma_a^2 K + sigma_e^2 I$ explicitly, factorizes the $n times n$ matrix
(Cholesky / direct inverse, $O(n^3)$), and takes the average-information Newton
step $theta arrow.l theta + (cal(A) + lambda I)^(-1) s$, clamped to
$sigma_i^2 gt.eq 0$, with

$ s_i = 1/2 (x^T K_i x - "tr"(V^(-1) K_i)), quad
  cal(A)_(i j) = 1/2 (K_i x)^T V^(-1) (K_j x) . $

The difference from the MC update is purely computational: here every piece is
exact. $x = V^(-1) y$ and $V^(-1)(K_j x)$ are direct solves against the
factorization, and each $"tr"(V^(-1) K_i)$ is the full contraction
$sum_(r s) (V^(-1))_(r s) (K_i)_(r s)$ of two dense $n times n$ matrices
($O(n^2)$) -- so the estimate is deterministic given $y$, with no Hutchinson
noise and no CG tolerance.

*Cost.* Per iteration, assembling $V$ from the prebuilt $K$ is $O(n^2)$, its
factorization $O(n^3)$, and each trace / $cal(A)_(i j)$ a further $O(n^2)$
contraction; the $O(n^3)$ factorization dominates. Exact AI-REML is thus
$O("iters" dot n^3)$ time and $O(n^2)$ memory, holding the dense $K, V, V^(-1)$
(with $K = Z Z^T \/ m$ built in $O(n^2 m)$ once). This is exactly the $O(n^3)$ /
$O(n^2)$ scaling that the matrix-free Monte-Carlo update below removes, at the
price of CG iterations and stochastic traces.

= MC AI-REML update

Writing $x = V^(-1) y$ and $K_1 = K$, $K_2 = I$, each iteration takes the
average-information Newton step
$ theta arrow.l theta + (cal(A) + lambda I)^(-1) s $ with

$ s_i = 1/2 (x^T K_i x - "tr"(V^(-1) K_i)), quad
  cal(A)_(i j) = 1/2 (K_i x)^T V^(-1) (K_j x) . $

Data quadratics are exact given $x$: $x^T K x = norm(Z^T x)^2 \/ m$,
$x^T I x = norm(x)^2$; the traces use $"Nmc"$ Rademacher probes
$U = [u_1, ..., u_"Nmc"] in RR^(n times "Nmc")$,
$"tr"(V^(-1) K_i) approx "Nmc"^(-1) sum_(b=1)^"Nmc" (V^(-1) u_b)^T K_i u_b$. Every
$V^(-1)$ is applied by batched CG whose only primitive is the $V$ mat-vec

$ V U = (sigma_a^2/m) Z (Z^T U) + sigma_e^2 U . $

The probes $U$ and their products $K U = Z(Z^T U) \/ m$ are drawn once and reused
across every iteration (common random numbers), which smooths the score and
stabilizes the Newton path.

= Cost

Let $n$ individuals, $m$ SNPs, $"Nmc"$ probes, $k = 2$ components, $"iters"$ REML
iterations, and $t_"cg"$ CG iterations per solve. `mc_reml` is fully matrix-free:
$K$ is never formed and the single CG primitive is the $V$ mat-vec on $c$ columns,

$ V U : quad O(n m c) quad ("additive " O(n m c) ", residual " O(n c)) . $

*Setup* (once per replicate):
+ Design $Z$ (column-standardize) -- $O(n m)$.
+ Draw Rademacher probes $U in {plus.minus 1}^(n times "Nmc")$ -- $O(n "Nmc")$.
+ Fixed-probe products $K U = Z(Z^T U) \/ m$ -- $O(n m "Nmc")$; built once, reused every iteration.

*Per REML iteration* ($times "iters"$):
+ Solve $x = V^(-1) y$, 1 RHS -- $O(t_"cg" thin n m)$.
+ Data quadratics $x^T K_i x$: $Z^T x$ at $O(n m)$, $norm(x)^2$ at $O(n)$.
+ Trace probes $P = V^(-1) U$, $"Nmc"$ RHS -- $O(t_"cg" thin n m thin "Nmc")$; then $"tr"(V^(-1) K_i) approx "Nmc"^(-1) sum_b P_b^T (K_i U)_b$ at $O(n "Nmc")$.
+ AI matrix: form $K_i x$ ($n times k$) at $O(n m)$; solve $G = V^(-1)(K_i x)$, $k$ RHS -- $O(t_"cg" thin n m thin k)$; then $cal(A) = 1/2 (K x)^T G$ at $O(n k^2)$.
+ AI-Newton step: solve the $k times k$ system -- $O(k^3)$, negligible.

The three solve groups carry $1 + "Nmc" + k = "Nmc" + 3$ RHS columns, so one
iteration costs $O(t_"cg" ("Nmc" + 3) thin n m)$ and the whole estimation is

$ O("iters" dot t_"cg" dot ("Nmc" + 3) dot n m) , $

with no $n^2$ or $n^3$ term and no GRM stored -- exactly the biobank-scale
scaling (linear in both $n$ and $m$) that the exact $O("iters" dot n^3)$
estimator above cannot reach.
