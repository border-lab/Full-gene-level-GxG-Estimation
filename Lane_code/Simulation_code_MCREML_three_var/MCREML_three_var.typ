#set page(margin: 1in)
#set text(size: 11pt)
#set par(justify: true)
#set heading(numbering: "1.1")

#align(center)[
  #text(size: 16pt, weight: "bold")[
    Monte-Carlo AI-REML for an Additive + Dominance + Epistasis Model
  ]
  #v(4pt)
  #text(size: 10pt)[Matrix-free simulation and estimation pipeline (`Simulation_code_MCREML_three_var`)]
]

= Model

With $n$ individuals, $m$ SNPs, and no fixed effects ($y$ mean-centred):

$ y = g_a + g_d + g_"gxg" + e, $
$ V = "Var"(y) = sigma_a^2 K_a + sigma_d^2 K_d + sigma_"gxg"^2 W + sigma_e^2 I, $

with parameter vector $theta = (sigma_a^2, sigma_d^2, sigma_"gxg"^2, sigma_e^2)$
so $k = 4$.

= Relationship matrices

*Additive* $K_a = Z_a Z_a^T \/ m$, with $Z_a$ the column-standardized dosages.

*Dominance* $K_d = Z_d Z_d^T \/ m$, with $Z_d$ the column-standardized GCTA
dominance coding $ {0,1,2} |-> {-p\/q, 1, -q\/p} $.

*Pairwise epistasis.* For each SNP pair $(a, b)$, $a < b$, form the element-wise
product $Z_a circle.small Z_b$ and standardize it to $H_(a b)$ (mean 0, variance
1). With $p = binom(m, 2) = m(m-1)\/2$ pairs,

$ W = 1/p sum_(a<b) H_(a b) H_(a b)^T . $


= Exact AI-REML update

With no fixed effects ($y$ mean-centred, as in the Model section) the reference
estimator (`gxg_reml/scripts/feasibility_sim.py`) fits the model *exactly*: every
$V^(-1)$ and every trace is formed directly from a dense factorization, with no
conjugate gradient and no stochastic probes. It is the ground truth used to check
the Monte-Carlo estimator (unbiasedness, sampling variance vs. the
Cramér--Rao bound) at small $n$.

Writing $x = V^(-1) y$ and $K_1 = K_a$, $K_2 = K_d$, $K_3 = W$, $K_4 = I$, each
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
$K_a, K_d, W, V, V^(-1)$ (with $K_a, K_d = Z Z^T \/ m$ built in $O(n^2 m)$ and $W$
in $O(n^2 p)$ once). This is exactly the $O(n^3)$ / $O(n^2)$ scaling that the
matrix-free Monte-Carlo update below removes, at the price of CG iterations and
stochastic traces.

= MC AI-REML update

Writing $x = V^(-1) y$ and $K_1 = K_a$, $K_2 = K_d$, $K_3 = W$, $K_4 = I$, each
iteration takes the average-information Newton step
$ theta arrow.l theta + (cal(A) + lambda I)^(-1) s $ with

$ s_i = 1/2 (x^T K_i x - "tr"(V^(-1) K_i)), quad
  cal(A)_(i j) = 1/2 (K_i x)^T V^(-1) (K_j x) . $

Data quadratics are exact given $x$:
$x^T K_a x = norm(Z_a^T x)^2\/m$, $x^T K_d x = norm(Z_d^T x)^2\/m$,
$x^T W x = x^T (W x)$, $x^T I x = norm(x)^2$; the traces use $"Nmc"$ Rademacher
probes $U = [u_1, ..., u_"Nmc"] in RR^(n times "Nmc")$,
$"tr"(V^(-1) K_i) approx "Nmc"^(-1) sum_(b=1)^"Nmc" (V^(-1) u_b)^T K_i u_b$. Every
$V^(-1)$ is applied by batched CG whose only primitive is the $V$ mat-vec

$ V U = (sigma_a^2/m) Z_a (Z_a^T U) + (sigma_d^2/m) Z_d (Z_d^T U)
        + sigma_"gxg"^2 (W U) + sigma_e^2 U . $


$W U$ is the standardized linear operator used in MoM.
= Cost


Let $n$ individuals, $m$ SNPs, $"Nmc"$ probes, $k = 4$ components, $"iters"$ REML
iterations, and $t_"cg"$ CG iterations per solve. `mc_reml` is fully matrix-free:
every epistasis apply $W v$ goes through `compute_WU` at $c_W = O(n m^2)$ -- a
dense $W$ is never formed or passed in. The single CG primitive is the $V$
mat-vec on $c$ columns,

$ V U : quad O((c_W + n m) thin c) quad
  ("additive/dominance " O(n m c) ", " W " " c_W thin c ", residual " O(n c)) . $

*Setup* (once per replicate):
+ Designs $Z_a, Z_d$ (column-standardize) -- $O(n m)$.
+ Weight matrices $S, R, T$ (`compute_weight_matrices`) -- $O(n m^2)$, stored $m times m$.
+ Draw Rademacher probes $U in {plus.minus 1}^(n times "Nmc")$ -- $O(n "Nmc")$.
+ Fixed-probe products $K_i U$, $i in {a, d, "gxg"}$ -- $O(n m "Nmc")$ (additive / dominance) $+ c_W "Nmc"$ (epistasis); built once, reused every iteration.

*Per REML iteration* ($times "iters"$):
+ Solve $x = V^(-1) y$, 1 RHS -- $O(t_"cg" (c_W + n m))$.
+ Data quadratics $x^T K_i x$: $Z_a^T x, Z_d^T x$ at $O(n m)$, $W x$ at $c_W$, $norm(x)^2$ at $O(n)$.
+ Trace probes $P = V^(-1) U$, $"Nmc"$ RHS -- $O(t_"cg" (c_W + n m) "Nmc")$; then $"tr"(V^(-1) K_i) approx "Nmc"^(-1) sum_b P_b^T (K_i U)_b$ at $O(n "Nmc")$.
+ AI matrix: form $K_i x$ ($n times k$) at $O(n m)$ ($W x$ reused); solve $G = V^(-1)(K_i x)$, $k$ RHS -- $O(t_"cg" (c_W + n m) k)$; then $cal(A) = 1/2 (K x)^T G$ at $O(n k^2)$.
+ AI-Newton step: solve the $k times k$ system -- $O(k^3)$, negligible.

The three solve groups carry $1 + "Nmc" + k = "Nmc" + 5$ RHS columns, so one
iteration costs $O(t_"cg" ("Nmc" + 5)(c_W + n m))$ and the whole estimation is

$ O("iters" dot t_"cg" dot ("Nmc" + 5) dot (c_W + n m))
  = O("iters" dot t_"cg" dot ("Nmc" + 5) dot n m^2) , $

to leading order, since the epistasis apply $c_W = O(n m^2)$ dominates the
additive / dominance terms $O(n m)$.

*A precomputed-$W$ speed-up is possible but not implemented here.* `compute_WU`
recomputes the $m$-by-$m$ contraction on every CG iteration, so at $m tilde n$ it
is $tilde m^2 \/ n$ times slower than applying a prebuilt dense $W$ would be
($O(n^2)$ per apply instead of $O(n m^2)$; $~2500 times$ at $n = m = 1000$). In
the current code `mc_reml` / `_v_matvec` take no $W$ argument and always rebuild
$(S, R, T)$ from $Z_a$, so estimation stays matrix-free. For $n gt.tilde m$ one
could form $W$ once (it is deterministic per genotype) and replace the
`compute_WU` in `_v_matvec` with a dense $W U$; the matrix-free apply remains
preferable in the large-$m$, memory-bound regime.
