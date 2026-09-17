#set page(
  paper: "a4",
  margin: (x: 2.2cm, y: 2.2cm),
  numbering: "1",
)
#set text(font: "New Computer Modern", size: 10.5pt)
#set par(justify: true, leading: 0.62em)
#set heading(numbering: "1.1")
#show heading: it => block(above: 1.1em, below: 0.6em)[#it]
#show raw.where(block: false): it => box(
  fill: luma(240), inset: (x: 3pt, y: 0pt), outset: (y: 3pt), radius: 2pt,
)[#it]
#show raw.where(block: true): it => block(
  fill: luma(245), inset: 8pt, radius: 3pt, width: 100%,
)[#it]

#align(center)[
  #text(size: 16pt, weight: "bold")[
    Average-Information REML for the Four-Component Model
  ]
  #v(0.3em)
  #text(size: 10pt)[Method note · #datetime.today().display("[year]-[month]-[day]")]
  #v(0.2em)
  #text(size: 10pt, style: "italic")[Simulation\_code\_AIREML pipeline · `Function_AIREML.py`]
]

#v(0.6em)
#line(length: 100%, stroke: 0.5pt)

This note derives, step by step, the exact average-information REML (AI-REML)
algorithm used to fit the four-component variance model, and gives the time and
space complexity of *every* step. The estimator is implemented in
`Function_AIREML.py` → `reml` (the solver) and `AI_REML` (the wrapper that
builds the GRMs). Throughout, $n$ = individuals, $m$ = SNPs, $k$ = number of
variance components (here $k = 4$), and `iters` = REML iterations.

= The model

We fit the linear mixed model with *no fixed effects* — the phenotype is
mean-centred up front ($bold(y) <- bold(y) - macron(y)$ in `simulate_remove_sampling_err`),
so $bb("E")[bold(y)] = bold(0)$ — and four independent random terms,

$ bold(y) = bold(g)_a + bold(g)_d + bold(g)_(g times g) + bold(e), $

whose covariance is the four-component structure

$ bold(V) = op("Var")(bold(y))
  = sigma_a^2 bold(K) + sigma_d^2 bold(D) + sigma_g^2 bold(W) + sigma_e^2 bold(I). $

To treat all components symmetrically, write the *parameter vector* and the
*component matrices* as

$ bold(theta) = (sigma_a^2, sigma_d^2, sigma_g^2, sigma_e^2)^top, quad quad
  (bold(K)_1, bold(K)_2, bold(K)_3, bold(K)_4) = (bold(K), bold(D), bold(W), bold(I)), $

so that $bold(V) = sum_(i=1)^k theta_i bold(K)_i$ and, crucially,

$ (partial bold(V)) / (partial theta_i) = bold(K)_i. $

== The three design matrices and their kernels

Every GRM follows the *same recipe* — a design matrix times its own transpose,
divided by the number of features it holds — applied to a different per-locus
coding. All start from the raw dosage matrix `SNP` $in {0, 1, 2}^(n times m)$,
with per-locus allele frequency $p_j = "mean"("SNP"_(dot j)) \/ 2$ and
$q_j = 1 - p_j$. The four components are:

#table(
  columns: (auto, auto, auto, auto),
  inset: 6pt,
  align: (left, left, left, left),
  stroke: 0.4pt + luma(180),
  table.header([*Component*], [*Design matrix*], [*Kernel*], [*Normalizer*]),
  [Additive $bold(K)$],
    [$bold(Z) in RR^(n times m)$ — std. dosage],
    [$bold(K) = bold(Z) bold(Z)' \/ m$], [$m$ = #box[SNPs]],
  [Dominance $bold(D)$],
    [$bold(Z)_D in RR^(n times L)$ — dominance recode],
    [$bold(D) = bold(Z)_D bold(Z)_D' \/ L$], [$L = m$ = #box[loci]],
  [G×G $bold(W)$],
    [$bold(H) in RR^(n times p)$ — pairwise products],
    [$bold(W) = bold(H) bold(H)' \/ p$], [$p = binom(m, 2)$ = #box[pairs]],
  [Residual $bold(I)$], [— (identity)], [$bold(I)$], [—],
)

The per-locus codings that build the columns are:

- *$bold(Z)$ (additive):* each column centred and scaled to unit variance,
  $bold(Z)_(dot j) = ("SNP"_(dot j) - macron(x)_j) \/ s_j$. One feature per SNP,
  so $m$ columns.

- *$bold(Z)_D$ (dominance):* the dominance *design matrix* — the analogue of
  $bold(Z)$, but each column uses a **heterozygote-contrast recode** (Zhu et al.
  2015), standardized with the theoretical HWE moments (centre $2 p_j^2$, scale
  $2 p_j q_j$; Hivert et al. 2021):
  $ (bold(Z)_D)_(i j)
    = (x_D - 2 p_j^2) / (2 p_j (1 - p_j)), quad
    x_D = cases(
      0 & "genotype " 0 " " (a a),
      2 p_j & "genotype " 1 " " (A a),
      4 p_j - 2 & "genotype " 2 " " (A A).
    ) $
  The recode is *flat-then-bumped*: it gives the heterozygote its own value
  rather than the homozygote average, capturing the non-additive within-locus
  signal. Under HWE it is uncorrelated with the matching additive column, so
  $bold(D)$ carries dominance and nothing else. One feature per locus, so
  *$L = m$*: the count is written $L$ (not $m$) only to signal "columns of
  $bold(Z)_D$," but here every SNP supplies exactly one. In code
  (`build_dominance_grm`), $bold(Z)_D$ is the array `Wd` and $bold(D)$ is
  `(Wd @ Wd.T) / L`.

- *$bold(H)$ (G×G):* one column per unordered SNP pair $(j, k)$, the
  element-wise product of the two standardized dosage columns, re-standardized:
  $bold(H)_(dot,(j k)) = "std"(bold(Z)_(dot j) circle.small bold(Z)_(dot k))$.
  There are $p = binom(m, 2) = m(m-1)/2$ pairs. Built by `build_W_batched`.

In code the four kernels are the list `Ks = [K, D, W, np.eye(n)]` passed to
`reml`; their construction costs are collected in §10.

= The log-likelihood (no fixed effects)

With no fixed effects to estimate, the phenotype is a mean-zero Gaussian,
$bold(y) tilde cal(N)(bold(0), bold(V))$. There is no $mu$ whose estimation
would spend a degree of freedom and bias the variance components downward, so no
REML error-contrast correction is needed — the maximum-likelihood estimator is
already unbiased for the (empty) fixed-effect part. The log-likelihood, up to an
additive constant, is the plain multivariate-normal one,

$ cal(L)(bold(theta)) = -1/2 (
  log|bold(V)| + bold(y)' bold(V)^(-1) bold(y)
). $

This is the restricted likelihood with an *empty* fixed-effect design
$bold(X)$: the correction term $log|bold(X)' bold(V)^(-1) bold(X)|$ and the
projection's rank-1 downdate both vanish. The average-information algorithm below
is therefore fitting *ML* variance components — still driven by the same
AI-REML machinery, only with the projection $bold(P)$ replaced everywhere by
$bold(V)^(-1)$.

= The covariance inverse (no projection)

Had the model carried a fixed-effect design $bold(X)$, the REML projection would
be $bold(P) = bold(V)^(-1) - bold(V)^(-1) bold(X) (bold(X)' bold(V)^(-1)
bold(X))^(-1) bold(X)' bold(V)^(-1)$. With $bold(X)$ empty the second term drops
and the projection collapses to a pure whitening,

$ bold(P) = bold(V)^(-1). $

So every $bold(P)$ that appears in the score and information below is simply
$bold(V)^(-1)$; there is no mean-removing rank-1 downdate. The code forms the
inverse once per iteration:

```python
Vi = np.linalg.inv(V + jitter*I)   # V^{-1}; the only cubic step
```

The inverse $bold(V)^(-1)$ is the sole genuinely cubic operation per sweep.

= The score (gradient)

Differentiate $cal(L)$ using two standard matrix identities,
$partial log|bold(V)| \/ partial theta_i = op("tr")(bold(V)^(-1) bold(K)_i)$ and
$partial bold(V)^(-1) \/ partial theta_i = - bold(V)^(-1) bold(K)_i bold(V)^(-1)$.
The $i$-th score is

$ s_i = (partial cal(L)) / (partial theta_i)
  = -1/2 (op("tr")(bold(V)^(-1) bold(K)_i) - bold(y)' bold(V)^(-1) bold(K)_i bold(V)^(-1) bold(y))
  = 1/2 (bold(y)' bold(V)^(-1) bold(K)_i bold(V)^(-1) bold(y) - op("tr")(bold(V)^(-1) bold(K)_i)). $

Two evaluation tricks make each score entry cheap. Let
$bold(u) = bold(V)^(-1) bold(y)$ be the whitened phenotype (one mat-vec). Then:

- The quadratic form is
  $bold(y)' bold(V)^(-1) bold(K)_i bold(V)^(-1) bold(y) = bold(u)' bold(K)_i bold(u) = (bold(K)_i bold(u))^top bold(u)$,
  needing only the mat-vec $bold(K)_i bold(u)$.

- The trace uses the element-wise identity (valid because $bold(V)^(-1)$,
  $bold(K)_i$ are symmetric):
  $ op("tr")(bold(V)^(-1) bold(K)_i) = sum_(a,b) (bold(V)^(-1))_(a b) (bold(K)_i)_(a b)
    = sum ( bold(V)^(-1) circle.small bold(K)_i ), $
  i.e. `np.sum(Vi * Ks[i])` — an $O(n^2)$ reduction, *not* an $O(n^3)$ matrix
  product. This is the key to keeping the per-iteration cost at one $O(n^3)$
  inverse rather than $k$ of them.

In code (`reml`):

```python
Viy   = Vi @ y                                   # u = V^{-1} y     O(n^2)
KiViy = [Ki @ Viy for Ki in Ks]                  # K_i u            O(k n^2)
score = np.array([
    0.5 * (Viy @ KiViy[i] - np.sum(Vi * Ks[i]))  # 0.5( u'K_i u - tr(V^-1 K_i) )
    for i in range(k)
])
```

= The information matrix: observed, expected, and *average*

A Newton-type update needs the second-order curvature. There are three choices.

*Observed information* (negative Hessian):

$ cal(H)_(i j) = -(partial^2 cal(L)) / (partial theta_i partial theta_j)
  = bold(y)' bold(V)^(-1) bold(K)_i bold(V)^(-1) bold(K)_j bold(V)^(-1) bold(y)
    - 1/2 op("tr")(bold(V)^(-1) bold(K)_i bold(V)^(-1) bold(K)_j). $

(Here $cal(H)$ is the negative Hessian — script, to avoid clashing with the G×G
design matrix $bold(H)$ of §1.1.)

*Expected (Fisher) information* (its expectation over $bold(y)$; the quadratic
term averages to half the trace):

$ cal(I)_(i j) = bb("E")[cal(H)_(i j)]
  = 1/2 op("tr")(bold(V)^(-1) bold(K)_i bold(V)^(-1) bold(K)_j). $

Both contain the term $op("tr")(bold(V)^(-1) bold(K)_i bold(V)^(-1) bold(K)_j)$,
which for each of the $k^2$ pairs costs an $O(n^3)$ matrix product (or two
mat-mats) — the expensive part of Fisher-scoring REML.

*Average information* is the arithmetic mean of the two. The awkward traces have
*opposite signs* and cancel exactly:

$ bold(A)_(i j) = 1/2 (cal(H)_(i j) + cal(I)_(i j))
  = 1/2 bold(y)' bold(V)^(-1) bold(K)_i bold(V)^(-1) bold(K)_j bold(V)^(-1) bold(y). $

This is the whole point of AI-REML: the curvature is now a *quadratic form* that
reuses vectors already computed for the score. With $bold(u) = bold(V)^(-1) bold(y)$
and $bold(K)_i bold(u)$ in hand,

$ bold(A)_(i j) = 1/2 (bold(K)_i bold(u))^top bold(V)^(-1) (bold(K)_j bold(u)), $

so each entry is an $O(n)$ dot product once $bold(V)^(-1)(bold(K)_j bold(u))$ is
formed ($k$ mat-vecs total). No $op("tr")(bold(V)^(-1) bold(K)_i bold(V)^(-1) bold(K)_j)$
is ever computed. In code (`reml`):

```python
ViKiViy = [Vi @ v for v in KiViy]                    # V^-1 K_i u       O(k n^2)
AI = np.array([
    [0.5 * (KiViy[i] @ ViKiViy[j]) for j in range(k)]  # 0.5 (K_i u)'(V^-1 K_j u)
    for i in range(k)
])
```

= The AI-REML update

Each iteration is a Newton / Fisher-scoring step with the average information as
the curvature:

$ bold(theta)^((t+1)) = bold(theta)^((t)) + bold(A)^(-1) bold(s), $

solved as a $k times k$ linear system (`np.linalg.solve(AI, score)`), then the
components are clamped to stay non-negative:

$ theta_i^((t+1)) <- max(theta_i^((t)) + [bold(A)^(-1) bold(s)]_i, " " epsilon), quad epsilon = 10^(-9). $

Convergence is declared when the largest absolute update falls below `tol`
($10^(-8)$). The `reml` defaults are `iters = 12`, `jitter = 1e-8`; the update is

```python
step = np.linalg.solve(AI + jitter*np.eye(k), score)   # k x k solve   O(k^3)
s    = np.clip(s + step, 1e-9, None)                    # non-negativity
if np.abs(step).max() < tol: break
```

= One full iteration, start to finish

Putting the pieces together, iteration $t$ of `reml` is:

+ *Assemble* $bold(V) = sum_i theta_i bold(K)_i$. #h(1fr) $O(k n^2)$
+ *Invert* $bold(V)^(-1) = (bold(V) + "jitter" bold(I))^(-1)$. #h(1fr) $O(n^3)$
+ *Whiten* $bold(u) = bold(V)^(-1) bold(y)$. #h(1fr) $O(n^2)$
+ *Mat-vecs* $bold(K)_i bold(u)$ and $bold(V)^(-1)(bold(K)_i bold(u))$ for all $i$. #h(1fr) $O(k n^2)$
+ *Score* $s_i = 1/2(bold(u)' bold(K)_i bold(u) - sum bold(V)^(-1) circle.small bold(K)_i)$. #h(1fr) $O(k n^2)$
+ *AI matrix* $bold(A)_(i j) = 1/2 (bold(K)_i bold(u))'(bold(V)^(-1) bold(K)_j bold(u))$. #h(1fr) $O(k^2 n)$
+ *Update* $bold(theta) <- op("clip")(bold(theta) + bold(A)^(-1) bold(s))$; test `tol`. #h(1fr) $O(k^3)$

Initialization (`reml`) splits the phenotypic variance equally,
$theta_i^((0)) = op("Var")(bold(y)) \/ k$, a robust starting point that keeps
$bold(V)$ positive-definite.

= Building the four GRMs (one-time, before the loop)

The kernels $bold(K), bold(D), bold(W)$ and their design matrices
$bold(Z), bold(Z)_D, bold(H)$ were defined in §1.1; here we record only their
one-time construction cost. They are built once in `AI_REML` and reused across
all `iters`.

#table(
  columns: (auto, auto, auto, auto),
  inset: 6pt,
  align: (left, left, center, center),
  stroke: 0.4pt + luma(180),
  table.header([*GRM*], [*Design matrix · normalizer*], [*Time*], [*Space*]),
  [Additive $bold(K)$], [$bold(Z) bold(Z)' \/ m$, #box[$m$ SNPs]], [$O(n^2 m)$], [$O(n^2)$],
  [Dominance $bold(D)$], [$bold(Z)_D bold(Z)_D' \/ L$, #box[$L = m$ loci]], [$O(n^2 m)$], [$O(n^2)$],
  [G×G $bold(W)$], [$bold(H) bold(H)' \/ p$, #box[$p = binom(m,2)$ pairs]], [$O(n^2 m^2)$], [$O(n^2 + n B)$],
)

The G×G kernel dominates the build: it sums $p tilde m^2 \/ 2$ rank-1 outer
products. `build_W_batched` processes the pairs in blocks of
$B$ = `pair_batch_size` columns, so peak memory is $O(n^2 + n B)$ instead of the
$O(n m^2)$ that materializing the full interaction matrix $bold(H)$ would need —
it trades passes over the data for memory.

= Time and space complexity — the whole estimator

*Per iteration.* The lone $O(n^3)$ operation is the inverse of $bold(V)$; every
other step is $O(k n^2)$ or cheaper. Hence

$ T_"iter" = O(n^3 + k n^2) = O(n^3) quad (n >> k), quad quad
  T_"REML" = O("iters" dot n^3). $

This is *tighter* than the naïve AI-REML bound $O("iters" dot k n^3)$: that bound
assumes each $bold(V)^(-1) bold(K)_i$ is formed as a dense matrix product. The
implementation avoids this by (i) computing $op("tr")(bold(V)^(-1) bold(K)_i)$ with
the element-wise sum ($O(n^2)$), and (ii) using only mat-*vec* products
$bold(K)_i bold(u)$, so the single factorization of $bold(V)$ is the only cubic
cost per sweep.

*One-time build.* Dominated by the G×G kernel, $O(n^2 m^2)$.

*Grand total.*

$ T_"total" = underbrace(O(n^2 m^2), "build" bold(W)) + underbrace(O("iters" dot n^3), "REML")
  quad+quad S_"total" = O(k n^2 + n B). $

The space is set by storing the $k$ dense GRMs simultaneously
(`Ks = [K, D, W, I]`), i.e. $O(k n^2)$, plus the batch buffer $O(n B)$ during the
$bold(W)$ build and the single working matrix $bold(V)^(-1)$ ($O(n^2)$).

#table(
  columns: (auto, auto, auto, auto),
  inset: 6pt,
  align: (left, left, center, center),
  stroke: 0.4pt + luma(180),
  table.header([*Phase*], [*Bottleneck*], [*Time*], [*Space*]),
  [Build $bold(K), bold(D)$], [dense matmul], [$O(n^2 m)$], [$O(n^2)$],
  [Build $bold(W)$], [$m^2 \/ 2$ pairs], [$O(n^2 m^2)$], [$O(n^2 + n B)$],
  [REML / iteration], [invert $bold(V)$], [$O(n^3)$], [$O(k n^2)$],
  [REML / total], [`iters` sweeps], [$O("iters" dot n^3)$], [$O(k n^2)$],
)

*Where each regime bites.* For the pipeline's simulation sizes ($n = 500$,
$m = 100$) the $O(n^2 m^2)$ G×G build dominates. As $n$ grows toward the
production panels ($n = 4000$), the $O("iters" dot n^3)$ REML inverse becomes the
wall. The two escape routes are the usual ones: build $bold(W)$ once and cache
it (already done), and — if $n$ grows further — replace the exact dense inverse
with stochastic-trace / conjugate-gradient AI-REML, which turns the per-iteration
$O(n^3)$ into $O(t_"cg" n^2)$ at the cost of approximate traces.

= Numerical safeguards (why the small constants are there)

- *`jitter` on $bold(V)$ and $bold(A)$* ($10^(-8)$): a ridge that keeps the
  inverse and the $k times k$ solve well-conditioned when a component collapses
  toward zero or two GRMs are nearly collinear (as $bold(D)$ and $bold(W)$ become
  under LD — see `four_variance_component.typ` §5).
- *`np.clip(..., 1e-9, None)`*: enforces $sigma_i^2 >= 0$. AI-REML is
  unconstrained, so without the clip a component can overshoot negative and
  render $bold(V)$ indefinite. Clipping is the simple projected-Newton fix.
- *Equal-share initialization*: guarantees $bold(V)^((0))$ is positive-definite
  regardless of the phenotype, avoiding a failed first inverse.
- *`tol` early stop*: AI-REML converges quadratically near the optimum, so the
  12-iteration cap is usually loose — the `max|step| < tol` test typically fires
  first.

= Summary

The estimator is exact dense average-information estimation for $bold(V) =
sigma_a^2 bold(K) + sigma_d^2 bold(D) + sigma_g^2 bold(W) + sigma_e^2 bold(I)$.
Because the model carries *no fixed effects* (the phenotype is centred so
$bb("E")[bold(y)] = bold(0)$), the REML projection $bold(P)$ collapses to
$bold(V)^(-1)$ and the fit is maximum likelihood. Each iteration
(i) inverts $bold(V)$, (ii) evaluates the score from one whitened residual
$bold(u) = bold(V)^(-1) bold(y)$ and $k$ mat-vecs, and (iii) builds the
average-information curvature as a quadratic form that *reuses* those mat-vecs —
sidestepping the $op("tr")(bold(V)^(-1) bold(K)_i bold(V)^(-1) bold(K)_j)$ that
makes Fisher scoring expensive. The result is $O("iters" dot n^3)$ estimation on
top of a one-time $O(n^2 m^2)$ G×G build, in $O(k n^2)$ memory.
