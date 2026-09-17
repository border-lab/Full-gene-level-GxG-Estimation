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
    Monte-Carlo Average-Information REML for the Additive Model
  ]
  #v(0.3em)
  #text(size: 10pt)[Method note · #datetime.today().display("[year]-[month]-[day]")]
  #v(0.2em)
  #text(size: 10pt, style: "italic")[
    Scalable variance-component estimation · BOLT-LMM–style stochastic AI-REML
  ]
]

#v(0.6em)
#line(length: 100%, stroke: 0.5pt)

This note derives, step by step, a *Monte-Carlo* average-information REML
(MC-AI-REML) algorithm for the additive-only variance model, and gives the time
and space complexity of *every* step. It is the scalable counterpart of the
exact dense estimator in `AIREML.typ`: the two fit the *same* model and take the
*same* Newton steps, but every dense $O(n^3)$ or $O(n^2)$ primitive is replaced
by a stochastic one that touches the genotypes only through matrix–vector
products. The design follows BOLT-LMM / BOLT-REML (Loh et al. 2015), which made
this feasible at biobank scale ($n tilde 10^5$–$10^6$). Throughout, $n$ =
individuals, $m$ = SNPs, $k$ = number of variance components (here $k = 2$),
$B$ = Monte-Carlo probe vectors, $t_"cg"$ = conjugate-gradient iterations per
solve, and `iters` = REML iterations. No code is given; the emphasis is the
mathematics and the cost of each step.

= The additive-only model

We fit the linear mixed model with *no fixed effects* — the phenotype is
mean-centred up front, so $bb("E")[bold(y)] = bold(0)$ — and two independent
random terms, an additive polygenic effect and residual noise,

$ bold(y) = bold(g)_a + bold(e), quad quad
  bold(g)_a tilde cal(N)(bold(0), sigma_a^2 bold(K)), quad
  bold(e) tilde cal(N)(bold(0), sigma_e^2 bold(I)), $

whose covariance is the two-component structure

$ bold(V) = op("Var")(bold(y)) = sigma_a^2 bold(K) + sigma_e^2 bold(I). $

Writing the *parameter vector* and *component matrices* symmetrically,

$ bold(theta) = (sigma_a^2, sigma_e^2)^top, quad quad
  (bold(K)_1, bold(K)_2) = (bold(K), bold(I)), $

gives $bold(V) = sum_(i=1)^k theta_i bold(K)_i$ and
$partial bold(V) \/ partial theta_i = bold(K)_i$, exactly as in the exact note.

The additive GRM is the *standardized-dosage cross-product*

$ bold(K) = bold(Z) bold(Z)^top \/ m, $

where $bold(Z) in RR^(n times m)$ holds the per-SNP standardized dosages
($bold(Z)_(dot j) = ("SNP"_(dot j) - macron(x)_j) \/ s_j$). With no fixed
effects the REML projection collapses to a pure whitening, $bold(P) = bold(V)^(-1)$
(see `AIREML.typ` §3), so every $bold(P)$ below is simply $bold(V)^(-1)$.

= The optimization objective

AI-REML is a numerical *optimizer*; before scaling it we state exactly what it
optimizes. With the mean-zero Gaussian model of §1,
$bold(y) tilde cal(N)(bold(0), bold(V)(bold(theta)))$, the estimator maximizes
the log-likelihood over the *non-negative* variance components,

$ hat(bold(theta)) = op("arg max")_(bold(theta) >= bold(0)) thin cal(L)(bold(theta)),
  quad
  cal(L)(bold(theta)) = -1/2 (log|bold(V)(bold(theta))|
                              + bold(y)^top bold(V)(bold(theta))^(-1) bold(y)),
  quad bold(V)(bold(theta)) = sigma_a^2 bold(K) + sigma_e^2 bold(I) $

(up to the constant $-n/2 log 2 pi$). This is a smooth but *non-convex* program on
the positive orthant $bold(theta) = (sigma_a^2, sigma_e^2) >= bold(0)$. Because the
model has no fixed effect this restricted likelihood is the plain one and
$bold(P) = bold(V)^(-1)$ (§1).

*Stationarity — the estimating equations.* A maximizer sets the *score*
(gradient) to zero:

$ s_i (bold(theta)) = (partial cal(L)) / (partial theta_i)
  = 1/2 (bold(y)^top bold(V)^(-1) bold(K)_i bold(V)^(-1) bold(y)
         - op("tr")(bold(V)^(-1) bold(K)_i)) = 0, quad i = 1, ..., k. $

Read one component at a time, each equation says: at the optimum the *observed*
generalized sum of squares $bold(y)^top bold(V)^(-1) bold(K)_i bold(V)^(-1) bold(y)$
must equal the *model-expected* one $op("tr")(bold(V)^(-1) bold(K)_i)$.

*The solver.* There is no closed form, so AI-REML drives
$bold(s)(bold(theta)) = bold(0)$ by a Newton-type iteration. It uses the
*average-information* matrix $bold(A)$ — the mean of the observed and Fisher
(expected) information — as a positive-(semi)definite curvature surrogate for
$-partial^2 cal(L)$ (derived in `AIREML.typ` §5),

$ bold(A)_(i j) = 1/2 bold(y)^top bold(V)^(-1) bold(K)_i bold(V)^(-1) bold(K)_j bold(V)^(-1) bold(y)
  approx - (partial^2 cal(L)) / (partial theta_i partial theta_j), quad quad
  bold(theta)^((t+1)) = bold(theta)^((t)) + bold(A)^(-1) bold(s). $

So the *entire* computational task at each iterate is to evaluate two things from
the current $bold(theta)$: the gradient $bold(s)$ and the curvature $bold(A)$. Both
are assembled from just two kinds of quantities,

- *quadratic forms in $bold(y)$* —
  $bold(y)^top bold(V)^(-1) bold(K)_i bold(V)^(-1) bold(y)$ and
  $bold(y)^top bold(V)^(-1) bold(K)_i bold(V)^(-1) bold(K)_j bold(V)^(-1) bold(y)$ —
  each obtained from the single linear solve $bold(u) = bold(V)^(-1) bold(y)$ plus
  a few mat-vecs, and

- *one trace* — $op("tr")(bold(V)^(-1) bold(K)_i)$.

Everything but the trace is a handful of linear solves against $bold(V)$; the
trace is the only term that couples to the *entire spectrum* of $bold(V)^(-1)$.
That asymmetry is the whole story of the next section.

= Why a stochastic method helps

Two observations about the objective of §2 turn an intractable computation into a
scalable one.

*(1) The gradient is "observed minus expected."* Under the model
$bold(y) tilde cal(N)(bold(0), bold(V))$,

$ bb("E")_bold(y) [bold(y)^top bold(V)^(-1) bold(K)_i bold(V)^(-1) bold(y)]
  = op("tr")(bold(V)^(-1) bold(K)_i bold(V)^(-1) bb("E")[bold(y) bold(y)^top])
  = op("tr")(bold(V)^(-1) bold(K)_i bold(V)^(-1) bold(V))
  = op("tr")(bold(V)^(-1) bold(K)_i). $

So the trace *is exactly the expectation* of the data quadratic term, and the
score $s_i = 1/2 ("observed" - "expected")$ pushes each observed sum of squares to
its model expectation. The one expensive term in the gradient is an *expectation*
— and expectations are precisely what Monte Carlo estimates cheaply.

*(2) The trace is the only piece that needs the whole inverse.* The quadratic
forms need $bold(V)^(-1)$ applied to a *few specific vectors* ($bold(y)$, and
$bold(K)_j bold(u)$), each a single linear solve. But
$op("tr")(bold(V)^(-1) bold(K)_i) = sum_a (bold(V)^(-1) bold(K)_i)_(a a)$ sums all
$n$ diagonal entries — naïvely the *full* inverse ($O(n^3)$ time, $O(n^2)$ memory)
or all $n$ eigenvalues. For the exact dense method that trace, together with the
GRM and inverse it rides on, is the wall:

- *Storage.* $bold(K)$ and $bold(V)^(-1)$ are $n times n$. At $n = 5 times 10^5$
  each is $2.5 times 10^11$ doubles ($tilde 2$ TB) — they cannot be formed at all.
- *Compute.* One dense inverse is $O(n^3) tilde 10^17$ flops per REML iteration.

*The stochastic fix.* Any trace is an expectation over random probes,

$ op("tr")(bold(M)) = bb("E")_bold(r) [bold(r)^top bold(M) bold(r)],
  quad bb("E")[bold(r) bold(r)^top] = bold(I), $

and an expectation is estimated by a *sample average*. Replacing the exhaustive
$n$-term diagonal sum with $B << n$ Hutchinson probes turns the trace into $B$
extra linear solves $bold(V)^(-1) bold(r)_b$, with error $O(1 \/ sqrt(B))$
*independent of $n$*. This is the single reason a stochastic method helps: it
converts the lone spectral quantity in the gradient — the trace — from an exact
$O(n^3)$ computation into a Monte-Carlo average of $O(n m)$ linear solves. The
symmetry is pleasing: the score is observed-minus-expected, and Monte Carlo
estimates the expected term by simulating the very randomness it averages over;
BOLT-REML sharpens this by drawing the probes from the fitted model so their
covariance matches $bold(V)$, cutting the variance at fixed $B$.

Once the trace is stochastic, the two remaining dense objects — the GRM $bold(K)$
and the inverse $bold(V)^(-1)$ — are removed by two *matrix-free* primitives, so
that every operation reduces to a pass over the thin genotype matrix $bold(Z)$
(kept in memory or streamed from disk):

#table(
  columns: (auto, auto, auto, auto),
  inset: 6pt,
  align: (left, left, left, center),
  stroke: 0.4pt + luma(180),
  table.header([*Exact primitive*], [*Dense cost*], [*MC replacement*], [*Cost*]),
  [Form GRM $bold(K) = bold(Z Z)^top \/ m$], [$O(n^2 m)$, $O(n^2)$ mem],
    [implicit mat-vec $bold(K) bold(b) = bold(Z)(bold(Z)^top bold(b)) \/ m$],
    [$O(n m)$],
  [Invert $bold(V)^(-1) bold(b)$], [$O(n^3)$],
    [conjugate gradient solve $bold(V) bold(x) = bold(b)$], [$O(t_"cg" n m)$],
  [Trace $op("tr")(bold(V)^(-1) bold(K)_i)$], [needs full $bold(V)^(-1)$],
    [Hutchinson probes $frac(1, B) sum_b bold(r)_b^top bold(V)^(-1) bold(K)_i bold(r)_b$],
    [$O(B t_"cg" n m)$],
)

The rest of this note develops the three ingredients (§4–§6), assembles the
Monte-Carlo score and average-information matrix (§7–§8), and totals the cost
(§10–§11).

= Ingredient 1 — the implicit GRM mat-vec (never form $bold(K)$)

$bold(K)$ enters the algorithm *only* as $bold(K) bold(b)$ for various vectors
$bold(b)$. Associativity turns this into two thin products against $bold(Z)$:

$ bold(K) bold(b) = frac(1, m) bold(Z) (bold(Z)^top bold(b)), quad quad
  bold(V) bold(b) = frac(sigma_a^2, m) bold(Z)(bold(Z)^top bold(b)) + sigma_e^2 bold(b). $

Computing $bold(Z)^top bold(b)$ (an $m$-vector) then $bold(Z)(dot)$ (an $n$-vector)
costs $O(n m)$ time and *no* $n times n$ storage — only the $n$- and $m$-vectors
plus $bold(Z)$ itself. Each $bold(V)$-multiply is thus **two streamed passes over
the genotypes**. This single reformulation removes both the $O(n^2)$ memory and
the $O(n^2 m)$ GRM build of the exact method; everything downstream is expressed
through it.

= Ingredient 2 — conjugate gradient in place of the inverse

Wherever the exact algorithm needs $bold(V)^(-1) bold(b)$, we instead solve the
symmetric positive-definite system

$ bold(V) bold(x) = bold(b) $

with *conjugate gradient* (CG), which uses only the $bold(V)$-mat-vec of §4 — one
mat-vec per CG iteration, so $O(t_"cg" n m)$ per solve. CG needs

$ t_"cg" = O(sqrt(kappa) log(1 \/ epsilon_"cg")) quad "iterations for relative residual " epsilon_"cg", $

where $kappa = kappa(bold(V))$ is the condition number. Here $kappa$ is *tame*
because the residual term $sigma_e^2 bold(I)$ acts as a Tikhonov ridge: the
eigenvalues of $bold(V)$ are $sigma_a^2 lambda_j (bold(K)) + sigma_e^2 >= sigma_e^2 > 0$, so

$ kappa(bold(V)) = frac(sigma_a^2 lambda_max(bold(K)) + sigma_e^2,
                        sigma_a^2 lambda_min(bold(K)) + sigma_e^2)
  <= 1 + frac(sigma_a^2, sigma_e^2) lambda_max(bold(K)). $

A moderate heritability keeps $sigma_a^2 \/ sigma_e^2 = O(1)$, so $t_"cg"$ is a
few tens of iterations rather than $n$. Two standard accelerations apply: a cheap
*Jacobi (diagonal) preconditioner* to shrink $kappa$, and **warm-starting** each
solve from the previous REML iteration's solution, since $bold(V)$ changes only
slightly between sweeps.

= Ingredient 3 — stochastic trace estimation (Hutchinson / BOLT probes)

The one quantity that genuinely seems to need the full inverse is the trace
$op("tr")(bold(V)^(-1) bold(K)_i)$ in the score. The *Hutchinson estimator*
removes it: for any square $bold(M)$ and random probe $bold(r)$ with
$bb("E")[bold(r) bold(r)^top] = bold(I)$,

$ op("tr")(bold(M)) = bb("E")[bold(r)^top bold(M) bold(r)]
  quad ("since " bb("E")[bold(r)^top bold(M) bold(r)]
        = op("tr")(bold(M) thin bb("E")[bold(r) bold(r)^top]) = op("tr")(bold(M))). $

Averaging $B$ independent probes gives the unbiased estimate, applied with
$bold(M) = bold(V)^(-1) bold(K)_i$:

$ hat(op("tr"))(bold(V)^(-1) bold(K)_i)
  = frac(1, B) sum_(b=1)^B bold(r)_b^top bold(V)^(-1) bold(K)_i bold(r)_b
  = frac(1, B) sum_(b=1)^B (bold(V)^(-1) bold(r)_b)^top (bold(K)_i bold(r)_b). $

Implementation, using §4–§5: solve $bold(w)_b = bold(V)^(-1) bold(r)_b$ once per
probe by CG ($B$ solves total, *shared across all $k$ components*), form
$bold(K)_i bold(r)_b$ by the implicit mat-vec, and dot. Key points:

- *Probe distribution.* Rademacher probes ($bold(r)_b in {plus.minus 1}^n$
  i.i.d.) satisfy $bb("E")[bold(r) bold(r)^top] = bold(I)$ and have *lower
  variance* than Gaussian probes for trace estimation (Hutchinson 1990).

- *Accuracy.* The estimator is unbiased with relative standard error
  $O(1 \/ sqrt(B))$; $B$ is a modest constant (tens), independent of $n$. This
  $1 \/ sqrt(B)$ floor is the *only* approximation the method introduces.

- *BOLT variance reduction.* BOLT-REML draws its stochastic probes to *match the
  model structure* — building random phenotypes from the current variance
  components so that the estimated quadratic forms track the quantities they
  approximate — which shrinks the Monte-Carlo variance at fixed $B$ relative to
  naïve i.i.d. probes.

= The Monte-Carlo score (gradient)

The score $bold(s)(bold(theta))$ of §2 is the object we must evaluate at each
iterate:

$ s_i = frac(1, 2) (bold(y)^top bold(V)^(-1) bold(K)_i bold(V)^(-1) bold(y)
                    - op("tr")(bold(V)^(-1) bold(K)_i)). $

Only its *evaluation* changes here. Let $bold(u) = bold(V)^(-1) bold(y)$ (*one* CG solve).

- *Data term — exact.* Reusing $bold(u)$,
  $ bold(y)^top bold(V)^(-1) bold(K)_i bold(V)^(-1) bold(y) = bold(u)^top bold(K)_i bold(u), quad
    bold(u)^top bold(K) bold(u) = frac(1, m) (bold(Z)^top bold(u))^top (bold(Z)^top bold(u)), quad
    bold(u)^top bold(I) bold(u) = bold(u)^top bold(u), $
  each an $O(n m)$ (additive) or $O(n)$ (residual) reduction — no extra solve.

- *Trace term — Monte-Carlo.* Replace $op("tr")(bold(V)^(-1) bold(K)_i)$ by the
  Hutchinson estimate of §6 ($B$ shared solves).

Thus the score is exact in its data-dependent part and stochastic *only* through
the trace: all of its Monte-Carlo noise has standard deviation $O(1 \/ sqrt(B))$.

= The Monte-Carlo average-information matrix

The average information reuses the AI algebra of the exact note,

$ bold(A)_(i j) = frac(1, 2) (bold(K)_i bold(u))^top bold(V)^(-1) (bold(K)_j bold(u)), $

which needs $bold(g)_j = bold(V)^(-1) (bold(K)_j bold(u))$ — one CG solve per
component, $k$ solves total — after which each entry is an $O(n)$ dot product:

$ bold(A)_(i j) = frac(1, 2) (bold(K)_i bold(u))^top bold(g)_j. $

Because this is a *quadratic form in the data* (not a trace), it needs **no**
Monte-Carlo probes: given $bold(u)$ it is computed exactly from $k$ solves. Here
$bold(A)$ is a $2 times 2$ matrix, and its inverse in the update is $O(k^3) = O(1)$.
(One may also estimate $bold(A)$ stochastically from the same probes, à la
BOLT; we keep it exact since $k$ extra solves are cheap.)

= The AI-REML update

Each iteration is the same Newton / Fisher-scoring step as the exact method,

$ bold(theta)^((t+1)) = op("clip")(bold(theta)^((t)) + bold(A)^(-1) bold(s), thin epsilon), quad epsilon = 10^(-9), $

a $k times k$ solve followed by the non-negativity clamp. Two MC-specific points:

- *Common random numbers.* Fix the probe seed (reuse the *same* $bold(r)_b$
  across REML iterations). This makes the stochastic score a smooth
  deterministic function of $bold(theta)$, so the Newton path does not jitter and
  convergence behaves like the exact method up to the $O(1 \/ sqrt(B))$ floor.

- *Noise floor.* With finite $B$ the estimate converges to a small neighbourhood
  of the true optimum whose radius scales as $1 \/ sqrt(B)$; increasing $B$ (or a
  final large-$B$ polishing step) tightens it.

= One full iteration, start to finish

Putting the pieces together, iteration $t$ of MC-AI-REML is:

+ *Assemble* $bold(V)$ *implicitly* — just store $(sigma_a^2, sigma_e^2)$; no
  matrix is formed. #h(1fr) $O(1)$
+ *Whiten* $bold(u) = bold(V)^(-1) bold(y)$ by CG. #h(1fr) $O(t_"cg" n m)$
+ *Data quadratics* $bold(u)^top bold(K)_i bold(u)$ via $bold(Z)^top bold(u)$. #h(1fr) $O(n m)$
+ *Trace probes:* for $b = 1..B$, CG-solve $bold(w)_b = bold(V)^(-1) bold(r)_b$,
  form $bold(K)_i bold(r)_b$, accumulate. #h(1fr) $O(B t_"cg" n m)$
+ *AI solves* $bold(g)_j = bold(V)^(-1)(bold(K)_j bold(u))$ for all $j$. #h(1fr) $O(k t_"cg" n m)$
+ *Score* $s_i = frac(1,2)(bold(u)^top bold(K)_i bold(u) - hat(op("tr"))_i)$. #h(1fr) $O(k n m)$
+ *AI matrix* $bold(A)_(i j) = frac(1,2)(bold(K)_i bold(u))^top bold(g)_j$. #h(1fr) $O(k^2 n)$
+ *Update* $bold(theta) <- op("clip")(bold(theta) + bold(A)^(-1) bold(s))$; test `tol`. #h(1fr) $O(k^3)$

The dominant term is the $(1 + B + k)$ CG solves — steps 2, 4, 5 — each a stack
of $t_"cg"$ genotype passes. Initialization splits the phenotypic variance
equally, $theta_i^((0)) = op("Var")(bold(y)) \/ k$, keeping $bold(V)$
positive-definite for the first CG solve.

= Time and space complexity — the whole estimator

*Per iteration.* Every step reduces to CG solves and thin genotype mat-vecs.
Counting $(1 + B + k)$ solves at $O(t_"cg" n m)$ each,

$ T_"iter" = O((B + k) thin t_"cg" thin n m), quad quad
  T_"MC-REML" = O("iters" dot (B + k) thin t_"cg" thin n m). $

There is *no* $n^3$ term and *no* $n^2$ term anywhere. Contrast the exact dense
method, $O("iters" dot n^3)$: the MC estimator wins whenever
$(B + k) thin t_"cg" thin m << n^2$, which holds decisively once $n$ reaches the
$10^5$–$10^6$ range — and, unlike the dense method, it *runs at all*, because it
never allocates an $n times n$ array.

*One-time build.* None: the GRM is never materialized. The genotypes $bold(Z)$
are the only large object, read (or streamed) once per mat-vec.

*Space.* $O(n m)$ if $bold(Z)$ is held in memory, or $O(n b)$ with a streamed
block of $b$ SNPs (BOLT's mode), plus $O(B n)$ for the probes/solutions and
$O(n)$ scratch — never $O(n^2)$.

#table(
  columns: (auto, auto, auto, auto),
  inset: 6pt,
  align: (left, left, center, center),
  stroke: 0.4pt + luma(180),
  table.header([*Phase*], [*Bottleneck*], [*Time*], [*Space*]),
  [GRM], [implicit — none built], [—], [$O(n m)$ or $O(n b)$],
  [1 CG solve], [$t_"cg"$ genotype passes], [$O(t_"cg" n m)$], [$O(n)$],
  [MC-REML / iteration], [$(1 + B + k)$ solves], [$O((B + k) t_"cg" n m)$], [$O(n m + B n)$],
  [MC-REML / total], [`iters` sweeps], [$O("iters" (B + k) t_"cg" n m)$], [$O(n m + B n)$],
)

*The BOLT headline.* Each CG iteration is one pass over the $n times m$
genotypes, so the whole estimator costs a fixed number of data passes —
roughly $"iters" times (1 + B + k) times t_"cg"$ of them — making the runtime
*linear in the data size* $n m$ rather than cubic in $n$. That linearity is what
lets BOLT-LMM / BOLT-REML fit variance components on hundreds of thousands of
individuals.

= Accuracy, variance reduction, and numerical safeguards

- *Probe count $B$.* Sets the trace variance: relative error $O(1 \/ sqrt(B))$.
  Modest $B$ (tens) suffices for the gradient; a larger $B$ only near convergence
  polishes the final estimate.
- *Rademacher probes* over Gaussian: same unbiasedness, smaller variance.
- *Common random numbers* across REML iterations: turns the noisy score into a
  smooth surrogate objective, stabilizing the Newton path.
- *CG tolerance $epsilon_"cg"$*: may be loose in early REML sweeps (the step is
  approximate anyway) and tightened as $bold(theta)$ settles — an inner/outer
  accuracy trade-off. Warm-starting each CG solve from the previous sweep cuts
  $t_"cg"$ further.
- *Conditioning.* The $sigma_e^2 bold(I)$ ridge bounds $kappa(bold(V))$ (§5), so
  CG converges quickly and a light Jacobi preconditioner usually suffices; no
  extra jitter is needed on $bold(V)$.
- *Non-negativity clip* $sigma_i^2 >= 10^(-9)$: as in the exact method, keeps
  $bold(V)$ positive-definite so the next CG solve is well-posed.

= Summary

MC-AI-REML fits the additive model $bold(V) = sigma_a^2 bold(K) + sigma_e^2 bold(I)$
by taking the *same* average-information Newton steps as exact AI-REML, but
replacing its three dense primitives — form $bold(K)$, invert $bold(V)$, trace
$bold(V)^(-1) bold(K)_i$ — with (i) an implicit genotype mat-vec
$bold(K) bold(b) = bold(Z)(bold(Z)^top bold(b)) \/ m$, (ii) conjugate-gradient
solves of $bold(V) bold(x) = bold(b)$, and (iii) a Hutchinson stochastic trace
over $B$ Rademacher probes. Each iteration costs $(1 + B + k)$ CG solves at
$O(t_"cg" n m)$ apiece, so the estimator is $O("iters" (B + k) t_"cg" n m)$ time
and $O(n m + B n)$ space — linear in the data and free of any $n times n$ object.
The single approximation is the $O(1 \/ sqrt(B))$ trace noise, controlled by $B$,
common random numbers, and a final polishing step. This is the algorithmic core
that lets BOLT-LMM estimate additive variance components at biobank scale.
