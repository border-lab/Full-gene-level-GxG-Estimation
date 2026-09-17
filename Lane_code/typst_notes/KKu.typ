#set page(margin: 1in)
#set text(size: 11pt)
#set par(justify: true)
#set heading(numbering: "1.1")

#align(center)[
  #text(size: 16pt, weight: "bold")[
    Fast Apply for the Hadamard-Squared Gram Matrix $(K circle.small K) u$
  ]
  #v(4pt)
]

*Problem.* Given $X in RR^(n times m)$ of full column rank $m$, let
$K = X X^T in RR^(n times n)$ be its Gram matrix. We want to apply the
*element-wise square*
$ (K circle.small K)_(i j) = K_(i j)^2 $
to a vector, $u |-> (K circle.small K) u$, without ever forming the $n times n$
matrix $K$ --- and, if $m$ is large, without paying $O(n m^2)$ per apply either.

*Notation.* $circle.small$ is the element-wise (Hadamard) product,
$⊗$ the Kronecker product, and $⊙$ the *row-wise
Khatri--Rao* product: $(A ⊙ B)$ has row $i$ equal to
$A_(i,:) ⊗ B_(i,:)$. We write $x_i^T$ for row $i$ of $X$,
$"dg"(dot)$ for the diagonal of a matrix read as a vector, and $"diag"(dot)$ for
the matrix built from a vector.

= The operator has an exact, small feature map

The whole method rests on one line. Since $K_(i j) = x_i^T x_j$,
$ K_(i j)^2 = (x_i^T x_j)^2 = lr(⟨x_i ⊗ x_i, thin x_j ⊗ x_j⟩) , $
so the Hadamard square is itself a Gram matrix --- of the *tensor-squared rows*.
With $P = X ⊙ X in RR^(n times m^2)$,
$ K circle.small K = P P^T . $
This is an identity, not an approximation.

Only $binom(m+1, 2) = m(m+1)\/2$ of the $m^2$ columns are distinct, since
columns $(a,b)$ and $(b,a)$ coincide. Folding them together, define the
*symmetric* feature matrix $P in RR^(n times d)$, $d = m(m+1)\/2$, with columns
$ P_(dot.c, (a a)) = X_(dot.c a) circle.small X_(dot.c a), quad
  P_(dot.c, (a b)) = sqrt(2) thin X_(dot.c a) circle.small X_(dot.c b) quad (a < b) , $
for which $P P^T = K circle.small K$ still holds, because
$ sum_a (x_(i a) x_(j a))^2 + 2 sum_(a < b) x_(i a) x_(i b) x_(j a) x_(j b)
  = (sum_a x_(i a) x_(j a))^2 = K_(i j)^2 . $

Two consequences:

+ *Rank.* $"rank"(K circle.small K) <= min(n, d)$ with $d = m(m+1)\/2$, and
  $K circle.small K succ.eq 0$ (Schur product theorem). The operator is a
  low-rank PSD object whenever $m^2 lt.tilde n$.
+ *Exact apply.* $(K circle.small K) u = P (P^T u)$ costs $2 n d approx n m^2$
  flops and needs $O(n d)$ storage for $P$, built once.

= Matrix-free form: no $P$, no $K$

If $n d$ floats will not fit, the same product has a storage-free form. For any
$u$,
$ [K "diag"(u) K]_(i i) = sum_j K_(i j) u_j K_(j i) = sum_j u_j K_(i j)^2 , $
that is, $(K circle.small K) u = "dg"(K "diag"(u) K)$. Substituting $K = X X^T$
and collapsing the middle,
$ (K circle.small K) u = "dg"(X A X^T), quad
  A = X^T "diag"(u) X in RR^(m times m) . $
The diagonal is read off as a row-sum, $"dg"(X A X^T) = (X circle.small (X A)) bold(1)$,
so the recipe is: build $A$ ($n m^2$), form $X A$ ($n m^2$), one row-sum ($n m$).

#table(
  columns: (auto, 1.6fr, auto, auto, auto),
  align: (left, left, center, center, center),
  inset: 6pt,
  stroke: 0.5pt + luma(180),
  table.header([*Step*], [*Calculation*], [*Shape*], [*Time*], [*Space*]),
  [1], [$A = X^T ("diag"(u) X)$], [$m times m$], [$O(n m^2)$], [$O(m^2)$],
  [2], [$X A$], [$n times m$], [$O(n m^2)$], [$O(n m)$],
  [3], [$(X circle.small (X A)) bold(1)$], [$n times 1$], [$O(n m)$], [$O(n)$],
)

Total $2 n m^2$ flops in $O(n m)$ memory. Note what this says structurally: the
map factors through the $m times m$ matrix $A$, which is the rank statement of
Section 1 seen from the other side.

= Why a low-rank approximation is the right move

Take the SVD $X = U Sigma V^T$ and let $Lambda = Sigma^2$ hold the eigenvalues
$lambda_1 >= dots.c >= lambda_m$ of $K$. Each tensor-squared row factors as
$x_i ⊗ x_i = (V ⊗ V)(Sigma ⊗ Sigma)(u_i ⊗ u_i)$,
so $P = (U ⊙ U)(Sigma ⊗ Sigma)(V ⊗ V)^T$ and,
using $(V ⊗ V)^T (V ⊗ V) = I$,
$ K circle.small K = (U ⊙ U) thin (Lambda ⊗ Lambda) thin (U ⊙ U)^T . $

The weights are the *pairwise products* $lambda_a lambda_b$. (The columns of
$U ⊙ U$ are not orthonormal, so these are not literally the eigenvalues
of $K circle.small K$ --- but they govern its decay.) The message: *squaring a
kernel squares its spectral decay.* If $lambda_a$ decays, $lambda_a lambda_b$
decays much faster, so $K circle.small K$ is far more compressible than $K$
itself. This is why we compress rather than sketch blindly.

= Recommended method: truncate $X$, not $K circle.small K$

Because $X |-> P$ is *quadratic*, a rank-$r$ truncation of $X$ produces a
rank-$binom(r+1,2)$ truncation of the operator --- the compression is squared,
for free.

*Setup (once).*
+ Randomized SVD of $X$ to rank $r$: $tilde(X) = U_r Sigma_r in RR^(n times r)$,
  so $tilde(K) = tilde(X) tilde(X)^T$ is the best rank-$r$ approximation of $K$.
  Cost $O(n m r)$, one pass over $X$.
+ Build the symmetric Khatri--Rao features of the *small* factor,
  $tilde(P) = tilde(X) ⊙ tilde(X) in RR^(n times binom(r+1,2))$.
  Cost and storage $O(n r^2 \/ 2)$.

*Apply (per vector).*
$ (K circle.small K) u approx tilde(P) (tilde(P)^T u), quad approx n r^2 " flops" . $
If even $tilde(P)$ is unwelcome, use the Section 2 form at the compressed level:
$"dg"(tilde(X) tilde(A) tilde(X)^T)$ with $tilde(A) = tilde(X)^T "diag"(u) tilde(X) in RR^(r times r)$
--- $2 n r^2$ flops in $O(n r)$ memory.

Taking $r = m\/4$ already gives a $16 times$ speedup over the exact apply.

== Error bound

The bound is rigorous and costs nothing to state. Let $E = K - tilde(K) succ.eq 0$
be the truncation residual, so $norm(E)_2 = lambda_(r+1)(K) = sigma_(r+1)(X)^2$.
Expanding the Hadamard square and regrouping,
$ K circle.small K - tilde(K) circle.small tilde(K)
  = tilde(K) circle.small E + E circle.small tilde(K) + E circle.small E
  = tilde(K) circle.small E + E circle.small K . $
Schur's inequality $norm(A circle.small B)_2 <= (max_i A_(i i)) norm(B)_2$ for
PSD $A, B$ applies to each piece with the diagonal taken from the *un-truncated*
factor, giving
$ norm(K circle.small K - tilde(K) circle.small tilde(K))_2
  <= 2 (max_i K_(i i)) thin lambda_(r+1)(K)
  = 2 (max_i norm(x_i)^2) thin sigma_(r+1)(X)^2 . $
Both quantities are free: $max_i norm(x_i)^2$ is a row-norm scan, and
$sigma_(r+1)(X)$ falls out of the same randomized SVD used for the setup. Choose
$r$ by the observed decay of $sigma(X)$.

= Alternatives

*RPCholesky / Nyström on $K circle.small K$.* Attractive here because the two
primitives it needs are unusually cheap:
- a column: $(K circle.small K)_(dot.c, j) = (X x_j) circle.small (X x_j)$ --- $O(n m)$;
- the diagonal: $(K circle.small K)_(i i) = norm(x_i)^4$ --- $O(n m)$ for all of it.

With $c$ adaptively chosen pivots $J$ one gets $C = (K circle.small K)_(dot.c, J)$,
$W = C_(J, dot.c)$, and $K circle.small K approx C W^dagger C^T$: setup
$O(n m c + c^3)$, apply $O(n c + c^2)$. Unlike the method of Section 4 this
adapts to the actual spectrum of $K circle.small K$ instead of assuming $X$
truncates well --- prefer it when $sigma(X)$ has a flat tail but the rows are
incoherent.

*TensorSketch (Pham--Pagh).* A data-oblivious sketch of the degree-2 polynomial
feature map: count-sketch each copy and convolve via FFT, $O(m + s log s)$ per
row for $s$ output dimensions. Gives an unbiased estimate of $K_(i j)^2$ with
apply cost $O(n s)$, streaming and single-pass. Use when $X$ arrives as a stream
or is too large to revisit; it is otherwise dominated by the $X$-truncation
method, which exploits the structure this one deliberately ignores.

= Choosing among them

#table(
  columns: (1.3fr, auto, auto, auto),
  align: (left, center, center, center),
  inset: 6pt,
  stroke: 0.5pt + luma(180),
  table.header([*Regime / method*], [*Setup*], [*Apply*], [*Memory*]),
  [$n lt.tilde m^2$: form $K$, then $K circle.small K$],
    [$O(n^2 m)$], [$O(n^2)$], [$O(n^2)$],
  [$m$ small: exact $P = X ⊙ X$],
    [$O(n m^2 \/ 2)$], [$n m^2$], [$O(n m^2 \/ 2)$],
  [memory-bound: $"dg"(X A X^T)$, matrix-free],
    [---], [$2 n m^2$], [$O(n m)$],
  [*$n, m$ both large: truncate $X$ to rank $r$*],
    [$O(n m r)$], [$n r^2$], [$O(n r^2 \/ 2)$],
  [adaptive black box: RPCholesky],
    [$O(n m c + c^3)$], [$O(n c)$], [$O(n c)$],
  [streaming / oblivious: TensorSketch],
    [$O(n(m + s log s))$], [$O(n s)$], [$O(n s)$],
)

= If the target is a trace

When $(K circle.small K) u$ is only a means to $"tr"(K circle.small K)$ (or
$"tr"(f(K circle.small K))$), do not silently discard the tail. Split the
estimator Hutch++-style:
$ "tr"(K circle.small K) = "tr"(tilde(P)^T tilde(P))
  + "tr"((K circle.small K) - tilde(P) tilde(P)^T) , $
evaluating the first term deterministically on the captured subspace and running
Hutchinson probes only on the residual. The estimate stays *unbiased*, and
because the residual carries only the fast-decaying tail of Section 3, its
variance is lower than plain Monte Carlo at a fraction of the cost.

= Relation to the pairwise-epistasis operator

This is the unweighted special case of the $W u$ apply. In `Wu_complexity.typ`
the leading term is
$t_1 = (Z circle.small (Z(V circle.small M))) bold(1)$ with
$M = Z^T "diag"(u) Z$ --- exactly the Section 2 form, with a weight matrix $V$
inserted between the two contractions. Under independent, column-standardized
SNPs one has $V = bold(1) bold(1)^T - I$, and with $K = Z Z^T$,
$ t_1 = (K circle.small K) u - hat(Z) hat(Z)^T u, quad
  hat(Z) = Z circle.small Z in RR^(n times m) , $
the Hadamard square of the additive GRM with the $a = b$ terms removed. The
general weighted case is recovered by absorbing $sqrt(V_(a b)) > 0$ into the
columns of $P$, which leaves every method above unchanged.
