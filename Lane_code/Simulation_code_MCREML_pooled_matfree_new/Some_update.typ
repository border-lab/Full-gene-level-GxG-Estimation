
= Basic Framework

== MC AI-REML 

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



== 1. Oringinal method (extended from single-gene case): $O(n m^2)$ 
\
This method build linear operator $W u$ implicitly using one random vector (matrix) per apply.
\
*Setup (once per gene).* From $Z$:
+ the weight matrices $B, R, T in RR^(m times m)$, symmetric with zero diagonal
$ B_(i j) = 1 / sigma^2_(i j), quad R_(i j) = mu_(i j) / sigma^2_(i j), quad T_(i j) = mu_(i j)^2 / sigma^2_(i j) $
where empirical $mu_(i j) = bb(E)[Z_i Z_j]$ and $sigma^2_(i j) = "Var"(Z_i Z_j) - O(n m^2)$ 
+ the vector $v_R = (Z circle.small (Z R)) bold(1) $ ---
  $O(n m^2)$;
+ the scalar $s_T = bold(1)^T T bold(1) in RR$ --- $O(m^2)$.
Only $Z$, $V$, $v_R$ and $s_T$ stay; $R$, $T$ and the $n times m$ intermediate
$Z R$ can all be discarded.
*Apply (per vector $u$).* Exactly one $u$-dependent object is built:
$ s_u = bold(1)^T u, quad M = Z^T "diag"(u) Z $
and the result is
$ W u = 1/(2 p) [
  (Z circle.small (Z (B circle.small M))) bold(1) - s_u thin v_R + (s_u s_T - u^T v_R) bold(1)
  ]. $

== 1.1 Opitimization: can we do it using LD score rather than weight matrices



= LD Score

Let the reference panel contain $n$ individuals and $m$ SNPs. For variant $i$, let $p_i$
denote the allele frequency and $g_i in {0, 1, 2}$ the raw genotype count. The
standardized genotype is

$ Z_i = (g_i - 2 p_i) / sqrt(2 p_i (1 - p_i)), $

so that $bb(E)[Z_i] = 0$ and $"Var"(Z_i) = 1$. LD between variants
$i$ and $j$ is measured by the genotypic correlation

$ r_(i j) = "Cor"(g_i, g_j) = bb(E)[Z_i Z_j], $

and we write $R in RR^(m times m)$ for the LD correlation matrix.


  The LD score of SNP $i$ is the sum of squared correlations with all variants:

  $ ell_i = sum_(j = 1)^m r_(i j)^2 = sum_(j = 1)^m bb(E)[Z_i Z_j]^2 =  sum_(j = 1)^m (R circle.small R)_(i j ) =  sum_(j = 1)^m mu^2_(i j ) $


- link LD score to $mu_(i j) = bb(E)[Z_i Z_j]$ and $sigma^2_(i j) = "Var"(Z_i Z_j)$:
  $ ell_i = sum_(j = 1)^m mu_(i j)^2 $
  $ sigma_(i j)^2 = bb(E)[Z_i^2 Z_j^2] - mu_(i j)^2
    quad <==> quad
    bb(E)[Z_i^2 Z_j^2] = sigma_(i j)^2 + mu_(i j)^2 $
  $ ell_i = sum_(j = 1)^m (bb(E)[Z_i^2 Z_j^2] - sigma_(i j)^2) $

$ sum^m_(i j) mu_(i j)^2 = sum^m_(i) ell_i $
$ sum_(i eq.not j)^m mu_(i j)^2 = sum_(i = 1)^m ell_i - m $

$ sum_(i eq.not j)^m sigma_(i j)^2 = ? $





= 2 Different way to build linear operator (remove $circle.small$), can we do it in  $O(n m)$?
\

make it easy, start from unstandardized way:

#align(left)[
$ bold(W) &= p^(-1) bold(H) bold(H)^T \
&= 1/(2p) ((bold(Z) bold(Z)^T) circle.small (bold(Z) bold(Z)^T) - (bold(Z) circle.small bold(Z))(bold(Z) circle.small bold(Z))^T) \
&= 1/(2p) (bold(K_w) circle.small bold(K_w) - bold(D) bold(D)^T) $
]



 $(bold(Z) circle.small bold(Z))(bold(Z)  circle.small bold(Z))^T$ is easy, $bold(Z) circle.small bold(Z)$ takes $O(m n)$ and same space like $bold(Z)$, The linear operator is  $(bold(Z) circle.small bold(Z))(bold(Z)  circle.small bold(Z))^T u $, takes  $O(m n)$. 
 \
 Then we should use way to remove $circle.small $ in $(bold(Z) bold(Z)^T) circle.small (bold(Z) bold(Z)^T)$. Some useful identities:
$ u^T (A circle.small B) v = tr (D_u A D_v B^T) $

$ (A circle.small B) v = "diag" (A D_v B^T) $


$ (A circle.small B) u = E [ (A D_u B^T z) circle.small z ] $
\
The linear operator:
$ ((bold(Z) bold(Z)^T) circle.small (bold(Z) bold(Z)^T)) u
  &= "diag" (bold(Z) bold(Z)^T D_u bold(Z) bold(Z)^T) \
  &= EE [ ((bold(Z) bold(Z)^T) D_u (bold(Z) bold(Z)^T) v) circle.small v ] \
  &approx 1/"Nmc" sum_(i=1)^"Nmc" (bold(Z) bold(Z)^T D_u bold(Z) bold(Z)^T v_i) circle.small v_i \
  &= 1/"Nmc" ( ( bold(Z) ( bold(Z)^T ( u circle.small ( bold(Z) ( bold(Z)^T bold(V) ) ) ) ) ) circle.small bold(V) ) bold(1)_"Nmc" $

where $v$ is a random vector and $V$  is the matix form. This is $O(m n "Nmc")$. Overall,
#align(left)[
$ bold(W) u
  &= 1/(2p) [ ((bold(K_w) circle.small bold(K_w)) u - bold(D) (bold(D)^T u) ] \
  &approx 1/(2p) [ 1/"Nmc" ( ( bold(Z) ( bold(Z)^T ( u circle.small ( bold(Z) ( bold(Z)^T bold(V) ) ) ) ) ) circle.small bold(V) ) bold(1)_"Nmc"
     - bold(D) (bold(D)^T u) ] $
]

where $bold(D) = bold(Z) circle.small bold(Z)$ and $bold(V) = [v_1, ..., v_"Nmc"] in RR^(n times "Nmc")$ has
i.i.d. Rademacher entries. The whole apply costs $O(m n "Nmc")$ against $O(n m^2)$ for the exact route.

== Can we build standardized $W u$?

Yes --- and the standardization does not change the *shape* of the problem at all.
Everything below is verified numerically to machine precision (`verify_std_operator.py`,
`final_formulas.py`).

=== 2.1 The reduction: a rank-$q$ weight factorization gives $q$ unstandardized applies

Write $bold(K)_b := bold(Z) bold(D)_b bold(Z)^T$ for the $b$-reweighted GRM, i.e.
$(bold(K)_b)_(t s) = sum_i b_i Z_(t i) Z_(s i)$. The only $u$-dependent object in the
standardized operator is the quartic term
$ A_(t s) = sum_(i,j) B_(i j) Z_(t i) Z_(t j) Z_(s i) Z_(s j) . $
Suppose the weight matrix admits a rank-$q$ factorization
$ bold(B) = sum_(a=1)^q s_a bold(b)_a bold(c)_a^T , quad s_a in {plus.minus 1} . $
Substituting and splitting the $i$- and $j$-sums,
$ A_(t s) = sum_(a=1)^q s_a
  underbrace((sum_i b_(a i) Z_(t i) Z_(s i)), (bold(K)_(b_a))_(t s))
  underbrace((sum_j c_(a j) Z_(t j) Z_(s j)), (bold(K)_(c_a))_(t s))
  quad ==> quad
  bold(A) = sum_(a=1)^q s_a thin bold(K)_(b_a) circle.small bold(K)_(c_a) . $

*This is the whole answer.* The standardized quartic is a sum of $q$ objects of exactly
the form $bold(K) circle.small bold(K)$ that section 2 already knows how to apply. The
unstandardized case is the special case $q = 1$, $bold(b) = bold(c) = bold(1)$,
$bold(K)_bold(1) = bold(K)_w$. Standardization therefore costs a *factor* $q$, never a
new structure: any method that applies $bold(K) circle.small bold(K)$ in $O(n m)$ applies
the standardized $bold(W)$ in $O(n m q)$.

=== 2.2 $v_R$ and $s_T$ are the same primitive evaluated at $u = bold(1)$

The $u$-independent parts need no separate machinery either. With $bold(G) = bold(Z)^T bold(Z) \/ n$
(so $mu_(i j) = G_(i j)$) and $bold(R) = bold(G) circle.small bold(B)$,
$ (v_R)_t = sum_(i j) mu_(i j) B_(i j) Z_(t i) Z_(t j)
  = 1/n sum_a s_a sum_s (bold(K)_(b_a))_(t s) (bold(K)_(c_a))_(t s)
  = 1/n [ sum_a s_a (bold(K)_(b_a) circle.small bold(K)_(c_a)) bold(1) ]_t , $
$ s_T = sum_(i j) mu_(i j)^2 B_(i j)
  = 1/n^2 sum_a s_a bold(1)^T (bold(K)_(b_a) circle.small bold(K)_(c_a)) bold(1)
  = 1/n bold(1)^T v_R^"full" . $
So $v_R$ is one apply of the *same* operator at $u = bold(1)$, and $s_T$ is one more inner
product. The $O(n m^2)$ setup that `compute_gene_weights` / `setup_pooled` currently pay
for $bold(R)$ and $bold(T)$ disappears entirely.

=== 2.3 The zero diagonal is an exact $O(n m)$ correction

$bold(B), bold(R), bold(T)$ carry a zero diagonal ($i eq.not j$) but the factorization does
not. With $bold(D) = bold(Z) circle.small bold(Z)$ and
$beta_i = B_(i i) = sum_a s_a b_(a i) c_(a i)$, and using $mu_(i i) = 1$ for standardized
columns (so $"diag"(bold(R)) = "diag"(bold(T)) = beta$):
$ bold(A) u = bold(A)^"full" u - bold(D) (beta circle.small (bold(D)^T u)) , quad
  v_R = v_R^"full" - bold(D) beta , quad
  s_T = s_T^"full" - bold(1)^T beta . $
All three corrections are exact and cost $O(n m)$. At $bold(b) = bold(c) = bold(1)$ the
first one is exactly the $- bold(D) (bold(D)^T u)$ of section 2, so the new formula
degenerates correctly.

=== 2.4 Assembled operator

$ bold(W) u = 1/(2 p) [
  sum_a s_a (bold(K)_(b_a) circle.small bold(K)_(c_a)) u
  - bold(D)(beta circle.small (bold(D)^T u))
  - s_u thin v_R + (s_u s_T - u^T v_R) bold(1) ] , quad s_u = bold(1)^T u , $
with each Hadamard product applied by the *symmetrized* estimator
$ (bold(K)_b circle.small bold(K)_c) u approx 1/(2"Nmc") sum_(l=1)^"Nmc" [
  v_l circle.small (bold(K)_b (u circle.small (bold(K)_c v_l)))
  + (bold(K)_c v_l) circle.small (bold(K)_b (v_l circle.small u)) ] . $
Symmetrization is not cosmetic: the one-sided estimator of section 2 is *asymmetric* at
finite Nmc, and both halves reuse the same intermediates, so it is free. Each
$bold(K)_b x = bold(Z)(b circle.small (bold(Z)^T x))$ is two gemms, $O(n m)$; batching over
$a$ makes the apply three $n times m times q$ gemms per probe.

*The probes must be frozen.* Drawing fresh $bold(V)$ per apply makes $hat(bold(W))$ a
different operator on every CG iteration, so CG and Lanczos have no fixed matrix to work on
and lose their convergence guarantees. One probe set must be drawn once and reused for
every apply in the whole REML run.

=== 2.5 The obstruction: the error law is $sqrt(n\/"Nmc")$, with no $m$ in it

Measured relative Frobenius error of $bold(K) circle.small bold(K)$ at $"Nmc" = 256$
(`scaling_and_weightmodel.py`):

#table(columns: 5, align: center,
 [$n$], [$m = 40$], [$m = 80$], [$m = 160$], [$"relFro" dot sqrt("Nmc"\/n)$],
 [200], [0.292], [0.351], [0.256], [$approx 0.33$],
 [400], [0.476], [0.482], [0.392], [$approx 0.36$],
 [800], [0.728], [0.670], [0.625], [$approx 0.38$],
)

The last column is flat in *both* $n$ and $m$:
$ "relFro" approx 0.35 sqrt(n \/ "Nmc") . $
The accuracy is governed by $n$, not by $m$. To hold a fixed relative error one needs
$"Nmc" prop n$, so the true fixed-accuracy cost of the stochastic route is $O(n^2 m)$ ---
*worse* than the exact $O(n m^2)$ whenever $n > m$. The $O(n m "Nmc")$ of section 2 is only
$O(n m)$ if Nmc is allowed to stay constant, which fixes nothing about the error.

The reason is spectral, not algebraic. A rank-Nmc sketch can only represent a matrix whose
effective rank is $lt.tilde$ Nmc, and here
$ "tr"(bold(W))^2 \/ "tr"(bold(W)^2) approx 73 "at" n = 300 , quad 82 "at" n = 400 $
--- $bold(W)$ has effective rank $O(n)$, because $bold(W) approx bold(I)$ plus structure
(this is the same near-identity fact that drives the weak-identifiability behaviour). No
$O(n m)$ sketch can approximate a near-full-rank $n times n$ matrix.

*Consequence for REML.* With frozen probes the sketched $hat(bold(W))$ is indefinite, and
profile REML on it is destroyed (`reml_sensitivity.py`, $n = 400$, $m = 60$, true $h^2 = 0.5$):

#table(columns: 4, align: center,
 [Nmc], [relFro], [$lambda_min (hat(bold(W)))$], [$hat(h)^2$],
 [32], [1.58], [$-13.69$], [1.000],
 [128], [0.73], [$-4.14$], [1.000],
 [512], [0.37], [$-0.86$], [0.308],
 [exact], [0], [$0.000$], [0.514],
)

$lambda_min < 0$ alone breaks the pipeline: $bold(V) = sigma^2_"gxg" hat(bold(W)) + sigma^2_e bold(I)$
is no longer positive definite, so CG diverges and the SLQ nodes leave $(0, infinity)$.

=== 2.6 What does work: $bold(B) = bold(1) bold(1)^T$

The useful question is not how to sketch $bold(B)$ but whether $bold(B)$ needs to be there.
Two facts:

+ *$bold(B)$ has almost no exploitable structure beyond its mean level.* Eigen-truncating
  the exact $bold(B)$ at rank $q$ ($n = 400$, $m = 200$) gives $||bold(W)_q u - bold(W) u|| \/ ||bold(W) u||$
  of $0.067$ at $q = 1$ and $0.047$ at $q = 64$ --- 64$times$ the work for a $1.4 times$
  gain. Nystrom on landmark SNPs (which gets the factors in $O(n m s)$ without ever forming
  $bold(B)$) is *non-monotone* in $s$, because $bold(B)$ is not PSD so $bold(B)_(S S)^+$ is
  ill-conditioned; it never beat $q = 1$.

+ *$q = 1$ is essentially $bold(1) bold(1)^T$.* Setting $bold(B) = bold(1) bold(1)^T$
  outright --- dropping the $1\/sigma^2_(i j)$ scaling, keeping the centring
  $bold(R) = bold(G)$, $bold(T) = bold(G) circle.small bold(G)$ --- gives $0.075$, against
  $0.067$ for the *optimal* rank one.

And the scaling barely moves the estimate. At $n = 500$, $m = 80$, $y$ simulated from the
fully standardized $bold(W)$, 20 replicates (`does_scaling_matter.py`):

#table(columns: 4, align: center,
 [kernel used to fit], [$"corr"$ of off-diag. with $bold(W)_"std"$], [$||dot - bold(W)_"std"||_F \/ ||bold(W)_"std"||_F$], [$hat(h)^2$ (sd)],
 [$bold(W)_"std"$ (full)], [1], [0], [0.490 (0.056)],
 [$bold(W)_"cen"$ ($bold(B) = bold(1)bold(1)^T$)], [0.989], [0.213], [0.458 (0.055)],
 [$bold(W)_"raw"$ (no centring)], [0.846], [0.797], [0.448 (0.063)],
)

The centring-only kernel is $0.99$-correlated with the standardized one and biases $h^2$ by
$0.03$ --- inside half a Monte-Carlo standard deviation. Dropping the centring *as well*
costs another $0.01$ in bias but $0.14$ in correlation, so centring is the part worth
keeping.

So the practical recommendation is $bold(B) = bold(1) bold(1)^T$, for which the $m times m$
weight matrices vanish completely:
$ bold(W) u = 1/(2p) [ (bold(K)_w circle.small bold(K)_w) u - bold(D)(bold(D)^T u)
  - s_u thin v_R + (s_u s_T - u^T v_R) bold(1) ] , $
$ v_R = 1/n (bold(K)_w circle.small bold(K)_w) bold(1) - bold(D) bold(1) , quad
  s_T = 1/n bold(1)^T (v_R + bold(D) bold(1)) - m . $
Verified against the explicitly centred $bold(H) bold(H)^T \/ p$ to $4 times 10^(-16)$.
This is *exactly* the operator of section 2 plus two rank-one corrections, so the
standardized apply costs the same as the unstandardized one --- $q = 1$, no weight storage,
and the `V_list` of $sum_g m_g^2$ doubles in `setup_pooled` is gone.

=== 2.7 Summary

#table(columns: 4, align: (left, center, center, left),
 [object], [exact cost], [$O(n m)$ possible?], [note],
 [$bold(A) u = (bold(B) circle.small bold(M))$ contraction], [$O(n m^2)$], [no], [$sqrt(n\/"Nmc")$ error; needs $"Nmc" prop n$],
 [zero-diagonal correction], [$O(n m)$], [yes], [exact],
 [$v_R$, $s_T$], [$O(n m^2)$ once], [same primitive], [inherits the same error law],
 [weight matrices $bold(B), bold(R), bold(T)$], [$O(n m^2)$, $O(m^2)$ store], [yes, $q = 1$], [set $bold(B) = bold(1)bold(1)^T$; costs $0.03$ in $h^2$],
)

The honest conclusion: the $O(n m)$ *standardization* is solved and free ($q = 1$), but the
$O(n m)$ *operator* is not --- the $O(n m^2)$ in section 1 comes from the Hadamard square
$bold(K)_w circle.small bold(K)_w$, which is present already in the unstandardized case, and
sketching it in the $n$-dimension trades a factor $m\/"Nmc"$ in cost for a factor
$sqrt(n\/"Nmc")$ in error. Use the stochastic form for the *scalar* functionals of section 3
(where the $n$-dimension is never sketched away and the errors are averages, not an
operator); keep the exact Khatri--Rao contraction for anything that goes inside CG or
Lanczos.

= We can do MoM in $O(n m)$

We estimate $tr[W]$, $tr[W^2]$ and $y^T W y$ without forming $K$, using the stochastic trace
estimator
$ hat(tr)[A] = 1/n_(m c) sum_(l=1)^n_(m c) u_l^T A u_l , $
where $n_(m c) <<< min(n, m)$ and each $u_l$ is a random vector with independent, mean-zero,
unit-variance elements. We denote this class by $U$.

*Pair probes.* For each row $t = 1, ..., n$ and $u in U^m$,
$ u^T X_(t.)^T X_(t.) u = 2 sum_(i<j)^m X_(t i) X_(t j) u_i u_j + sum_i^m X_(t i)^2 u_i^2
  = 2 H_(t.) v + X_(t.) X_(t.)^T , $
where $v_k = u_i u_j$ for $k <-> (i < j)$ satisfies $v in U^p$. Hence $z = H v$ is available
elementwise from a single $X u$:
$ z_t = 1/2 ( u^T X_(t.)^T X_(t.) u - X_(t.) X_(t.)^T ) , quad O(n m) . $
Since $EE[v v^T] = I_p$, we have $EE[z z^T] = H H^T = p W$, so one $z$ serves both traces.

== $tr[W]$

$ tr[W] = p^(-1) EE[z^T z]
  = 1/(4p) EE sum_t ( u^T X_(t.)^T X_(t.) u - X_(t.) X_(t.)^T )^2 , $
$ hat(tr)[W] = 1/(4 p n_(m c)) sum_(l=1)^n_(m c) sum_(t=1)^n
  ( tilde(u)_l^T X_(t.)^T X_(t.) tilde(u)_l - X_(t.) X_(t.)^T )^2 , quad tilde(u)_l in U^m . $

== $tr[W^2]$

$ tr[W^2] &= p^(-2) EE[z^T H H^T z]
  = p^(-2) EE sum_(i<j) ( sum_(t=1)^n X_(t i) X_(t j) z_t )^2 \
&= 1/(2p^2) ( sum_(i,j) sum_(t,l) X_(l i) X_(l j) X_(t i) X_(t j) z_t z_l
   - sum_(i=1)^m sum_(t,l) X_(l i)^2 X_(t i)^2 z_t z_l ) \
&= 1/(2p^2) tr[K D_z K D_z] - 1/(2p^2) z^T (X compose X)(X compose X)^T z , $
where $D_z$ is the diagonal matrix generated by $z$. Averaging pair probes $u_l in U^m$,
$ hat(z)_t = 1/(2 sqrt(n_(m c))) sum_(l=1)^n_(m c)
  ( u_l^T X_(t.)^T X_(t.) u_l - X_(t.) X_(t.)^T ) , $
the $sqrt(n_(m c))$ keeping $EE[hat(z) hat(z)^T] = p W$. Given $hat(z)$,
$ hat(tr)[W^2] = 1/(2 p^2 n_(m c)) sum_(l=1)^n_(m c)
  tilde(u)_l^T X X^T D_hat(z) X X^T D_hat(z) tilde(u)_l
  - 1/(2p^2) hat(z)^T (X compose X)(X compose X)^T hat(z) . $



Every term above is applied right to left, so $K = X X^T$ is never formed and each probe costs
four matrix--vector products with $X$, i.e. $O(n m)$.


== how to solve that

The cross terms cannot be dropped. Each $W_g$ is positive semidefinite, so
$tr[W_g W_h] = angle(W_g, W_h)_F >= 0$, and since $EE[W_g] = I_n$ for standardized
genotypes, independence between two genes gives $EE tr[W_g W_h] = tr[I_n] = n$ rather than zero.
Their total therefore grows like $G(G-1) n$ against $G n$ for the diagonal part, so they dominate
$tr[W^2]$ for any realistic $G$.

They can, however, be made implicit. Draw one pair probe $u_g in U^(m_g)$ per gene,
*independently across genes*, and form
$ z_(g,t) = 1/(2 sqrt(p_g)) ( u_g^T X_(g,t.)^T X_(g,t.) u_g - X_(g,t.) X_(g,t.)^T ) ,
  quad z = sum_(g=1)^G z_g . $
Then $EE[z_g z_g^T] = p_g^(-1) H_g H_g^T = W_g$, and independence kills every
$EE[z_g z_h^T]$, so
$ EE[z z^T] = sum_(g=1)^G W_g = W
  quad ==> quad tr[W^2] = EE[z^T W z] = EE[ sum_(g=1)^G z^T W_g z ] . $
The cross terms are carried by the covariance of $z$ and never have to be enumerated.

*The same $z$ must enter every term.* Expanding $z^T W_g z = sum_(h, h') z_h^T W_g z_(h')$ shows
where they live; replacing $z$ by the gene's own $z_g$ in the $g$-th term would return
$sum_g tr[W_g^2]$ and discard exactly the part that dominates.

Each summand is then the single-gene quadratic form, evaluated with an inner probe
$c_l in U^(m_g)$ drawn separately per gene:
$ z^T W_g z = 1/(2 p_g) ( tr[M_(g,z)^2] - norm(D_g^T z)^2 ) , quad M_(g,z) = X_g^T D_z X_g , $
$ tr[M_(g,z)^2] approx 1/n_"in" sum_(l=1)^n_"in" norm( X_g^T ( z compose (X_g c_l) ) )^2 . $

Since $EE[z z^T] = W$ already holds for a single draw, enlarging the probe count inside one $z$
does not reduce the variance; the estimator must be averaged over $R$ independent draws,
$ hat(tr)[W^2] = 1/R sum_(r=1)^R sum_(g=1)^G
  ( 1/(2 p_g n_"in") sum_(l=1)^n_"in" norm( X_g^T ( z^((r)) compose (X_g c_l) ) )^2
    - 1/(2 p_g) norm(D_g^T z^((r)))^2 ) , $
converging at $O(R^(-1/2))$. Building one $z$ costs $O(n M)$ and contracting it costs
$O(n M n_"in")$, so the total is $O(n M R thin n_"in")$ --- linear in the number of variants and
with no $binom(G,2)$ term anywhere.

== how to solve that

Draw an independent pair probe $u_g in U^(m_g)$ per gene and set
$ z_(g,t) = 1/(2 sqrt(p_g)) ( u_g^T X_(g, t.)^T X_(g, t.) u_g - X_(g, t.) X_(g, t.)^T ) ,
  quad z = sum_(g=1)^G z_g , $
so that $EE[z_g z_g^T] = p_g^(-1) H_g H_g^T = W_g$ and, by independence across genes,
$ EE[z z^T] = sum_(g=1)^G W_g = W
  quad ==> quad tr[W^2] = EE[z^T W z] = EE[ sum_(g=1)^G z^T W_g z ] . $
The cross terms never appear: they are carried by the covariance of $z$, and the quadratic form is
evaluated gene by gene with the single-gene machinery,
$ z^T W_g z = 1/(2 p_g) ( tr[M_(g,z)^2] - norm(D_g^T z)^2 ) , quad M_(g,z) = X_g^T D_z X_g ,
  quad tr[M_(g,z)^2] approx 1/n_"in" sum_(l=1)^n_"in" norm(X_g^T (z compose (X_g c_l)))^2 , $
with $c_l in U^(m_g)$. Independence of the $u_g$ is essential: reusing one probe across genes
would make $EE[z_g z_h^T] eq.not 0$ and bias the estimate.

== Averaging

$EE[z z^T] = W$ holds for a single draw, so raising the number of pair probes inside one $z$ does
not reduce the variance of $tr[W^2]$; the estimator must be averaged over $R$ independent draws
of $z$:
$ hat(tr)[W^2] = 1/R sum_(r=1)^R sum_(g=1)^G
  ( 1/(2 p_g n_"in") sum_(l=1)^n_"in" norm(X_g^T (z^((r)) compose (X_g c_l)))^2
    - 1/(2 p_g) norm(D_g^T z^((r)))^2 ) , $
which converges at $O(R^(-1/2))$. Each draw costs $O(n M)$ to build and $O(n M n_"in")$ to
contract, so the total is $O(n M R thin n_"in")$ with no dependence on $G$ beyond $M$.

== Per-gene variance components

If instead each gene carries its own $sigma_g^2$, the moment equations need the full
$G times G$ matrix $tr[W_g W_h]$, which the shared probe cannot deliver --- it returns only the
sum. Drawing $G$ separate probes and using $tr[W_g W_h] = EE[z_g^T W_h z_g]$ restores them at
$binom(G, 2)$ quadratic forms, which is the practical reason a single pooled component is usually
estimated instead.





== 2.3 SVD
low rank approximate factor A ~~ La Ra, B ~~ Lb Rb.  then (A*B)u approximately expands as matvecs and vec*vec hadamard products. for r unique svs of A, k unique svs of B, define w_rk := l_Ar *l_Bk,  z_rk := r_Ar *r_Bk where R, L are as above
(A*B)u ~~ sum_r,k  w_rk z_rk^T u
svs = singular values

FASTGWA
