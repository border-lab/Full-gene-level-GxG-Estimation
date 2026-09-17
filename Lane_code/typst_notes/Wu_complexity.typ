#set page(margin: 1in)
#set text(size: 11pt)
#set par(justify: true)
#set heading(numbering: "1.1")

#align(center)[
  #text(size: 16pt, weight: "bold")[
    Time and Space Complexity of the Matrix-Free $W u$ Apply
  ]
  #v(4pt)
  
]

*Input*:
$Z in RR^(n times m)$ is the genotype
($n$ individuals, $m$ SNPs), $u in RR^n$ is the vector to which $W$ is applied,
and $p = binom(m, 2) = m(m-1)\/2$ is the number of within-set SNP pairs. The
apply is analysed one vector at a time; a block of $c$ probe vectors is written
$U in RR^(n times c)$ and costs $c$ times the single-vector figure throughout.

= The Math

*Compute three weight matrices.*
To express these summations in matrix form, we defined three symmetric $m times m$ matrices with zero diagonal ($a != b$):
$ V_(a b) = 1 / sigma^2_(a b), quad R_(a b) = mu_(a b) / sigma^2_(a b), quad T_(a b) = mu_(a b)^2 / sigma^2_(a b) $

where $mu_(a b) = bb(E)[Z_a Z_b]$ and $sigma^2_(a b) = "Var"(Z_a Z_b)$ are computed directly from the individual-level genotype data as the empirical mean and variance of $Z_(dot.c a) circle.small Z_(dot.c b)$ across all $n$ individuals. The time complexity is $O(n m^2)$ and the space complexity is $O(n m + m^2)$ (the $n times m$ element-wise square $Z circle.small Z$ is needed for $bb(E)[(Z_a Z_b)^2]$).\

*Formula for the implicit $W u$.*
Let $bold(1)$ be the all-ones vector, $s_u = bold(1)^T u = sum_i u_i$ the sum of
$u$, and
$ M = Z^T "diag"(u) Z in RR^(m times m), quad
  M_(a b) = sum_i u_i Z_(i a) Z_(i b) , $
the single $u$-dependent contraction ($u$ enters *only* here, through
$"diag"(u)$). Expanding
$W u = 1/p sum_(a<b) h_(a b) (h_(a b)^T u)$ with the standardized interaction
column $h_(a b) = (d_(a b) - mu_(a b) bold(1)) \/ sigma_(a b)$,
$d_(a b) = Z_(dot.c a) circle.small Z_(dot.c b)$, and collecting each of the four
resulting pieces under one weight matrix gives the form below. Each bracketed
term is a sum over *ordered* pairs $a != b$ --- that is what the zero diagonal of
$V, R, T$ buys us --- so every unordered pair is counted twice and the $1 \/ p$
becomes $1 \/ (2 p)$:
$ W u = 1/(2 p) [
     (Z circle.small (Z (V circle.small M))) bold(1)
  - s_u (Z circle.small (Z R)) bold(1)
  -  (bold(1)^T (R circle.small M) bold(1)) bold(1)
  + s_u (bold(1)^T T bold(1)) bold(1)
  ] , $
where $A bold(1)$ is the row-sum of $A$ (a vector in $RR^n$) and
$bold(1)^T A bold(1)$ the sum of all entries of $A$ (a scalar). Every term is
built from $Z$, $u$ and the weight matrices $V, R, T$ alone --- no $n times n$
GRM $W$ and no $n times p$ interaction matrix $H$ is ever formed.

*Degenerate pairs.* A pair whose product column $d_(a b)$ is (near-)constant ---
a near-monomorphic SNP, or a perfectly co-inherited pair --- has
$sigma^2_(a b) approx 0$, so $h_(a b)$ is undefined and $V_(a b) = 1 \/ sigma^2_(a b)$
would overflow. The convention is to *mask* such pairs: if
$sigma^2_(a b) <= epsilon$ (the pipeline uses $epsilon = 10^(-20)$) then
$ V_(a b) = R_(a b) = T_(a b) = 0 , $
which is exactly $h_(a b) := bold(0)$ --- the pair contributes nothing to $W$.
It is still counted in the normalizer $p$, matching the explicit construction,
which zeroes the same columns of $H$ and divides by $p$ regardless.



= Time and space complexity for $W u$

We read the cost straight off the formula, evaluating $W u$ for a single vector
$u in RR^n$ one operation at a time. Write the four brackets as
$ t_1 = (Z circle.small (Z (V circle.small M))) bold(1), quad
  t_2 = s_u (Z circle.small (Z R)) bold(1), quad
  t_3 = (bold(1)^T (R circle.small M) bold(1)) bold(1), quad
  t_4 = s_u (bold(1)^T T bold(1)) bold(1), $
so that $W u = (t_1 - t_2 - t_3 + t_4) \/ (2 p)$. 
- $A bold(1)$ is the row-sum of a matrix $A$
- $bold(1)^T A bold(1)$ its total sum. 


#table(
  columns: (auto, 1.6fr, auto, auto, auto),
  align: (left, left, center, center, center),
  inset: 6pt,
  stroke: 0.5pt + luma(180),
  table.header(
    [*Term*], [*Calculation*], [*Shape*], [*Time*], [*Space*],
  ),

  [shared], [$s_u = bold(1)^T u$], [scalar], [$O(n)$], [$O(1)$],
  [shared], [$"diag"(u) Z$ ], [$n times m$], [$O(n m)$], [$O(n m)$],
  [shared], [$M = Z^T ("diag"(u) Z)$], [$m times m$], [$O(n m^2)$], [$O(m^2)$],

  [$t_1$], [$V circle.small M$ ], [$m times m$], [$O(m^2)$], [$O(m^2)$],
  [], [$Z (V circle.small M)$], [$n times m$], [$O(n m^2)$], [$O(n m)$],
  [], [$(Z circle.small (Z(V circle.small M))) bold(1)$ ], [$n times 1$], [$O(n m)$], [$O(n)$],

  [$t_2$], [$Z R$], [$n times m$], [$O(n m^2)$], [$O(n m)$],
  [], [$s_u (Z circle.small (Z R)) bold(1)$ ], [$n times 1$], [$O(n m)$], [$O(n)$],

  [$t_3$], [$bold(1)^T (R circle.small M) bold(1)$ ], [scalar], [$O(m^2)$], [$O(1)$],

  [$t_4$], [$s_u thin bold(1)^T T bold(1)$ ], [scalar], [$O(m^2)$], [$O(1)$],

  [combine], [$(t_1 - t_2 - t_3 + t_4) \/ (2 p)$], [$n$], [$O(n)$], [$O(n)$],
)


$ "Time" = O(n m^2), quad "Space" = O(n m + m^2) . $



= How to reduce the constant factor

Three $n m^2$ products dominate the time:

#table(
  columns: (auto, auto, auto),
  align: (left, left, center),
  inset: 6pt,
  stroke: 0.5pt + luma(180),
  table.header([*Product*], [*Part*], [*Cost*]),
  [$M = Z^T "diag"(u) Z$], [$t_1, t_3$], [$1 times n m^2$],
  [$Z (V circle.small M)$], [$t_1$], [$1 times n m^2$],
  [$Z R$], [$t_2$], [$1 times n m^2$],
)
The linear operator cost $approx 3 n m^2$ per apply.

*Precompute the elements independent of $u$ *

In $t_2$, the $(Z circle.small (Z R)) bold(1)$ can be computed once, which is a vector, and be denoted as the $v_R$, so does $bold(1)^T T bold(1)$, can be precomputed as a scalar $s_T$:

$ v_R = (Z circle.small (Z R)) bold(1) in RR^n, quad
  s_T = bold(1)^T T bold(1) in RR . $
Then $t_2 = s_u thin v_R$ ($O(n)$) and $t_4 = s_u thin s_T$ ($O(1)$).

The vector $v_R$ is worth a second look: its $i$-th entry is
$ (v_R)_i = sum_(a, b) Z_(i a) R_(a b) Z_(i b) = z_i^T R z_i = (Z R Z^T)_(i i) , $
writing $z_i^T$ for row $i$ of $Z$. So $v_R = "dg"(Z R Z^T)$, where $"dg"(dot)$
extracts a diagonal into a vector (as opposed to $"diag"(dot)$, which builds a
matrix from one).

Also $t_3$ can be represented as $(u^T v_R) bold(1)$, as:

$ bold(1)^T (R circle.small M) bold(1) = lr(⟨R, M⟩) = "tr"(R M)
  = "tr"(R Z^T "diag"(u) Z) = "tr"("diag"(u) Z R Z^T)
  = u^T "dg"(Z R Z^T) = u^T v_R , $

using $R = R^T$ in the second equality and cyclicity of the trace in the fourth.
So $t_3 = (u^T v_R) bold(1)$ reuses the *same* precomputed $v_R$ at $O(n)$ --- no
new object is needed. The product $Z R$ leaves the apply entirely, and with it
one of the three $n m^2$ passes.



= The new formula for $W u$

The precomputation splits the work cleanly into a *setup* phase that runs once
per gene and depends on $Z$ alone, and an *apply* phase that runs once per
vector and is the only place $u$ appears.

*Setup (once per gene).* From $Z$:
+ the weight matrices $V, R, T in RR^(m times m)$, symmetric with zero diagonal,
  from the empirical $mu_(a b)$ and $sigma^2_(a b)$ --- $O(n m^2)$;
+ the vector $v_R = (Z circle.small (Z R)) bold(1) $ ---
  $O(n m^2)$;
+ the scalar $s_T = bold(1)^T T bold(1) in RR$ --- $O(m^2)$.

Only $Z$, $V$, $v_R$ and $s_T$ stay; $R$, $T$ and the $n times m$ intermediate
$Z R$ can all be discarded.

*Apply (per vector $u$).* Exactly one $u$-dependent object is built:
$ s_u = bold(1)^T u, quad M = Z^T "diag"(u) Z $
and the result is


$ W u = 1/(2 p) [
  (Z circle.small (Z (V circle.small M))) bold(1) - s_u thin v_R + (s_u s_T - u^T v_R) bold(1)
  ]  . $




*Cost.*

#table(
  columns: (auto, auto, auto),
  align: (left, center, center),
  inset: 6pt,
  stroke: 0.5pt + luma(180),
  table.header([*Phase*], [*Time*], [*Space*]),
  [setup (once per gene)], [$O(n m^2)$], [$O(n m + m^2)$],
  [apply, one vector], [$2 n m^2$ (build $M$, then $Z (V circle.small M)$)],
    [$O(n m)$ scratch],
  [apply, block of $c$], [$2 c thin n m^2$], [$O(c m^2 + n m)$ scratch],
)



