#set page(margin: 1in)
#set text(size: 11pt)
#set par(justify: true)
#set heading(numbering: "1.1")

#align(center)[
  #text(size: 15pt, weight: "bold")[
    Why the four variance components must be modelled jointly
  ]
  #v(3pt)
  #text(size: 10pt)[Additive + dominance + additive-by-additive epistasis + residual
  (after Hivert et al., #emph[AJHG] 2021)]
]

= The model

One linear mixed model fits all genetic effects as random effects, giving a
phenotypic covariance that is a sum of four terms:

$ V = "Var"(y) = sigma_a^2 K_a + sigma_d^2 K_d + sigma_"gxg"^2 W + sigma_e^2 I, $

with $K_a$ the additive GRM, $K_d$ the dominance GRM, $W$ the additive-by-additive
(epistasis) GRM, and $I$ the residual. REML estimates
$theta = (sigma_a^2, sigma_d^2, sigma_"gxg"^2, sigma_e^2)$. The question is why all
four enter *one* model rather than being fitted one component at a time.

= Why joint, not separate

*1. Omitted components bias the ones you keep.* The kernels are not orthogonal in
finite data. The epistasis GRM is literally built from the additive one,
$W prop K_a compose K_a$ (Hadamard square), so additive and epistatic relatedness
are correlated. Fit $K_a$ alone and the dominance and epistatic variance have
nowhere to go but into $hat(sigma)_a^2$ -- it inflates. Only the joint model gives
each source its own kernel, so each component is partitioned to where it belongs.

*2. Orthogonality is a design goal, not a guarantee.* The GCTA dominance coding
${0,1,2} |-> {-p\/q, 1, -q\/p}$ is chosen precisely so $K_d$ is orthogonal to
$K_a$ *under Hardy--Weinberg equilibrium* (the naive $0,1,0$ coding gives
$"Cov"(x_A, x_D) = 2p(1-2p) != 0$). Likewise $"Cov"(K_(a,i j), K_("gxg",i j)) =
E[K_(a,i j)^3] = 0$ only in expectation for unrelated, outbred samples;
dominance--epistasis orthogonality holds only empirically. Real data depart from
these ideals, and the joint fit absorbs the residual non-orthogonality that a
sequence of separate fits cannot.

*3. The estimates are statistically coupled.* The sampling variance of one
component depends on the others through their covariance, so a per-component
standard error computed in isolation is simply wrong. Epistatic variance carries
a very large sampling variance and is nearly collinear with the additive (and
residual) terms; it can only be *separated* from them -- and its uncertainty
correctly stated -- when they sit in the model together. This is exactly why the
paper needs $N tilde 10^5$--$10^6$ to pin down $sigma_"gxg"^2$.

*4. The quantities of interest are ratios over the total.* Every reported number
-- $h_"SNP"^2 = sigma_a^2 \/ (sigma_a^2 + sigma_d^2 + sigma_"gxg"^2 + sigma_e^2)$,
the dominance and epistatic fractions, and broad-sense $H^2$ -- is normalized by
the *sum* of all four. Claiming "variance is predominantly additive" is only
meaningful once dominance and epistasis are estimated alongside it, on the same
scale, in the same model.

= Takeaway

The four components are not independent nuisances to be peeled off one by one;
their kernels overlap, their estimates covary, and their meaning is defined
relative to the whole. Joint REML is what makes the additive estimate unbiased,
the epistatic estimate identifiable, and every heritability fraction
interpretable.
