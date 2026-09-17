== Simulation pipeline

Each experiment is a Monte-Carlo study of the *pooled within-gene pairwise-epistasis*
model. The phenotype is drawn from the two-component covariance

$ y = g_(g times g) + e, quad
  V = "Var"(y) = sigma^2_(g times g) W + sigma_e^2 I, quad
  W = 1/P sum_(g=1)^G H_g H_g^T, quad P = sum_(g=1)^G binom(m_g, 2), $

where the $m$ SNPs are cut into $G$ *contiguous* genes (equal, or in the per-gene
proportions given by `RATIO`), each within-gene SNP pair $(a, b)$ contributes the
column-standardized product $H_(a b)$, and the single $1\/P$ division weights every
pair equally --- so a larger gene carries proportionally more heritability
($prop m_g^2$). This is the pooled kernel of `generative_model.typ`; the whole run is
parameterized by $(n, m, G, sigma^2_(g times g), sigma_e^2, "mode", "RATIO")$ plus a
replicate count.

The pipeline `MCREML_pipeline.sh` submits *one SLURM job chain of four steps*, each
waiting on the previous (`--dependency=afterok`). The design decision that shapes
everything is that $W$ is *deterministic* given the genotype and the gene split, so it
is built and Cholesky-factored *once* up front, cached to disk, and thereafter merely
*loaded* --- estimation never touches the genotype.

+ *Step 1 --- Cholesky (single job, `Simulate_Cholesky.py`).* Reads the genotype
  ($n times m$), column-standardizes it, and splits it into the $G$ contiguous gene
  blocks. It builds the pooled kernel $W$ in pair-batches, forms the lower Cholesky
  factor $L_(g times g)$ with $L_(g times g) L_(g times g)^T = sigma^2_(g times g) W$
  (a tiny jitter added for positive-definiteness), and *caches* both $W$ and
  $L_(g times g)$ keyed by a `split` tag (e.g. `G5_r0.1-0.2-0.3-0.2-0.2`). When an
  estimation subset is requested it builds a *second* kernel $W_"est"$ in the same
  pass, pooled over the selected genes only, and records the subset's share of
  within-gene pairs $P_"est"\/P_"full"$ to `info/`.

+ *Step 2 --- Phenotype (array job, `Simulate_Phenotype.py`).* Each array task is one
  replicate. It loads the cached $L_(g times g)$ and draws
  $g_(g times g) = L_(g times g) u_1 tilde N(0, sigma^2_(g times g) W)$ with
  $e = sqrt(sigma_e^2) u_2$, then *rescales each component to its exact target
  variance* (removing sampling error) and mean-centres $y = g_(g times g) + e$ so the
  fitted model carries no fixed effect. The phenotype depends only on the *full*
  kernel, so a single phenotype set is reused across every estimation subset at the
  same split.

+ *Step 3 --- MC-AI-REML (array job, `Simulate_MCREML.py`).* Each task loads the
  pre-computed dense kernel ($W$ for the well-specified fit, or $W_"est"$ for a subset
  fit) and the matching replicate's $y$, then runs Monte-Carlo average-information REML
  on $V = sigma^2_(g times g) W + sigma_e^2 I$: conjugate-gradient solves against $V$,
  a Hutchinson stochastic trace over `NMC` Rademacher probes, and a damped, bounded
  AI-Newton step (Levenberg--Marquardt ridge + trust region + box clamp, needed because
  $W$ is often near-collinear with $I$). The kernel enters only as the dense mat-vec
  $W B$ ($O(n^2 c)$ per CG iteration) --- no genotype is needed. The probe seed is fixed
  to the replicate index, and each task writes its $(hat(sigma)^2_(g times g),
  hat(sigma)^2_e)$ estimate and its wall-clock time.

+ *Step 4 --- Combine and clean up (single job).* Concatenates the per-replicate
  estimates into one result file, averages the timing records, then deletes the
  intermediate result/time directories and the large cached $L_(g times g)$ and
  phenotype set. `calc_stats.py` reports the across-replicate mean, median, SD and a
  95% CI for each component.

*Well-specified vs. misspecified fits.* The phenotype is *always* simulated from the
full $G$-gene kernel, but the `ESTIMATE` setting decouples estimation from it. With
`ESTIMATE=""` the REML fit uses the full $W$ (correctly specified). With a strict
subset (e.g. `ESTIMATE="1,2"`) the fit uses $W_"est"$ pooled over those genes only,
while the truth stays the full panel --- a deliberately misspecified fit. Because
$W_"est"$ is normalized by its *own* pair total, the estimate is attenuated toward
$sigma^2_(g times g) times P_"est"\/P_"full"$, the subset's share of within-gene pairs;
multiplying $hat(sigma)^2_(g times g)$ by $P_"full"\/P_"est"$ rescales it to the full
panel. The next sections read this attenuation off the simulations, and show how LD
between kept and dropped genes inflates the subset estimate above the pair-count target.

== Simulate using all genes and estimate all

Well-specified baseline: the phenotype is drawn from the pooled kernel built over all
$G = 5$ genes.
Genotypes are real, contiguous chr1 SNPs ($m = 10000$),
split into five genes of unequal size in the proportions
$ 0.1 : 0.2 : 0.4 : 0.2 : 0.1$., the true $sigma^2_(g times g) $= 0.2

#figure(
image("chr1_10ksnp_s2gxg0.2_s2e0.8_G5_r0.1-0.2-0.4-0.2-0.1_est1-2-3-4-5_fixed_m10000_gxg_estimate.pdf", width: 100%),caption: [Well-specified case: phenotype simulated from all five genes and estimated with the same five-gene kernel. The estimate $hat(sigma)^2_(g times g)$ is unbiased at every sample size (means $0.200$--$0.202$, all `ns`) and its SD contracts from $0.070$ at $n = 1000$ to $0.011$ at $n = 16000$. Dashed line: truth $sigma^2_(g times g) = 0.2$.]
) <exactMoM>


== Simulate using all genes but estimate using a subset of genes

Here the phenotype is still generated from all genes, but the estimation kernel is
restricted to a *subset* of them.

=== Random case

Independently simulated SNPs ($m = 1000$), $G = 2$ equal-size genes ($r = 0.5 : 0.5$),
with the estimation kernel using gene 2 only.

#figure(
image("RandomSNP_s2gxg0.2_s2e0.8_G2_r0.5-0.5_est2_fixed_m1000_gxg_estimate.pdf", width: 100%),caption: [Random (no-LD) case: two equal-size genes simulated, estimated from gene 2 only. The estimate is significantly biased downward at every $n$ (`***`) and converges to $approx 0.10$ (means $0.139, 0.115, 0.102, 0.099, 0.101$) --- exactly half the truth, the share of pairs carried by the single retained gene. Dashed line: target $0.1$.]
)

As we only estimate the $h^2$ from gene 2, the target is its share of the within-gene pairs
times the truth. The two genes are equal-size ($m_g = 500$ each), so they carry the same
pair count $binom(500, 2) = 124750$ and gene 2's share is exactly $124750 \/ 249500 = 0.5$;
the target is therefore $0.5 times 0.2 = 0.1$ (dashed line). Because the SNPs are independent
(no LD), the dropped gene 1 is invisible to the estimator, so the estimate cannot borrow any
of its variance and converges cleanly to the target: means $0.139, 0.115, 0.102, 0.099,
0.101$ across $n = 1000 dots 16000$, settling at $approx 0.10$. The estimate is significantly
below the truth $0.2$ at every $n$ (`***`), and the mild upward inflation at small $n$
($0.139$ at $n = 1000$) reflects the weak identifiability of the epistatic component in the
small-sample regime and disappears as $n$ grows.

=== Contiguous (LD) case

Real, contiguous chr1 SNPs ($m = 10000$), the same $G = 5$ genes with
$r = 0.1 : 0.2 : 0.4 : 0.2 : 0.1$, but the estimation kernel uses genes 1, 3, 5 only
(`est1-3-5`), dropping genes 2 and 4.

//calculate the pair-wise pair in gene 1,3,5 and devide tby the whole paris in 5 genes and multiply by 0.2 is the target h2(dash line). but with the LD, it infalted

Because the pooled kernel weights each within-gene SNP *pair* equally, the epistatic
variance carried by a subset of genes is proportional to that subset's share of the
within-gene pairs. With $m = 10000$ and $r = 0.1 : 0.2 : 0.4 : 0.2 : 0.1$ the gene sizes are
$m_g = 1000, 2000, 4000, 2000, 1000$, giving pair counts
$binom(m_g, 2) = 499500, thick 1999000, thick 7998000, thick 1999000, thick 499500$ and a
total of $P = 12995000$. The retained genes 1, 3, 5 hold

$ 499500 + 7998000 + 499500 = 8997000 quad "pairs," quad
  8997000 \/ 12995000 approx 0.692, $

so the target for the subset estimate is $0.692 times 0.2 approx 0.138$ (dashed line). The
estimate instead converges to $approx 0.17$ (means $0.170, 0.167, 0.168, 0.169, 0.172$),
*inflated above* the pair-count target: linkage disequilibrium between the retained genes
$(1, 3, 5)$ and the dropped genes $(2, 4)$ lets the subset kernel tag part of the omitted
epistatic variance, so it recovers more than the pairs it explicitly contains. The bias is
significant at every $n$ (`***`) and the SD shrinks from $0.065$ to $0.011$, confirming the
inflation is a stable feature of the LD structure rather than a small-sample artefact.


#figure(
image("chr1_10ksnp_s2gxg0.2_s2e0.8_G5_r0.1-0.2-0.4-0.2-0.1_est1-3-5_fixed_m10000_gxg_estimate.pdf", width: 100%),caption: [Contiguous (LD) case: five genes simulated, estimated from genes 1, 3, 5 only (genes 2, 4 dropped). The estimate converges to $approx 0.17$, above the pair-count target $0.138$ (dashed line); LD between the kept and dropped genes inflates the subset estimate. Significant at every $n$ (`***`).]
)

