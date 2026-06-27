# Verified baseline ingredients — REML within-gene cis-AxA note

Run date: 2026-06-23. Threads capped (`OMP_NUM_THREADS=6`). Re-verified in-repo, not transcribed.

## Task 1 — ingredients

### within_gene_pooling.py (G=25, m=8)
```
M_e^gene   =    8.4   (predict ~m = 8)
M_e^within =   42.1   (predict sqrt(G)*Me_gene = 42.1)   # EXACT
M_e^genome =  200.3   (predict G*Me_gene = 210.5)        # ~4.8% low
gain genome/within = 4.76x   (predict sqrt(G) = 5.00)
```
Note: the within↔gene √G law is exact; the genome = G·gene identity is approximate because the
genome kernel's between-gene AxA pairs make it marginally more concentrated than G independent
genes (small at G=25). The note leans on the exact within↔gene law.

### effective_markers.py
All `M_e` predictions match empirical (`1/M_e=1/n+1/M+ρ̄²`); `r4(Σ)=r4(R)` identity holds;
floors at h2=0.9 match (MoM `√(2/r4)`, REML `√(2/r⋆)`); feature LD decorrelation confirmed.

### check_knobs.py — the constant c (SE_AxA = M_e/(cN))
```
geno     n     m    M_e   SEaxa pred(c=1)  SEaxa obs(full-info)   c = pred/obs
indep   2500  1500  1501      0.6005            0.3036             1.98
LD .6   2500  1500  1099      0.4396            0.2314             1.90
LD .9   2500  1500   368      0.1471            0.1060             1.39
```
=> **c ∈ [1.4, 2.0]**, mid-value ~1.6. Required N ∝ 1/c, so this is a ±~20% band on the headline N.
Additive SE prediction matches observed exactly (`√(2 M_e)/N`).

### verify_all.py
Reproduced everything, exit 0, no drift (centering bias contig m=80: uncentered 0.50 bias -0.40,
corrected 0.88; Hivert HE=0.291 REMLjoint=0.178 meta=0.250; additive-vs-AxA ratio √(Me/2)=192).

## Task 2 — cost of partialling out additive (partialling_additive.py)

**CORRECTION to the design's §3.** The brainstorm claimed partialling out additive is "nearly free"
because `tr(AW) ≈ Σ A_ij³ ≈ 0`. **That is wrong for real genotypes.** Penalty = Var_3comp/Var_2comp
of σ̂²_AA in the detection limit, well-predicted by `1/(1-r²)`, r = off-diag corr(A,W):

```
(A) WORST CASE  cis-additive vs cis-AxA (same SNPs):
    contiguous (high-LD gene)  v3/v2 = 1.38   corr_od +0.53   mean A_ij^3 = +0.13
    random     (dispersed)     v3/v2 = 1.13   corr_od +0.39   mean A_ij^3 = +4.9e-4
(B) Gaussian AR(1) genes: v3/v2 = 1.000 for ALL rho up to 0.9, ALL G  -> the penalty is a
    GENOTYPE-SKEW effect: 0 for symmetric features, >1 only for real (skewed) 0/1/2 genotypes
    under LD. Also confirms the penalty is G-INDEPENDENT (pooled = per-gene).
(C) DILUTION (real data): W = high-LD cis gene AxA; A = 150 cis SNPs + T trans SNPs:
    cis frac 1.00 -> 1.38 ; 0.50 -> 1.37 ; 0.25 -> 1.34 ; 0.15 -> 1.31
    Barely dilutes — high-LD cis SNPs dominate the GRM off-diagonal variance even at 15% count.
```

**Corrected §3 result:** jointly fitting the genome-wide additive GRM inflates the cis-AxA
detection **variance by ~1.1× (low-LD genes) to ~1.4× (high-LD genes)** — SE by √ of that,
≈1.05–1.18× — a genotype-skewness effect, G-independent, robust to trans dilution. NOT free.
Feeds a ~6–18% inflation of the required N in §6.

## Task 3 — pooled kernel vs per-gene aggregation (pooled_vs_pergene.py)

**CORRECTION to the design's §5** (which predicted pooled would dilute/bias under cross-gene sparsity).
```
INDEPENDENT genes (G=40, n=2500, true h2=0.05):
  signal           pool(mean,SD)     per-gene(mean,SD)
  uniform          0.0494, 0.0180    0.0494, 0.0180     <- IDENTICAL
  sparse(5 genes)  0.0506, 0.0200    0.0506, 0.0200     <- IDENTICAL
  sparse(1 gene)   0.0499, 0.0259    0.0499, 0.0259     <- IDENTICAL
BETWEEN-GENE LD (contiguous block split into 40 adjacent windows, cross/within overlap=5.88):
  pooled  : mean 0.0504  SD 0.033   (unbiased, robust)
  per-gene: mean 0.384   SD 0.242   (biased UP 7.7x: naive sum double-counts shared signal)
```
**Corrected §5 result:** for INDEPENDENT genes the pooled-HE and per-gene-sum-HE are the *same*
estimator (cross-kernel off-diag products vanish -> den_W = δ/G -> he_W = Σ he_g), identical bias
and SD in every regime (uniform or sparse). They diverge only under BETWEEN-gene LD, where the
naive per-gene SUM double-counts and is biased upward while the single pooled kernel stays
unbiased. So the pooled kernel is the ROBUST choice; per-gene aggregation needs LD-separated/pruned
genes (or a joint G-component fit). The remaining real difference is computational (per-gene = G
small parallel problems; pooled = one matvec-based fit).

### Task 3b — direct REML check + HE/REML correction (pooled_vs_pergene_reml.py)
CORRECTION: an earlier draft said "detection-regime HE equals REML" -- WRONG. HE and REML are
different estimators. What is true: full-information MoM and REML share efficiency only as h2->0; the
off-diagonal HE of eq.(10) is the LESS-efficient variant (discards the AxA diagonal), used in the note
only as an algebraically transparent proxy. The pooled-vs-per-gene conclusions are verified for ACTUAL
REML directly:
```
(i)  INDEPENDENT genes (REML): pooled mean 0.048 SD 0.021 | per-gene-sum mean 0.053 SD 0.021  -> agree
(ii) BETWEEN-gene LD (REML, 8 windows): pooled 0.056 | per-gene-sum 0.128 (biased up) (true 0.05)
```
Structural conclusions hold for REML; only the "HE=REML" justification was wrong. The proxy is faithful
because the off-diagonal's only loss (the AxA diagonal) carries little info once additive is partialled
out, so proxy c=1 sits just below the measured 3-comp REML c=1.12.

## Task 5 — GO/NO-GO end-to-end 3-component AI-REML (feasibility_sim.py)  => **GO**

Design n=1500, G=50, Me_within=45.4, Me_within/n=0.030 (matches full-scale operating point 0.029).
200 reps, true h2_A=0.30, h2_AA=0.05:
```
detection CRB (sigma->0) = 0.0302 (c=1.00); EXACT CRB at h2=0.05 = 0.0270 (just past boundary)
mean h2_AA_hat = 0.0460  median 0.0452  bias -0.0040 (~2 MC-SE)
empirical SD   = 0.0261  vs exact CRB 0.0270  -> ratio 0.97  (REML ATTAINS the CRB)
boundary hits  = 8/200 (4%)
empirical power= 0.58    vs analytic 0.58  (exact)
```
**Effective c (3-component within-gene REML) = Me_within/(n*SE) = 1.12**, sim-verified.
RECONCILIATION: check_knobs.py's c~2 is a 2-component single-AxA-kernel MoM where the AxA DIAGONAL
carries info comparable to the off-diagonal (Var(Q_ij)~2Var(G_ij)^2), doubling information. That c
does NOT transfer: adding the additive component partials info away, and a POOLED within-gene kernel
has a different diagonal/off-diagonal split. c is regime- and model-dependent; the sim gives the
right value (1.12) for the actual model. So 1.6 was the wrong plug-in.

## Task 4 — required N, regime, sensitivity (required_n.py, sim-calibrated c=1.12)

k = z(.95)+z(.80) = 2.486. **Headline (Me_genome=74k, h2_AA=0.05): N ~ 23,200 (low-LD) to 27,500
(high-LD, with sqrt of the 1.1-1.4 skew partialling penalty).** Genome-wide 3.7M; gain sqrt(G)=141.
- h2_AA=0.02 -> N ~ 58k-65k (still inside UKB 500k).
- Knob band (Me_genome 20k-250k, penalty 1.0-1.4): N ~ 6k-93k. All inside UK Biobank.
- Regime ratio Me_within/N ~ 0.023 < h2=0.05 -> just past detection boundary; sim confirms REML
  attains the exact CRB there (so the c=1.12 / exact-CRB scaling is the right basis, not c=1.6).
- CAVEAT: Me_genome~74k is a UKB common-SNP literature/data quantity, NOT recomputed here. Reported
  as a band, with provenance to be cited in the note.

## Adversarial scaling stress-test (adversarial_scaling.py) — feasibility SURVIVES

```
(1) N-scaling: exact CRB*n at h2_AA=0.05, n=500..3000 = 41.5,40.9,40.5,40.6,41.3,42.0  -> SE ∝ 1/N. ✓
(2) G-scaling: Me_within/(sqrt(G)*Me_gene), G=25..1600:
      independent       1.006, 0.999, 1.000, 0.999   (algebraic identity)
      between-LD rho=.3 0.996, 0.994, 0.995, 0.997   (LD does NOT break the sqrt(G) law) ✓
(3a) REML vs exact CRB:  n=1500 ratio 0.97 ; n=2500 ratio 1.13 (unbiased, mean 0.0509)
     -> REML attains CRB to ~10% (finite-sample inefficiency); CRB-based N is mildly optimistic.
(3b) NULL type-I (LRT, h2_AA=0, crit 2.706, 300 reps) = 0.063 (target 0.05) -> power is REAL.
```
**Refined headline: N ≈ 25–30k** for 80% power at total cis h2_AA=0.05 (CRB-based 23–27k + ~10% REML
inefficiency). Still inside UK Biobank. The numbers DO scale; not too good to be true.

---
## PART I VERIFICATION COMPLETE — feasibility CONFIRMED (N ~ 25-30k, inside UK Biobank).
Three design claims corrected by the checks: (1) §3 partialling NOT free (genotype-skew, 1.1-1.4x);
(2) §5 the per-gene SUM is the fragile one (double-counts under between-gene LD), pooled is robust;
(3) §6 effective c=1.12 not 1.6, and REML ~10% above CRB -> N~25-30k not 18k.
Adversarial: SE∝1/N and sqrt(G) law both confirmed to large scale; null type-I 0.063 (power real).
All headline numbers re-derived in-repo.

## Architecture / sparsity power sims (power_sparsity.py) — 2026-06-24

Figure-1 variation: detection power vs N for the within-gene pooled estimator under varying AxA
architecture, variance-matched at total h2_AA=0.05. Causal fraction pi swept (k = pi*P causal pairs).
Test uses the infinitesimal kernel W; signal g=H*gamma with sparse gamma (Bernoulli(pi)).

P=750 (m=6, within_gene_power_sparsity.png):
```
n     typeI  k=750  k=75   k=15   k=4
500   .056   .125   .138   .162   .182
1200  .055   .419   .417   .399   .395
2500  .056   .941   .898   .810   .590
4000  .046  1.000   .995   .919   .721
```
P=2250 (m=10, within_gene_power_sparsity_largeP.png):
```
n     typeI  k=2250 k=225  k=45   k=11
1200  .056   .160   .184   .191   .222
2500  .049   .517   .486   .479   .453
4000  .056   .886   .853   .796   .651
```
Findings: (i) NULL is exactly architecture-free (type-I ~0.05 for all pi). (ii) Power curves CROSS
near power~0.4: variance-matching equalizes the test's non-centrality (mean), so only the statistic's
tail differs -- heavier tail (sparse) helps when underpowered, hurts when well-powered. (iii) The
binding quantity is the causal COUNT k, not the fraction pi: k~4 lags badly at large N, k in the dozens
tracks infinitesimal. So detection is architecture-robust at/below the operating point for realistically
polygenic genes; only a handful of causal pairs costs power, and only past the detection boundary.
Refines note caveat (iv).
