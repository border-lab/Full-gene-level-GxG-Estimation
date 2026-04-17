# Thesis Feedback: Scalable Computational Approaches for Epistatic Variance Component Estimation


---

## Fixes (do immediately)

- Heritability formula (p.12) is wrong: denominator should be $\sigma^2_g + \sigma^2_{g \times g} + \sigma^2_e$, not $\sigma^2_g + \sigma^2_{g \times g}$. (fixed)
- Table of contents: all Results subsections point to "page 2." (check page finnaly)

## Framing

- The space complexity improvement ($O(nm)$ vs $O(n^2)$) is the strongest selling point. An $n \times n$ GRM can't even be stored at biobank scale. The abstract mentions $O(nm^2)$ time but buries the space result. Foreground it. (add to abstract page v; emphasize in the text)
- The abstract and intro oversell. The method works under independence; for correlated SNPs (i.e., real data), SE doesn't decrease with n. Frame honestly: "works under independence, adjusted version removes bias for LD but consistency under strong LD remains open. (add to abstract page v; how to fix intro as we put forward the object)"
- The introduction doesn't engage with the counter-argument that statistical epistatic variance is expected to be small in outbred populations. If the quantity being estimated is near zero in humans, the SE requirements change. (I have add this at the section of Epistasis as a Contributor to Missing Heritability, it is a strong evidence why we need consistent estimation)
- specific about traits  wrt missing h2, gxg （add more example ? if it is, add at the end of The Challenge of Complex Trait Architecture and Epistasis as a Contributor to Missing Heritability）
- bological framing of epistasis

## Simulations: things to add

- **Realistic effect sizes.** Everything uses $\sigma^2_{g \times g} = 0.9$. Rerun at $\sigma^2_{g \times g} = 0.05$ or $0.1$. Can the estimator detect it? If not, say so. -- how big of a sample would you need for 80% power as a function of sample size and h2_gxg 
- **Include additive effects.** The model drops the additive component with a one-line justification. Run at least one simulation with $y = Z\beta + H\gamma + e$ and show what happens to the epistatic estimate when additive variance is present but unmodeled. 
- adding h2 across multiple genes 

## Results: things to explain

- The non-monotonic bias in sparse-SNP results is unexplained. At minimum offer a hypothesis.
- The counterintuitive finding that SE decreases with increasing m for contiguous SNPs is flagged but left without even speculation. A likely explanation: larger m dilutes LD structure by including more distantly-spaced SNPs. If true, this tells you something important about how the method interacts with LD blocks. Can you propose a potential experiment?

## Method concerns
- The adjusted estimator requires $\mu_{ij}$ and $\sigma^2_{ij}$ for all SNP pairs. Clarify: are these estimated from data or computed from known allele frequencies? If from data, what sample size is needed for stable pairwise LD estimates?

## Discussion

- missing contextualization within existing literature 
  - esp wrt to hivert findings, 
  - is hivert se approximation accurate wrt what you observe?
- Restructure with subheadings: contributions, limitations, future work.
- Spend more time on the LD/SE problems, which is the main open question
- Propose future real data analysis:
  - what data
  - what traits
  - computational feasibility

## Overall
need more hypotheses to explain mysterious findings, propose experiments to test these hypotheses


## meeting
1. The cholesky problem
2. Simulation code
3. add cite to support gene-level episatis