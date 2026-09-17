=== Estimating the four components

@realizedFour was a property of the simulator. This one is the estimator: the same
run, but now the *fitted* components, each centred on the expected realized variance
it should recover. Columns 1, 2, 4 and 5 of the result file are plotted;
column 3, the raw $hat(sigma)^2_(g times g)$, is not, because it lives on the
uncorrected $H$-component scale and only $c dot hat(V)_gamma$ is on a variance scale
(the same reason `calc_stats.py` reports no paired difference for gxg).

#figure(
image("estimate_four_components_m1000_G10_cexact.pdf", width: 100%),caption: [The four estimated variance components ($m = 1000$, $G = 10$, $R = 200$ replicates per box). Panels, top to bottom: $hat(sigma)^2_a$, $hat(sigma)^2_d$, $hat(sigma)^2_e$ and $c dot hat(V)_gamma$ --- columns 1, 2, 4 and 5 of the result file. The raw $hat(sigma)^2_(g times g)$ (column 3) is not shown: it lives on the uncorrected $H$-component scale, and only $c dot hat(V)_gamma$ is comparable to a variance. Each box is the deviation from the expected realized variance of its component ($sigma^2_a$, $sigma^2_d$, $sigma^2_e$, and $c dot sigma^2_(g times g)$ with the exact $c$); the dashed line is zero, and the label above each box is a one-sample $t$-test of that deviation against zero (`ns` / `*` / `**` / `***`) with its mean and SD. The estimator is unbiased throughout: $19$ of the $20$ boxes are `ns`, the exception being epistasis at $n = 4000$ (mean $+0.006$, `*`). SDs contract with $n$ in every panel --- additive $0.053 arrow.r 0.022$, dominance $0.039 arrow.r 0.013$, environment $0.038 arrow.r 0.009$, epistasis $0.073 arrow.r 0.023$ from $n = 1000$ to $n = 16000$. Six replicates across the grid have a component at `mc_reml`'s lower clamp and are retained, as in `calc_stats.py`; the panel scales are per-component, since the environment errors are about a third of the epistasis errors at the same $n$.]
) <estFour>

The contrast with @realizedFour is the point. There, only $hat(V)_e$ concentrated as
$n$ grew and the three genetic components did not shrink at all; here *every* panel
contracts, so the estimator does convert extra samples into precision even for the
components whose realized values do not settle. The rates differ, though. Fitting
$log "SD"$ on $log n$ gives $"SD" prop n^(-0.53)$ for the environment --- the $sqrt(n)$
rate --- against $n^(-0.29)$, $n^(-0.39)$ and $n^(-0.43)$ for additive, dominance and
epistasis. The genetic components buy precision more slowly than $sqrt(n)$, which is
consistent with the LD-limited effective dimension read off @realizedFour: with only
$approx 18$ independent additive directions in $m = 1000$ contiguous SNPs, adding rows
to $X$ eventually stops adding information about $sigma^2_a$.
