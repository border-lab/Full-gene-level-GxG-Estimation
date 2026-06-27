"""
Verify the within-gene-pooling result (note sec:localize):
  M_e^genome = G * M_e^gene   and   M_e^within = sqrt(G) * M_e^gene,
so pooling within-gene AxA gains a factor sqrt(G) over the genome-wide all-pairs SE.
M_e of an AxA kernel: Var(offdiag) = 2/M_e^2 -> M_e = sqrt(2/Var).
"""
import numpy as np
rng = np.random.default_rng(0)
n, G, m, maf = 2000, 25, 8, 0.3            # G genes, m independent SNPs each
P = np.eye(n) - np.ones((n, n))/n
cen = lambda W: P @ W @ P
norm = lambda K: K * (n/np.trace(K))
od = ~np.eye(n, dtype=bool)
Me = lambda K: np.sqrt(2/K[od].var())

X = rng.binomial(2, maf, size=(n, G*m)).astype(float); gene = np.repeat(np.arange(G), m)
X = X[:, X.std(0) > 0]; gene = gene[:X.shape[1]]
Z = (X - X.mean(0))/X.std(0)
Kg = [norm(cen((lambda Zg: (Zg@Zg.T/Zg.shape[1])**2)(Z[:, gene == g]))) for g in range(G)]
Kwithin = norm(sum(Kg)/G)
Kgenome = norm(cen((Z@Z.T/Z.shape[1])**2))

Me_gene = np.mean([Me(K) for K in Kg])
print(f"M_e^gene   = {Me_gene:6.1f}   (predict ~m = {m})")
print(f"M_e^within = {Me(Kwithin):6.1f}   (predict sqrt(G)*Me_gene = {np.sqrt(G)*Me_gene:.1f})")
print(f"M_e^genome = {Me(Kgenome):6.1f}   (predict G*Me_gene = {G*Me_gene:.1f})")
print(f"gain genome/within = {Me(Kgenome)/Me(Kwithin):.2f}x   (predict sqrt(G) = {np.sqrt(G):.2f})")
