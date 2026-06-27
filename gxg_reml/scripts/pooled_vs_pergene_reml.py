"""
Does the pooled-vs-per-gene story hold for ACTUAL REML (not just the off-diagonal HE of eq.10)?
  pooled REML      : fit 3-component (A, W, I), h2_AA = s_W/sum(s)
  per-gene REML sum: for each gene g fit 3-component (A, K_g, I), sum_g s_{K_g}/sum(s_g)
Configs: (i) independent genes, (ii) between-gene LD (one contiguous LD block split into windows).
Report mean (bias) and SD of total h2_AA for each. HE claimed they are identical for independent
genes and that the per-gene SUM is biased up under between-gene LD; check both for REML directly.
"""
import numpy as np, pandas as pd
rng = np.random.default_rng(5)
norm = lambda K: K * (K.shape[0] / np.trace(K))
h2_A, h2_AA, reps = 0.30, 0.05, 40

def reml(y, Ks, n, iters=14):
    one = np.ones(n); I = np.eye(n); s = np.full(len(Ks), y.var() / len(Ks))
    for _ in range(iters):
        V = sum(si * Ki for si, Ki in zip(s, Ks)); Vi = np.linalg.inv(V + 1e-8 * I)
        Vio = Vi @ one; Pm = Vi - np.outer(Vio, Vio) / (one @ Vio); Py = Pm @ y
        sc = np.array([0.5 * (Py @ Ki @ Py - np.sum(Pm * Ki)) for Ki in Ks])
        AI = np.array([[0.5 * (Py @ Ki @ Pm @ Kj @ Py) for Kj in Ks] for Ki in Ks])
        s = np.clip(s + np.linalg.solve(AI + 1e-8 * np.eye(len(Ks)), sc), 1e-9, None)
    return s

def run(A, Kg, n, tag):
    W = norm(sum(Kg) / len(Kg)); I = np.eye(n)
    LA = np.linalg.cholesky(A + 1e-8 * I); LW = np.linalg.cholesky(W + 1e-8 * I)
    pool, perg = [], []
    for _ in range(reps):
        y = np.sqrt(h2_A) * (LA @ rng.standard_normal(n)) + np.sqrt(h2_AA) * (LW @ rng.standard_normal(n)) \
            + np.sqrt(1 - h2_A - h2_AA) * rng.standard_normal(n)
        y -= y.mean(); y /= y.std()
        s = reml(y, [A, W, I], n); pool.append(s[1] / s.sum())
        tot = 0.0
        for K in Kg:
            sg = reml(y, [A, K, I], n); tot += sg[1] / sg.sum()
        perg.append(tot)
    pool, perg = np.array(pool), np.array(perg)
    print(f"{tag}: pooled mean={pool.mean():.4f} SD={pool.std():.4f} | "
          f"per-gene-sum mean={perg.mean():.4f} SD={perg.std():.4f}  (true {h2_AA})", flush=True)

# (i) independent genes
n, G, m, maf = 1000, 8, 6, 0.3
X = rng.binomial(2, maf, (n, G * m)).astype(float); gene = np.repeat(np.arange(G), m)
keep = X.std(0) > 0; X, gene = X[:, keep], gene[keep]; Z = (X - X.mean(0)) / X.std(0)
P = np.eye(n) - np.ones((n, n)) / n
A = norm(Z @ Z.T / Z.shape[1])
Kg = [norm(P @ ((Zg := Z[:, gene == g]) @ Zg.T / Zg.shape[1])**2 @ P) for g in range(G)]
print("=== (i) INDEPENDENT genes (REML) ===", flush=True)
run(A, Kg, n, "indep")

# (ii) between-gene LD: contiguous block split into windows
n2, G2 = 1200, 8
Xc = pd.read_csv("/tmp/geno_contig/chr1_Contiguous1KSNP.raw", sep=r"\s+", nrows=n2).iloc[:, 6:].to_numpy(float)
col = np.nanmean(Xc, 0); bad = np.where(np.isnan(Xc)); Xc[bad] = np.take(col, bad[1])
Xc = Xc[:, Xc.std(0) > 0]; Zc = (Xc - Xc.mean(0)) / Xc.std(0)
P2 = np.eye(n2) - np.ones((n2, n2)) / n2
A2 = norm(Zc @ Zc.T / Zc.shape[1])
win = np.array_split(np.arange(Zc.shape[1]), G2)
Kg2 = [norm(P2 @ ((Zw := Zc[:, idx]) @ Zw.T / Zw.shape[1])**2 @ P2) for idx in win]
print("=== (ii) BETWEEN-GENE LD: contiguous block in 8 windows (REML) ===", flush=True)
run(A2, Kg2, n2, "between-LD")
