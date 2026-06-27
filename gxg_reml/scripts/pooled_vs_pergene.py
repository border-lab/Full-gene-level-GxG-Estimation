"""
Two estimators of total cis h2_AA, detection regime (HE = REML at small h2):
  (a) pooled : single AxA component, kernel W = mean_g K_g (assumes UNIFORM per-gene scaling)
  (b) per-gene aggregate: HE on each K_g, sum the per-gene estimates (adapts to any weighting)
Signal is generated from S = norm(sum_g w_g K_g): uniform weights vs sparse (few active genes).
Report bias (mean - h2) and SD for each estimator.
Claim: uniform -> both ~unbiased, comparable SD; sparse -> pooled BIASED (misspecified weights),
per-gene unbiased.
"""
import numpy as np
rng = np.random.default_rng(0)
n, G, m, maf = 2500, 40, 6, 0.3
P = np.eye(n) - np.ones((n, n)) / n
norm = lambda K: K * (n / np.trace(K))
od = ~np.eye(n, dtype=bool)

X = rng.binomial(2, maf, size=(n, G * m)).astype(float); gene = np.repeat(np.arange(G), m)
keep = X.std(0) > 0; X, gene = X[:, keep], gene[keep]
Z = (X - X.mean(0)) / X.std(0)
Kg = [norm(P @ ((Zg := Z[:, gene == g]) @ Zg.T / Zg.shape[1])**2 @ P) for g in range(G)]
W = norm(sum(Kg) / G)                      # pooled kernel (uniform), fixed

def he_factory(K):
    # off-diagonal HE slope: (y'Ky - sum_i K_ii y_i^2) / (sum K^2 - sum K_ii^2), no outer product
    d = np.diag(K).copy(); den = np.sum(K * K) - np.sum(d**2)
    return lambda y: (y @ K @ y - d @ (y * y)) / den

he_W = he_factory(W); he_g = [he_factory(K) for K in Kg]

def run(weights, h2=0.05, reps=400):
    w = np.array(weights, float); w = w / w.sum() * G
    S = norm(sum(wi * K for wi, K in zip(w, Kg)))
    LS = np.linalg.cholesky(h2 * S / (np.trace(S) / n) + 1e-9 * np.eye(n))
    ep, eg = [], []
    for _ in range(reps):
        y = LS @ rng.standard_normal(n) + np.sqrt(1 - h2) * rng.standard_normal(n)
        y -= y.mean()
        ep.append(he_W(y)); eg.append(sum(h(y) for h in he_g))
    return np.mean(ep), np.std(ep), np.mean(eg), np.std(eg)

print(f"=== INDEPENDENT genes (G={G}, n={n}, true h2=0.05) ===")
print(f"{'signal':16s} {'pool mean':>9s} {'pool SD':>8s} {'perg mean':>9s} {'perg SD':>8s}")
for label, w in [("uniform", np.ones(G)),
                 ("sparse(5 genes)", np.r_[np.ones(5), np.zeros(G - 5)]),
                 ("sparse(1 gene)", np.r_[np.ones(1), np.zeros(G - 1)])]:
    mp, sp, mg, sg = run(w)
    print(f"{label:16s} {mp:9.4f} {sp:8.4f} {mg:9.4f} {sg:8.4f}")
print("=> for independent genes, pooled-HE and per-gene-sum-HE are the SAME estimator (identical).")

# --- between-gene LD: split the contiguous block (one LD region) into adjacent windows ---
print(f"\n=== BETWEEN-GENE LD: contiguous block split into G2 adjacent windows ===")
import pandas as pd
Xc = pd.read_csv("/tmp/geno_contig/chr1_Contiguous1KSNP.raw", sep=r"\s+", nrows=n).iloc[:, 6:].to_numpy(float)
col = np.nanmean(Xc, 0); bad = np.where(np.isnan(Xc)); Xc[bad] = np.take(col, bad[1])
Xc = Xc[:, Xc.std(0) > 0]; Zc = (Xc - Xc.mean(0)) / Xc.std(0)
G2 = 40; win = np.array_split(np.arange(Zc.shape[1]), G2)
Kg2 = [norm(P @ ((Zw := Zc[:, idx]) @ Zw.T / Zw.shape[1])**2 @ P) for idx in win]
W2 = norm(sum(Kg2) / G2)
heW2 = he_factory(W2); heg2 = [he_factory(K) for K in Kg2]
# cross-gene off-diagonal overlap that breaks the identity:
cross = sum(np.sum(Kg2[a] * Kg2[b]) - np.diag(Kg2[a]) @ np.diag(Kg2[b])
            for a in range(G2) for b in range(G2) if a != b)
within = sum(np.sum(K * K) - np.diag(K) @ np.diag(K) for K in Kg2)
print(f"  cross/within off-diag overlap = {cross/within:.3f} (0 => identical estimators)")
S2 = norm(sum(Kg2) / G2)
LS2 = np.linalg.cholesky(0.05 * S2 / (np.trace(S2) / n) + 1e-9 * np.eye(n))
ep, eg = [], []
for _ in range(400):
    y = LS2 @ rng.standard_normal(n) + np.sqrt(0.95) * rng.standard_normal(n); y -= y.mean()
    ep.append(heW2(y)); eg.append(sum(h(y) for h in heg2))
ep, eg = np.array(ep), np.array(eg)
print(f"  pooled : mean={ep.mean():.4f} SD={ep.std():.4f}")
print(f"  per-gene: mean={eg.mean():.4f} SD={eg.std():.4f}  (differ when cross-overlap > 0)")
