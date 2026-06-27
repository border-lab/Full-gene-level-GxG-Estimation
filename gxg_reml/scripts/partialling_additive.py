"""
Detection-limit (sigma^2 -> 0) CRB for sigma^2_AA in a 3-component model
(additive A, within-gene AxA W, residual I) vs the 2-component model (W, I).
At sigma^2 -> 0, I_kl ∝ tr(G_k G_l), so CRB_W = (Gram^{-1})_{WW}.
Predict penalty ≈ 1/(1 - r^2), r = partial corr of W,A off-diagonals (I handles the mean/diagonal).
The cross term tr(AW) ≈ sum A_ij^3 (third moment of relatedness) ≈ 0 -> penalty ≈ 1.

Uses Ziyan's chr1 panels: contiguous = one high-LD block (stress case), random = dispersed.
"""
import numpy as np, pandas as pd

def load_raw(path, n_sub=3000, m=150):
    df = pd.read_csv(path, sep=r"\s+", nrows=n_sub)
    X = df.iloc[:, 6:6 + m].to_numpy(float)
    col = np.nanmean(X, 0); bad = np.where(np.isnan(X)); X[bad] = np.take(col, bad[1])
    return X[:, X.std(0) > 0]

def kernels(X):
    n = X.shape[0]
    Z = (X - X.mean(0)) / X.std(0); M = Z.shape[1]
    A = Z @ Z.T / M
    P = np.eye(n) - np.ones((n, n)) / n
    W = P @ (A * A) @ P                      # centered Hadamard square = implicit AxA, centered
    A = A * (n / np.trace(A)); W = W * (n / np.trace(W))
    return A, W

def crb_W(Gs, idx):
    Gram = np.array([[np.sum(Gi * Gj) for Gj in Gs] for Gi in Gs])
    return np.linalg.inv(Gram)[idx, idx]

print("=== (A) WORST CASE: cis-additive vs cis-AxA (A,W from the SAME SNPs) ===")
print(f"{'panel':11s} {'v3/v2':>8s} {'1/(1-r^2)':>10s} {'corr_od':>8s} {'mean A_ij^3':>12s}")
for name, path in [("contiguous", "/tmp/geno_contig/chr1_Contiguous1KSNP.raw"),
                   ("random",     "/tmp/geno_random/chr1_Random1KSNP.raw")]:
    A, W = kernels(load_raw(path)); n = A.shape[0]; I = np.eye(n)
    v3 = crb_W([A, W, I], 1); v2 = crb_W([W, I], 0)
    od = ~np.eye(n, dtype=bool)
    r = np.corrcoef(A[od], W[od])[0, 1]
    print(f"{name:11s} {v3/v2:8.4f} {1/(1-r**2):10.4f} {r:+8.3f} {np.mean(A[od]**3):+12.2e}")

print("\n=== (B) POOLED penalty vs WITHIN-GENE LD and G ===")
print("    genes = AR(1) LD blocks (corr rho within gene, independent across genes);")
print("    A = additive GRM over all G genes, W = pooled cis AxA. Claim: penalty is")
print("    G-independent and set by within-gene LD (= per-gene penalty).")
rng = np.random.default_rng(0)
P_ = lambda n: np.eye(n) - np.ones((n, n)) / n
norm = lambda K: K * (K.shape[0] / np.trace(K))

def ar1_genes(n, G, m, rho):
    C = np.linalg.cholesky(np.array([[rho**abs(i - j) for j in range(m)] for i in range(m)]))
    Z = np.hstack([rng.standard_normal((n, m)) @ C.T for _ in range(G)])
    Z = (Z - Z.mean(0)) / Z.std(0)
    return Z, np.repeat(np.arange(G), m)

def pen(Z, gene, G):
    n = Z.shape[0]; Pn = P_(n); I = np.eye(n); od = ~np.eye(n, dtype=bool)
    A = norm(Z @ Z.T / Z.shape[1])
    W = norm(sum(norm(Pn @ ((Zg := Z[:, gene == g]) @ Zg.T / Zg.shape[1])**2 @ Pn)
                 for g in range(G)) / G)
    v3 = crb_W([A, W, I], 1); v2 = crb_W([W, I], 0)
    return v3 / v2, np.corrcoef(A[od], W[od])[0, 1]

print(f"{'rho':>5s} {'G':>4s} {'pooled v3/v2':>13s} {'corr_od':>8s}   {'per-gene(G=1) v3/v2':>20s}")
for rho in (0.0, 0.6, 0.9):
    p1, r1 = pen(*ar1_genes(2500, 1, 6, rho), 1)          # single gene
    for G in (20, 80):
        pG, rG = pen(*ar1_genes(2500, G, 6, rho), G)
        print(f"{rho:5.1f} {G:4d} {pG:13.4f} {rG:+8.3f}   {p1:20.4f}")
print("  => penalty is a GENOTYPE-SKEW effect: 1.0 for symmetric (Gaussian) features,")
print("     >1 only for real (skewed) genotypes (case A). G-independent.")

print("\n=== (C) DILUTION on real data: cis-gene AxA vs additive GRM = cis SNPs + T trans ===")
print("    W = cis AxA of a high-LD gene (contiguous block, m=150).")
print("    A = additive GRM from those 150 cis SNPs + T trans SNPs (random panel).")
def stdize(X):
    col = np.nanmean(X, 0); bad = np.where(np.isnan(X)); X = X.copy(); X[bad] = np.take(col, bad[1])
    X = X[:, X.std(0) > 0]; return (X - X.mean(0)) / X.std(0)
Zc = stdize(pd.read_csv("/tmp/geno_contig/chr1_Contiguous1KSNP.raw", sep=r"\s+", nrows=3000
                        ).iloc[:, 6:156].to_numpy(float))
Rand = pd.read_csv("/tmp/geno_random/chr1_Random1KSNP.raw", sep=r"\s+", nrows=3000)
n = Zc.shape[0]; Pn = P_(n); I = np.eye(n); od = ~np.eye(n, dtype=bool)
Acis = Zc @ Zc.T / Zc.shape[1]
W = norm(Pn @ (Acis * Acis) @ Pn)
print(f"{'T trans':>8s} {'cis frac':>9s} {'v3/v2':>8s} {'corr_od':>8s}")
for T in (0, 150, 450, 850):
    if T == 0:
        Zall = Zc
    else:
        Zt = stdize(Rand.iloc[:, 6:6 + T].to_numpy(float)); Zall = np.hstack([Zc, Zt])
    A = norm(Zall @ Zall.T / Zall.shape[1])
    v3 = crb_W([A, W, I], 1); v2 = crb_W([W, I], 0); r = np.corrcoef(A[od], W[od])[0, 1]
    print(f"{T:8d} {Zc.shape[1]/(Zc.shape[1]+T):9.2f} {v3/v2:8.4f} {r:+8.3f}")
