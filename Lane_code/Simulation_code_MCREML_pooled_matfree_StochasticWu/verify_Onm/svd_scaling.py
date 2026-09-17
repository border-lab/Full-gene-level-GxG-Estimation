"""Does the truncation rank s needed for a fixed accuracy grow like m or like sqrt(m)?

  s ~ sqrt(m)  ->  apply O(n s^2) = O(n m)      the goal
  s ~ c m      ->  apply O(n s^2) = O(c^2 n m^2)  only a constant-factor win

Also: why truncation is not subject to the rank-Nmc obstruction that kills the
sketch.  rank(K_s .* K_s) <= s(s+1)/2, so a rank-s truncation of K can still
produce a FULL-RANK n x n operator whenever s(s+1)/2 >= n.  The Hutchinson
sketch is hard-capped at rank Nmc.
"""
import numpy as np
from svd_truncation import sim_geno, W_centred, W_trunc

TOL = 0.15

print("=" * 88)
print(f"(A)  smallest s with relFro(W_s) <= {TOL}, as m grows at fixed n")
print("=" * 88)
print(f"{'rho':>5} {'n':>6} {'m':>6} {'s*':>6} {'s*/m':>7} {'s*/sqrt(m)':>11} "
      f"{'s*(s*+1)/2':>11} {'>=n?':>5} {'cost s*^2/m^2':>14}")
for rho in (0.99, 0.9):
    for n in (400,):
        for m in (50, 100, 200, 400):
            Z, _ = sim_geno(n, m, rho=rho, seed=1000 + m)
            Wex = W_centred(Z)
            nrm = np.linalg.norm(Wex)
            lo, hi = 1, m
            while lo < hi:                      # relFro is monotone in s
                mid = (lo + hi) // 2
                if np.linalg.norm(W_trunc(Z, mid) - Wex) / nrm <= TOL:
                    hi = mid
                else:
                    lo = mid + 1
            s = lo
            print(f"{rho:>5} {n:>6} {m:>6} {s:>6} {s/m:>7.3f} {s/np.sqrt(m):>11.2f} "
                  f"{s*(s+1)//2:>11} {'yes' if s*(s+1)//2 >= n else 'no':>5} "
                  f"{s*s/(m*m):>14.3f}")

print()
print("=" * 88)
print("(B)  same, as n grows at fixed m   (the sketch needed Nmc ~ n here)")
print("=" * 88)
print(f"{'rho':>5} {'m':>6} {'n':>6} {'s*':>6} {'s*/m':>7}  <- does s* track n?")
for rho in (0.99, 0.9):
    for m in (150,):
        for n in (200, 400, 800, 1600):
            Z, _ = sim_geno(n, m, rho=rho, seed=2000 + n)
            Wex = W_centred(Z)
            nrm = np.linalg.norm(Wex)
            lo, hi = 1, m
            while lo < hi:
                mid = (lo + hi) // 2
                if np.linalg.norm(W_trunc(Z, mid) - Wex) / nrm <= TOL:
                    hi = mid
                else:
                    lo = mid + 1
            print(f"{rho:>5} {m:>6} {n:>6} {lo:>6} {lo/m:>7.3f}")

print()
print("=" * 88)
print("(C)  effective rank of W_s vs the rank cap of an Nmc-probe sketch")
print("     (this is why truncation escapes the obstruction in 2.5)")
print("=" * 88)
n, m = 400, 200
Z, _ = sim_geno(n, m, rho=0.9, seed=7)
Wex = W_centred(Z)
lam = np.linalg.eigvalsh(Wex)
print(f"  exact W:  eff.rank tr(W)^2/tr(W^2) = "
      f"{lam.sum()**2/(lam**2).sum():.1f}   (n = {n})")
print(f"\n{'s':>5} {'algebraic rank cap':>20} {'eff.rank(W_s)':>15} {'relFro':>9}")
for s in (10, 20, 40, 80, 160):
    Ws = W_trunc(Z, s)
    l = np.linalg.eigvalsh(Ws)
    print(f"{s:>5} {min(s*(s+1)//2, n):>20} {l.sum()**2/(l**2).sum():>15.1f} "
          f"{np.linalg.norm(Ws-Wex)/np.linalg.norm(Wex):>9.4f}")
print(f"\n  a Nmc-probe Hutchinson sketch is capped at algebraic rank Nmc, and needs")
print(f"  Nmc >~ eff.rank(W) ~ O(n); truncation reaches rank n from s ~ sqrt(2n) = "
      f"{int(np.sqrt(2*n))}.")
