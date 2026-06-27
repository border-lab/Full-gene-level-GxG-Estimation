"""
Detection power for total cis h2_AA. SE(h2_AA) = Me_within/(c N), one-sided test H0: h2_AA=0.
Power beta at level alpha needs h2_AA/SE = z(1-alpha)+z(beta) =: k  =>  N = k Me_within/(c h2).
Me_within = Me_genome/sqrt(G). Me_genome is a UK Biobank common-SNP quantity (NOT recomputed
here): we report a sensitivity band and a partialling-penalty factor from partialling_additive.py.
"""
import numpy as np
from scipy.stats import norm
alpha, beta = 0.05, 0.80
k = norm.ppf(1 - alpha) + norm.ppf(beta)          # ~2.487
G = 20000
# c is NOT universal. feasibility_sim.py measured the effective constant for the ACTUAL
# 3-component within-gene REML at the operating point (h2_AA=0.05, just past detection boundary):
#   SE = 0.0270 at Me_within=45.4, n=1500  =>  c = Me_within/(n*SE) = 1.12  (REML attains CRB, ratio 0.97)
# check_knobs.py's c~2 was a 2-component single-AxA-kernel MoM (diagonal-information regime) -- not this model.
C_SIM = 1.12
reqN = lambda h2, Mew, c, pen: k * Mew * np.sqrt(pen) / (c * h2)   # SE inflates by sqrt(LD-skew penalty)

print(f"k = z(.95)+z(.80) = {k:.3f}; G = {G}, sqrt(G) = {np.sqrt(G):.1f}; c (sim-calibrated) = {C_SIM}")
print("\nrequired N for 80% power (c=1.12 from feasibility_sim; LD-skew partialling penalty 1.0-1.4 var):")
print(f"{'Me_genome':>9s} {'Me_within':>9s} {'h2_AA':>6s} {'N(pen1.0)':>10s} {'N(pen1.25)':>11s} {'genome-wide':>12s}")
for Meg in (20000, 74000, 250000):
    Mew = Meg / np.sqrt(G)
    for h2 in (0.05, 0.02):
        N0 = reqN(h2, Mew, C_SIM, 1.0); N1 = reqN(h2, Mew, C_SIM, 1.25)
        Ngw = reqN(h2, Meg, C_SIM, 1.25)
        print(f"{Meg:9d} {Mew:9.0f} {h2:6.2f} {N0:10,.0f} {N1:11,.0f} {Ngw:12,.0f}")

print("\nheadline (Me_genome=74k, h2=0.05): N =", f"{reqN(0.05, 74000/np.sqrt(G), C_SIM, 1.0):,.0f}",
      "(low-LD) to", f"{reqN(0.05, 74000/np.sqrt(G), C_SIM, 1.4):,.0f}", "(high-LD).")
print("knob band (Me_genome 20k-250k, penalty 1.0-1.4): N =",
      f"{reqN(0.05, 20000/np.sqrt(G), C_SIM, 1.0):,.0f}", "to",
      f"{reqN(0.05, 250000/np.sqrt(G), C_SIM, 1.4):,.0f}", "-- all inside UK Biobank (500k).")
print(f"regime: Me_within/N ~ {(74000/np.sqrt(G))/reqN(0.05,74000/np.sqrt(G),C_SIM,1.0):.3f} < h2=0.05 ->")
print("  just past detection boundary; feasibility_sim confirms REML attains the exact CRB there.")
