"""
Figure: SE of the epistatic variance estimate vs the effective rank the estimator
exploits.  Master floor SE = sqrt(2/r) (phenotypic variance standardized to 1).
The two estimators land at DIFFERENT r for the same data/N: REML reaches the true
resolved rank; MoM(HE) is throttled to r4 << rank under spectral skew.
Existing/projected results are placed at their (effective rank, SE).

All SE values are on the h2_AA scale (sigma_y^2 = 1).
Sources:
  - Hivert et al. 2021 (UKB, N=254,679): SE from their own Eq 6 (HE) and Eq 7 (REML),
    with b=Var(Gij)=1/74k (anchors SE_HE=0.29) and a=Var(Gii)~2b. r = 2/SE^2.
    HE joint=0.29; REML meta 8x32k (what they did)=0.25; REML joint/CRB=0.18.
  - REML SE scales ~1/N (detection regime), so r ~ N^2.  Projections at 0.5/1/2/5 M.
  - Gene-level array (Ziyan contiguous): r4~5, MoM floored at SE~0.4-0.6, no N helps.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

floor = lambda r: np.sqrt(2.0/r)          # SE floor for r independent dimensions
r_from_se = lambda se: 2.0/se**2

fig, ax = plt.subplots(figsize=(9.2, 6.4))
rr = np.logspace(0, 5.2, 400)
ax.plot(rr, floor(rr), 'k-', lw=2.2, zorder=5,
        label=r'fundamental floor  SE $=\sqrt{2/r}$  (attained by REML)')

# --- detectability thresholds (80% power, one-sided alpha=5%: SE < effect/2.49) ---
for eff, lab, y in [(0.05, r'detect $h^2_{AA}=0.05$', 0.05/2.49),
                    (0.02, r'detect $h^2_{AA}=0.02$', 0.02/2.49)]:
    ax.axhline(y, ls=':', color='gray', lw=1.1)
    ax.text(1.15, y*1.07, lab+f'  (SE<{y:.3f})', fontsize=8.5, color='gray', va='bottom')
ax.axhspan(0.05/2.49, 1.5, color='red', alpha=0.045, zorder=0)
ax.text(5.5e4, 0.74, 'cannot detect realistic epistasis (per trait)',
        fontsize=8.5, color='#aa0000', va='center', ha='right')

# --- existing & projected points ---
# Hivert genome-wide AxA at N=254,679 (from their Eq 6/7; a=Var(Gii)~2 Var(Gij))
H_he   = (r_from_se(0.29), 0.29)   # HE, joint on all pairs
H_meta = (r_from_se(0.25), 0.25)   # REML meta-analysed 8x32k -- what they actually did
H_reml = (r_from_se(0.18), 0.18)   # REML joint = Cramer-Rao bound at N=255k (NOT run)
ax.scatter(*H_he,   s=90, marker='s', color='#1f77b4', zorder=6)
ax.scatter(*H_meta, s=85, marker='D', facecolor='none', edgecolor='#1f77b4', lw=1.4, zorder=6)
ax.scatter(*H_reml, s=160, marker='*', color='#1f77b4', zorder=7, edgecolor='k', lw=0.6)
ax.annotate('', xy=H_reml, xytext=H_he,
            arrowprops=dict(arrowstyle='->', color='#1f77b4', lw=1.4))
ax.text(H_he[0]*0.78, H_he[1]*0.80, 'Hivert\nHE joint, N=255k', fontsize=8, ha='center', va='top', color='#1f77b4')
ax.text(H_meta[0]*0.62, H_meta[1]*1.16, 'REML meta 8$\\times$32k\n(what they did)', fontsize=7.3, ha='center', color='#1f77b4')
ax.text(H_reml[0]*1.15, H_reml[1]*1.16, 'REML joint = CRB\n(not run: needs\nmatrix-free solver)', fontsize=7.6, color='#1f77b4')

# REML projections vs N (r ~ N^2 in detection regime; anchored at the JOINT REML CRB)
N0, r0 = 254679, H_reml[0]
for N in (0.5e6, 1e6, 2e6, 5e6):
    r = r0*(N/N0)**2; se = floor(r)
    ax.scatter(r, se, s=70, marker='o', color='#2ca02c', zorder=6)
    ax.text(r*1.05, se*1.18, f'{N/1e6:g}M', fontsize=8.5, color='#2ca02c', ha='left')
ax.text(3.5e4, 0.0033, 'genome-wide AxA, REML,\nfuture biobanks  (SE $\\propto 1/N$)',
        fontsize=8.5, color='#2ca02c', ha='center')

# Gene-level array, HE (Ziyan) -- floored, on the curve at r4~5
gx = (5, floor(5))
ax.scatter(*gx, s=110, marker='X', color='#d62728', zorder=6)
ax.text(gx[0]*1.35, gx[1]*1.02, 'gene-level array, HE\n($r_4\\approx5$ — floored,\nno $N$ helps)',
        fontsize=8, color='#d62728', va='center', ha='left')
ax.annotate('', xy=(5e3, 0.020), xytext=(9, 0.52),
            arrowprops=dict(arrowstyle='->', color='#d62728', lw=1.3, ls='--', alpha=0.8))
ax.text(330, 0.115, 'WGS + REML\n(rank $\\uparrow$)', fontsize=8.2, color='#d62728', rotation=-21)

ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(1, 1.6e5); ax.set_ylim(2.5e-3, 1.0)
ax.set_xlabel('effective rank exploited,  $r$  (independent genetic dimensions)', fontsize=11)
ax.set_ylabel(r'standard error of $\hat h^2_{AA}$   (phenotypic var. $=1$)', fontsize=11)
ax.set_title('Precision of epistatic variance estimators vs. effective rank', fontsize=12)
ax.grid(True, which='both', ls='-', alpha=0.13)

handles = [Line2D([0],[0],color='k',lw=2.2,label=r'floor SE$=\sqrt{2/r}$ (REML-attainable)'),
           Line2D([0],[0],marker='*',color='w',markerfacecolor='#1f77b4',markeredgecolor='k',markersize=13,label='Hivert REML joint = CRB (not run)'),
           Line2D([0],[0],marker='D',color='w',markeredgecolor='#1f77b4',markersize=9,label='Hivert REML meta 8$\\times$32k (run)'),
           Line2D([0],[0],marker='s',color='w',markerfacecolor='#1f77b4',markersize=9,label='Hivert HE joint (all pairs)'),
           Line2D([0],[0],marker='o',color='w',markerfacecolor='#2ca02c',markersize=9,label='genome-wide AxA REML, projected $N$'),
           Line2D([0],[0],marker='X',color='w',markerfacecolor='#d62728',markersize=10,label='gene-level array (floored)')]
ax.legend(handles=handles, loc='lower left', fontsize=8.6, framealpha=0.95)

plt.tight_layout()
plt.savefig('rank_precision_figure.pdf')
plt.savefig('rank_precision_figure.png', dpi=160)
print("saved rank_precision_figure.pdf / .png")
# print the placed points for the writeup
print(f"Hivert HE joint   N=255k: r={H_he[0]:.0f}, SE=0.29")
print(f"Hivert REML meta  N=255k: r={H_meta[0]:.0f}, SE=0.25")
print(f"Hivert REML joint N=255k: r={H_reml[0]:.0f}, SE=0.18 (CRB)")
for N in (0.5e6,1e6,2e6,5e6):
    r=r0*(N/N0)**2; print(f"REML N={N/1e6:g}M: r={r:.0f}, SE={floor(r):.4f}")
