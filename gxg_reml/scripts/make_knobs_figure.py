"""
Figure: SE of h2_AA detection vs N, with the marker/LD/relatedness dependence explicit.
  genome-wide AxA, REML:  SE = M_e/(2.2 N)      (Hadamard-square -> Var(G_AA)=2/M_e^2)
  additive,        REML:  SE = sqrt(2 M_e)/(1.22 N)
  1/M_e = 1/N + 1/M + rho2_LD + Var(relatedness)
Hivert Eq6/7: M_e=74k, N=255k -> AxA HE(joint) 0.29 / REML(joint CRB) 0.18 / REML(meta 8x32k) 0.25; additive ~0.001.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

N = np.logspace(4, 7, 400)
# c = M_e/(N*SE): joint full-information / REML at biobank scale.  From Hivert Eq 7
# (a=Var(Gii)~2 Var(Gij)) the joint CRB gives c~1.6 (SE 0.18 at N=255k, M_e=74k).
se_axa  = lambda Me: Me/(1.6*N)
se_add  = lambda Me: np.sqrt(2*Me)/(1.22*N)
Me0 = 74000                      # UKB common SNPs, unrelated

fig, ax = plt.subplots(figsize=(9.4, 6.4))

# AxA band over M_e range (LD/relatedness lower it; more common markers raise it)
ax.fill_between(N, se_axa(2.0e4), se_axa(2.5e5), color='#1f77b4', alpha=0.12, zorder=1)
ax.plot(N, se_axa(Me0), '-',  color='#1f77b4', lw=2.4, zorder=4,
        label=r'genome-wide AxA detection  ($M_e\approx74$k)')
ax.plot(N, se_axa(2.0e4), '--', color='#1f77b4', lw=1.2, alpha=0.8, zorder=3)
ax.plot(N, se_axa(2.5e5), ':',  color='#1f77b4', lw=1.2, alpha=0.8, zorder=3)
ax.plot(N, se_add(Me0), '-', color='gray', lw=2.0, zorder=4,
        label=r'additive, REML  ($M_e\approx74$k, reference)')

# Hivert point: joint CRB (on the curve) vs what they actually ran (meta 8x32k, SE 0.25)
ax.scatter(254679, 74000/(1.6*254679), s=170, marker='*', color='#1f77b4',
           edgecolor='k', lw=0.6, zorder=8)
ax.scatter(254679, 0.25, s=80, marker='D', facecolor='none', edgecolor='#1f77b4', lw=1.4, zorder=8)
ax.text(254679*0.92, 0.182*1.5, 'joint REML = CRB\n(N=255k)', fontsize=8.0,
        ha='right', color='#1f77b4')
ax.text(254679*1.12, 0.25, 'Hivert (meta 8$\\times$32k)', fontsize=7.6, va='center', color='#1f77b4')
ax.scatter(254679, np.sqrt(2*74000)/(1.22*254679), s=60, marker='o', color='gray', zorder=8)
ax.text(254679*1.1, 0.00124*0.78, 'additive REML\n(SE~0.001)', fontsize=8, color='gray')

# detection thresholds (80% power, one-sided a=5%: SE < effect/2.49)
for h2, y in [(0.05, 0.05/2.49), (0.02, 0.02/2.49)]:
    ax.axhline(y, ls=':', color='#aa0000', lw=1.0)
    ax.text(1.05e4, y*1.1, f'detect $h^2_{{AA}}={h2}$', fontsize=8.5, color='#aa0000')

# band labels: what moves M_e
ax.text(8.5e6, se_axa(2.5e5)[-1]*1.05, 'more common markers $\\uparrow M_e$\n(+ better tagging via WGS)',
        fontsize=7.6, color='#1f77b4', ha='right', va='bottom')
ax.text(8.5e6, se_axa(2.0e4)[-1]*0.80, 'LD / relatedness $\\downarrow M_e$',
        fontsize=7.6, color='#1f77b4', ha='right', va='top')

# formula box (mathtext-safe: plain slashes, no \dfrac/\underbrace)
ax.text(1.1e4, 0.0019,
        r'$\mathrm{SE}(\hat h^2_{AA}) \approx M_e\,/\,(c\,N)$,   $c\approx1.6$' '\n'
        r'$1/M_e = 1/N + 1/M + \overline{\rho^2} + \mathrm{Var}_{\mathrm{rel}}$' '\n'
        'MoM $=$ REML here (small $h^2$);\n'
        'LD / relatedness lower $M_e$, lower SE',
        fontsize=9.2, va='top', ha='left',
        bbox=dict(boxstyle='round', fc='white', ec='0.6', alpha=0.95))

# N annotations on the AxA curve (c=1.6, matching se_axa)
for Nv, lab in [(5e5,'0.5M'),(1e6,'1M'),(2e6,'2M'),(5e6,'5M')]:
    ax.scatter(Nv, 74000/(1.6*Nv), s=28, color='#1f77b4', zorder=6)
    ax.text(Nv, 74000/(1.6*Nv)*0.62, lab, fontsize=7.5, color='#1f77b4', ha='center')

ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(1e4, 1e7); ax.set_ylim(5e-4, 1.0)
ax.set_xlabel('sample size  $N$', fontsize=11)
ax.set_ylabel(r'standard error of $\hat h^2_{AA}$  (phenotypic var. $=1$)', fontsize=11)
ax.set_title(r'Detection SE vs $N$, and its dependence on markers ($M$), LD, relatedness (via $M_e$)', fontsize=11.5)
ax.grid(True, which='both', ls='-', alpha=0.13)
ax.legend(loc='upper right', fontsize=8.8, framealpha=0.95)
plt.tight_layout()
plt.savefig('se_knobs_figure.pdf'); plt.savefig('se_knobs_figure.png', dpi=160)
print('saved se_knobs_figure.pdf/.png')
