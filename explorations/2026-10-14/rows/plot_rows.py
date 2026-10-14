"""Stability boundary of vortex rows on curves of geodesic curvature kappa in H^2 (tanh model).
Reads rows_boundary.txt (hypercycle branch, from rows.py) and circle_branch.out (certified ring thresholds).
Dashed: the proved sufficient bound of stability_bound.py (Theorem A), l_A = arccosh(1/(2k^2-1))/sqrt(1-k^2).
Usage: python plot_rows.py -> rows_boundary.png"""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt

hyp = np.loadtxt('rows_boundary.txt')
cir = np.loadtxt('circle_branch.out', skiprows=1)
INK, MUTED, GRID = '#2b2b29', '#6b6a64', '#e4e3df'
C1, C2, C3 = '#2a78d6', '#eb6834', '#1baf7a'

fig, ax = plt.subplots(figsize=(8.6, 5.4))
ax.fill_between([0.5, 1/np.sqrt(2)], 0.3, 6, color='#f1f0ec', zorder=0)
ax.axvline(1/np.sqrt(2), color=MUTED, lw=1, ls='--', zorder=1)
kk = np.linspace(0.7075, 0.99999, 400)
ax.plot(kk, np.arccosh(1/(2*kk**2-1))/np.sqrt(1-kk**2), '--', color=C1, lw=1.4, alpha=0.75, zorder=2,
        label='proved: stable above this (hypercycles)')
ax.plot(hyp[:, 0], hyp[:, 1], '-o', color=C1, lw=2, ms=6, zorder=3, label='hypercycle rows (computed)')
ax.plot([1.0], [0.8439085], 's', color=C2, ms=9, zorder=4, mec='white', mew=1.5, label='horocycle row (note: l* = 0.8439)')
ax.plot(cir[:, 1], cir[:, 2], '-^', color=C3, lw=2, ms=7, zorder=3, label='rings N = 8..64 (certified thresholds)')
for N, k, l in cir:
    if int(N) in (8, 10, 12, 16, 64):
        ax.annotate(f'N={int(N)}', (k, l), textcoords='offset points', xytext=(6, -12), fontsize=8, color=INK)
ax.text(0.515, 4.6, 'proved unstable\nat every spacing\n(kappa <= 1/sqrt 2)', fontsize=9, color=INK, va='top')
ax.text(1/np.sqrt(2)+0.006, 0.32, 'kappa = 1/sqrt 2  (exact critical curvature)', fontsize=8.5, color=MUTED, rotation=90, va='bottom')
ax.text(1.35, 2.2, 'stable above the curve', fontsize=10, color=INK)
ax.text(1.35, 0.36, 'unstable below', fontsize=10, color=INK)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(0.5, 2.2); ax.set_ylim(0.3, 6)
ax.set_xticks([0.5, 0.6, 0.7, 0.8, 1, 1.25, 1.5, 2]); ax.set_xticklabels(['0.5', '0.6', '0.7', '0.8', '1', '1.25', '1.5', '2'])
ax.set_yticks([0.3, 0.5, 1, 2, 4]); ax.set_yticklabels(['0.3', '0.5', '1', '2', '4'])
from matplotlib.ticker import NullFormatter
ax.yaxis.set_minor_formatter(NullFormatter()); ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xlabel('geodesic curvature kappa of the curve carrying the row  (0 geodesic, 1 horocycle, >1 circle)', color=INK, fontsize=9.5)
ax.set_ylabel('critical spacing l* along the curve (curvature radii)', color=INK, fontsize=9.5)
ax.grid(True, which='major', color=GRID, lw=0.8); ax.set_axisbelow(True)
for s in ['top', 'right']: ax.spines[s].set_visible(False)
for s in ['left', 'bottom']: ax.spines[s].set_color(MUTED)
ax.tick_params(colors=MUTED, labelsize=8.5)
ax.legend(frameon=False, fontsize=8.5, loc='upper right', labelcolor=INK)
ax.set_title('Vortex rows on H$^2$ (tanh model): one stability boundary across geodesics, hypercycles, horocycles and circles',
             fontsize=10, color=INK, loc='left')
plt.tight_layout(); plt.savefig('rows_boundary.png', dpi=140)
