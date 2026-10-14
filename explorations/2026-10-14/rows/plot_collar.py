"""Necklace thresholds around the Bolza systole on the row stability diagram.
Reads rows_boundary.txt and collar.out (from collar.py). Usage: python plot_collar.py -> collar.png"""
import numpy as np, re, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

hyp = np.loadtxt('rows_boundary.txt')
pts = []
for line in open('collar.out'):
    m = re.search(r'^\s*(\d+).*kappa\* = ([\d.]+), spacing ([\d.]+)', line)
    if m: pts.append((int(m.group(1)), float(m.group(2)), float(m.group(3))))
INK, MUTED, GRID = '#2b2b29', '#6b6a64', '#e4e3df'
C1, C2 = '#2a78d6', '#eb6834'
kc_collar = np.sqrt(2) - 1

fig, ax = plt.subplots(figsize=(8.6, 5.2))
ax.axvspan(0.3, kc_collar, color='#f1f0ec', zorder=0)
ax.axvline(1/np.sqrt(2), color=MUTED, lw=1, ls='--', zorder=1)
ax.plot(hyp[:, 0], hyp[:, 1], '-o', color=C1, lw=2, ms=5, zorder=3, label='infinite row: critical spacing')
P = np.array(pts)
ax.plot(P[:, 1], P[:, 2], 's', color=C2, ms=8, mec='white', mew=1.2, zorder=4, label='necklace of N around the systole:\nfirst stable point')
for N, k, l in pts:
    if N in (2, 3, 4, 6, 8, 12):
        ax.annotate(f'N={N}', (k, l), textcoords='offset points', xytext=(7, 4), fontsize=8, color=INK)
ax.text(0.31, 5.2, 'inside the systole collar\n(kappa <= sqrt2 - 1 = 0.414):\nno stable necklace', fontsize=9, color=INK, va='top')
ax.text(1/np.sqrt(2)+0.006, 0.82, 'kappa = 1/sqrt 2', fontsize=8.5, color=MUTED, rotation=90, va='bottom')
ax.set_yscale('log'); ax.set_xlim(0.3, 1.0); ax.set_ylim(0.8, 6)
ax.set_yticks([1, 2, 3, 4, 5]); ax.set_yticklabels(['1', '2', '3', '4', '5']); ax.yaxis.set_minor_formatter(NullFormatter())
ax.set_xlabel('geodesic curvature kappa = tanh h  (h = distance from the closed geodesic)', color=INK, fontsize=9.5)
ax.set_ylabel('spacing along the necklace (curvature radii)', color=INK, fontsize=9.5)
ax.grid(True, which='major', color=GRID, lw=0.8); ax.set_axisbelow(True)
for s in ['top', 'right']: ax.spines[s].set_visible(False)
for s in ['left', 'bottom']: ax.spines[s].set_color(MUTED)
ax.tick_params(colors=MUTED, labelsize=8.5)
ax.legend(frameon=False, fontsize=8.5, loc='upper left', bbox_to_anchor=(0.2, 0.6), labelcolor=INK)
ax.set_title('Vortex necklaces around the Bolza systole (collar model, tanh kernel): every necklace inside the collar is unstable',
             fontsize=9.5, color=INK, loc='left')
plt.tight_layout(); plt.savefig('collar.png', dpi=140)
