"""Vortex rings and rows on the sphere (K = +1), the plane (K = 0) and the hyperbolic plane (K = -1, tanh model), in one
pair of coordinates: geodesic curvature kappa of the carrying curve, and arc spacing l between neighbours (curvature
radius 1). Sphere: latitude rings, kappa = cot(theta), l = 2 pi/(N sqrt(1 + kappa^2)), stable iff kappa >= kappa_N
(sphere_formula.py). Plane: rings, kappa l = 2 pi/N, stable iff N <= 7 (Thomson; N = 7 marginal). H^2: rings
(kappa = coth r > 1, l = 2 pi/(N sqrt(kappa^2 - 1))): N <= 7 stable at every radius, N >= 8 stable for kappa <= kappa*_N
(certified thresholds, ../rows/circle_branch.out); hypercycle rows (kappa < 1): stable iff l >= l*(kappa) (root of the exact
criterion F(pi), ../rows/qpi_proof.py), never for kappa <= 1/sqrt 2.
Usage: python plot_three_geometries.py -> three_geometries.png"""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from scipy.optimize import brentq

def Fpi(kappa, lp):
    a = np.arccos(kappa); s2 = np.sin(a)**2; n = np.arange(1, int(np.ceil(90/lp)) + 6); r = np.exp(-n*lp)
    return np.sum((np.sin(2*a)*np.sin(2*n*a) - 2*n*s2*np.sin(n*a)**2)*r/(1-r) - 2*n*s2*np.cos(n*a)**2*r/(1+r))
def lstar(kappa):
    g = np.geomspace(0.005, 40, 600); v = np.array([Fpi(kappa, x) for x in g]); i = np.where(np.diff(np.sign(v)) > 0)[0][-1]
    return brentq(lambda x: Fpi(kappa, x), g[i], g[i+1])/np.sqrt(1 - kappa**2)

INK, MUTED, GRID, SHADE = '#2b2b29', '#6b6a64', '#e4e3df', '#eef3fb'
C = {'S': '#eb6834', 'P': '#6b6a64', 'H': '#2a78d6'}
XL, YL = (0.25, 12), (0.08, 12)
fig, axs = plt.subplots(1, 3, figsize=(14, 4.9), sharey=True)
kk = np.geomspace(XL[0], XL[1], 600)

# sphere
ax = axs[0]; kN = {2: 0, 3: 0, 4: 1/np.sqrt(2), 5: 1.0, 6: 2.0}
for N in range(2, 9):
    l = 2*np.pi/(N*np.sqrt(1 + kk**2))
    ax.plot(kk, l, '--', color=C['S'], lw=0.9, alpha=0.5)
    if N in kN:
        m = kk >= kN[N]; ax.plot(kk[m], l[m], '-', color=C['S'], lw=2.4)
    i = 5 if N <= 6 else 5
    ax.annotate(f'N={N}', (kk[i], 2*np.pi/(N*np.sqrt(1 + kk[i]**2))), textcoords='offset points', xytext=(-2, 3), fontsize=7.5, color=INK)
for N, k in [(4, 1/np.sqrt(2)), (5, 1.0), (6, 2.0)]:
    ax.plot([k], [2*np.pi/(N*np.sqrt(1 + k*k))], 'o', color=C['S'], ms=5, mec='white')
ax.set_title('Sphere (K = +1): latitude rings\nstable only near the poles, only N <= 6', fontsize=9.5, color=INK, loc='left')

# plane
ax = axs[1]
for N in range(2, 11):
    l = 2*np.pi/(N*kk)
    ax.plot(kk, l, '-' if N <= 6 else (':' if N == 7 else '--'), color=C['P'], lw=2.2 if N <= 6 else 1.2, alpha=1 if N <= 7 else 0.5)
    ax.annotate(f'N={N}', (kk[30], 2*np.pi/(N*kk[30])), textcoords='offset points', xytext=(-2, 3), fontsize=7.5, color=INK)
ax.set_title('Plane (K = 0): rings, no length scale\nstable iff N <= 7 at every radius', fontsize=9.5, color=INK, loc='left')

# hyperbolic
ax = axs[2]
kh = np.concatenate([np.linspace(0.7085, 0.8, 30), np.linspace(0.8, 0.9995, 40)[1:]])
lh = np.array([lstar(k) for k in kh])
ax.fill_between(kh, lh, YL[1], color=SHADE, zorder=0)
ax.plot(kh, lh, '-', color=C['H'], lw=2.4, label='hypercycle rows: stable above')
ax.axvline(1/np.sqrt(2), color=MUTED, lw=1, ls='--')
ax.plot([1.0], [0.8439085], 's', color=INK, ms=6, label='horocycle row, l* = 0.8439')
cir = np.loadtxt('../rows/circle_branch.out', skiprows=1)
kc = np.geomspace(1.0005, XL[1], 500)
for N in range(2, 17):
    l = 2*np.pi/(N*np.sqrt(kc**2 - 1))
    if N <= 7:
        ax.plot(kc, l, '-', color=C['H'], lw=1.6, alpha=0.8)
    else:
        row = cir[cir[:, 0] == N]
        if not len(row): continue
        ks = row[0, 1]; m = kc <= ks
        ax.plot(kc[~m], l[~m], '--', color=C['H'], lw=0.9, alpha=0.45); ax.plot(kc[m], l[m], '-', color=C['H'], lw=1.6, alpha=0.8)
        ax.plot([ks], [row[0, 2]], 'o', color=C['H'], ms=4, mec='white')
ax.text(0.27, 0.11, 'hypercycle rows,\nkappa <= 1/sqrt2:\nnever stable', fontsize=8, color=INK)
ax.legend(frameon=False, fontsize=8, loc='upper right', labelcolor=INK)
ax.set_title('Hyperbolic plane (K = -1, tanh model): rings and rows\nstable once the spacing is large enough (kappa > 1/sqrt 2)', fontsize=9.5, color=INK, loc='left')

for ax in axs:
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(*XL); ax.set_ylim(*YL)
    ax.set_xticks([0.3, 0.5, 1/np.sqrt(2), 1, 2, 5, 10]); ax.set_xticklabels(['0.3', '0.5', '1/sqrt2', '1', '2', '5', '10'])
    ax.set_yticks([0.1, 0.3, 1, 3, 10]); ax.set_yticklabels(['0.1', '0.3', '1', '3', '10'])
    ax.xaxis.set_minor_formatter(NullFormatter()); ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel('geodesic curvature kappa of the curve', color=INK, fontsize=9)
    ax.grid(True, which='major', color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ['top', 'right']: ax.spines[s].set_visible(False)
    for s in ['left', 'bottom']: ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
axs[0].set_ylabel('arc spacing l between neighbours (curvature radii)', color=INK, fontsize=9)
fig.text(0.01, 0.005, 'Solid: linearly stable; dashed: unstable. Plane: units arbitrary (only kappa*l = 2pi/N matters).', fontsize=8, color=MUTED)
plt.tight_layout(rect=(0, 0.03, 1, 1)); plt.savefig('three_geometries.png', dpi=140)
