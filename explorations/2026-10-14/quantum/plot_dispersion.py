"""Figure: chiral phonon dispersion E(q)/eps of a quantised hypercycle row. Left: kappa = 0.9 approaching the threshold
(soft mode at q = pi, Doppler-tilted). Right: near the critical curvature (kappa = 0.7075), where half the zone is almost
flat, with the leading-order kappa_c form sqrt8 s^3 e^{-l'} (|sin q| - sin q). Usage: python plot_dispersion.py"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rows'))
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from symbol import symbol
from exact_symbol import lstar
COL = ['#2a78d6', '#eb6834', '#1baf7a']; INK, MUTED, GRID = '#1f1f1e', '#6b6a64', '#e6e5e0'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False})
def E(kappa, l, q):
    s = np.sqrt(1 - kappa**2); B, C, dl = symbol(kappa, l*s, q); return np.sqrt(np.maximum(-B*C, 0)) + dl
q = np.linspace(-np.pi, np.pi, 1201)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
ls = lstar(0.9)
for c, r, xq in zip(COL, [1.001, 1.2, 2.0], [0.42, 0.8, 0.55]):
    y = E(0.9, ls*r, q); ax[0].plot(q/np.pi, y, color=c, lw=2, label=f'$\\ell/\\ell^*$ = {r}')
    ax[0].annotate(f'{r}', (xq, np.interp(xq*np.pi, q, y)), color=INK, fontsize=8, xytext=(4, 2), textcoords='offset points')
ax[0].set_title('$\\kappa$ = 0.9: the zigzag mode $q=\\pi$ softens; the branch is tilted (chiral)', fontsize=9, color=INK)
k = 0.7075; s = np.sqrt(1 - k*k); l = lstar(k)*1.2; lp = l*s
y = E(k, l, q)/np.exp(-lp)
ax[1].plot(q/np.pi, y, color=COL[0], lw=2, label=f'exact, $\\kappa$ = {k}, $\\ell = 1.2\\,\\ell^*$')
ax[1].plot(q/np.pi, np.sqrt(8)*s**3*(np.abs(np.sin(q)) - np.sin(q)), color=COL[1], lw=1.5, ls='--',
           label='leading order at $\\kappa_c$: $\\sqrt{8}s^3(|\\sin q| - \\sin q)$')
ax[1].set_title('near $\\kappa_c = 1/\\sqrt{2}$: half the zone is (almost) flat', fontsize=9, color=INK)
for a, yl in zip(ax, ['$E(q)/\\hbar_{\\rm eff}$  ($\\Gamma = 2\\pi$, $R = 1$)', '$E(q)/(\\hbar_{\\rm eff}\\,e^{-\\ell\'})$']):
    a.set_xlabel('Bloch wavenumber $q/\\pi$'); a.set_ylabel(yl); a.grid(color=GRID, lw=0.6); a.axhline(0, color=MUTED, lw=0.8)
    a.legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig('dispersion.png', dpi=150); print('ok')
