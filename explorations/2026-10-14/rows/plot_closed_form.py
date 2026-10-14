"""Exact hypercycle-row threshold vs the elementary closed form F0 = 0 (dual_series.py) and the proved bound l_A.
Usage: python plot_closed_form.py -> closed_form.png"""
import numpy as np, matplotlib, mpmath as mp
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from scipy.optimize import brentq
from dual_series import Fpi_dual, F0 as F0mp

def Fpi(kappa, lp):
    a = np.arccos(kappa); s2 = np.sin(a)**2; n = np.arange(1, int(np.ceil(90/lp)) + 6); r = np.exp(-n*lp)
    return np.sum((np.sin(2*a)*np.sin(2*n*a) - 2*n*s2*np.sin(n*a)**2)*r/(1-r) - 2*n*s2*np.cos(n*a)**2*r/(1+r))
def F0(kappa, lp):
    a = np.arccos(kappa); eta = np.pi*a/lp
    return (np.pi*np.sin(2*a)/(2*lp))/np.tanh(2*eta) - (np.pi**2*np.sin(a)**2/(4*lp**2))/np.tanh(eta)**2 - a*np.sin(2*a)/lp - np.cos(2*a)/4
def root(f, kappa, lo=0.005, hi=40):
    g = np.geomspace(lo, hi, 800); v = np.array([f(kappa, x) for x in g]); i = np.where(np.diff(np.sign(v)) > 0)[0]
    return brentq(lambda x: f(kappa, x), g[i[-1]], g[i[-1]+1], xtol=1e-14) if len(i) else np.nan

ks = np.concatenate([np.linspace(0.7085, 0.8, 40), np.linspace(0.8, 0.9995, 50)[1:]])
ex, c0 = [], []
for k in ks:
    s = np.sqrt(1-k**2); ex.append(root(Fpi, k)/s); c0.append(root(F0, k, hi=12)/s)
ex, c0 = np.array(ex), np.array(c0)
# high-precision gap where it is below double precision
gap = []
for k, e in zip(ks, ex):
    s = mp.sqrt(1-mp.mpf(k)**2)
    if k < 0.86: gap.append(abs(e - c0[list(ks).index(k)])/e); continue
    E = mp.findroot(lambda l: Fpi_dual(k, l), (e*s*0.98, e*s*1.02), solver='anderson')
    R = mp.findroot(lambda l: F0mp(k, l), (E*0.98, E*1.0000001), solver='anderson')
    gap.append(float(abs(E-R)/E))
gap = np.array(gap)
lA = np.arccosh(1/(2*ks**2-1))/np.sqrt(1-ks**2)

INK, MUTED, GRID = '#2b2b29', '#6b6a64', '#e4e3df'; C1, C2, C3 = '#2a78d6', '#eb6834', '#1baf7a'
fig, (ax, bx) = plt.subplots(1, 2, figsize=(11.5, 4.8), gridspec_kw={'width_ratios': [1.25, 1]})
ax.plot(ks, ex, '-', color=C1, lw=2.4, label='exact threshold (root of F(pi))')
ax.plot(ks, c0, '--', color=C2, lw=1.8, label='elementary closed form F0 = 0  (F < F0)')
ax.plot(ks, lA, ':', color=C3, lw=2, label='proved upper bound l_A')
ax.axvline(1/np.sqrt(2), color=MUTED, lw=1, ls='--')
ax.plot([1.0], [0.8439085], 's', color=INK, ms=7, label='horocycle row, l* = 0.8439')
ax.set_yscale('log'); ax.set_ylim(0.7, 8); ax.set_xlim(0.7, 1.005)
ax.set_yticks([0.8, 1, 2, 4, 8]); ax.set_yticklabels(['0.8', '1', '2', '4', '8']); ax.yaxis.set_minor_formatter(NullFormatter())
ax.set_xlabel('geodesic curvature kappa of the hypercycle', color=INK, fontsize=9.5)
ax.set_ylabel('critical arc spacing l* (curvature radii)', color=INK, fontsize=9.5)
ax.legend(frameon=False, fontsize=8.5, labelcolor=INK, loc='upper right')
ax.set_title('Stable above the blue curve; the closed form is the dashed curve', fontsize=9.5, color=INK, loc='left')
m = ks <= 0.985
bx.semilogy(ks[m], np.maximum(gap[m], 1e-60), '-', color=C2, lw=2)
bx.set_xlabel('kappa', color=INK, fontsize=9.5); bx.set_ylabel('relative gap (l* - l*_F0)/l*', color=INK, fontsize=9.5)
bx.set_xlim(0.7, 0.985); bx.set_ylim(1e-38, 1)
bx.set_title('Error of the closed form: about exp(-2 pi (pi - alpha)/l\')', fontsize=9.5, color=INK, loc='left')
for a_ in (ax, bx):
    a_.grid(True, which='major', color=GRID, lw=0.8); a_.set_axisbelow(True)
    for sp_ in ['top', 'right']: a_.spines[sp_].set_visible(False)
    for sp_ in ['left', 'bottom']: a_.spines[sp_].set_color(MUTED)
    a_.tick_params(colors=MUTED, labelsize=8.5)
plt.tight_layout(); plt.savefig('closed_form.png', dpi=140)
print('saved; relative gap at kappa = 0.8:', gap[np.argmin(abs(ks-0.8))])
