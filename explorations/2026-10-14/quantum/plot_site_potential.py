"""Figure: the exact single-site potential of a relaxed hypercycle row (nearest-neighbour model),
h(phi) = 2P ln sin(phi) + V cot(phi), V = P sin(2 alpha), in units of P, for several curvatures kappa = cos(alpha);
minimum on the row's hypercycle (phi = alpha), barrier on the conjugate one (phi = pi/2 - alpha); and the barrier
Delta E / P = 2[ln cot(alpha) - cos(2 alpha)] vs kappa with its cubic law near kappa_c. Usage: python plot_site_potential.py"""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 9}); COL = ['#2a78d6', '#1baf7a', '#eda100', '#eb6834']
fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.9))
ph = np.linspace(0.12, 1.55, 600)
for c, kappa in zip(COL, [0.95, 0.85, 0.75, 1/np.sqrt(2)]):
    a = np.arccos(kappa); h = 2*np.log(np.sin(ph)) + np.sin(2*a)/np.tan(ph); h0 = 2*np.log(np.sin(a)) + np.sin(2*a)/np.tan(a)
    lab = f'$\\kappa$ = {kappa:.3f}' + (' ($\\kappa_c$)' if abs(kappa - 2**-0.5) < 1e-9 else '')
    ax[0].plot(ph/np.pi*180, h - h0, color=c, lw=2, label=lab)
    ax[0].plot([a/np.pi*180], [0], 'o', color=c, ms=5); ax[0].plot([90 - a/np.pi*180], [2*(np.log(1/np.tan(a)) - np.cos(2*a))], 's', color=c, ms=5)
ax[0].set_xlabel('$\\varphi$ (degrees; row at $\\alpha$ = o, barrier at $90° - \\alpha$ = s)'); ax[0].set_ylabel('$[h(\\varphi) - h(\\alpha)]/P$')
ax[0].set_ylim(-0.3, 1.2); ax[0].set_title('single-site potential: minimum and barrier sit at conjugate curvatures', fontsize=9)
k = np.linspace(0.7072, 0.99, 300); a = np.arccos(k)
dE = 2*(np.log(1/np.tan(a)) - np.cos(2*a))
ax[1].plot(k, dE, color=COL[0], lw=2, label='exact $2[\\ln\\cot\\alpha - \\cos 2\\alpha]$')
ax[1].plot(k, (16/3)*(np.pi/4 - a)**3, color=COL[3], lw=1.5, ls='--', label='$(16/3)(\\pi/4 - \\alpha)^3$')
ax[1].set_xlabel('geodesic curvature $\\kappa$'); ax[1].set_ylabel('barrier $\\Delta E / P$ per vortex')
ax[1].set_title('per-vortex barrier vanishes cubically at $\\kappa_c = 1/\\sqrt{2}$', fontsize=9)
for a_ in ax: a_.grid(color='#e6e5e0', lw=0.6); a_.spines[['top', 'right']].set_visible(False); a_.legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig('site_potential.png', dpi=150); print('ok')
