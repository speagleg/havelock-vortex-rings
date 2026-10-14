"""Figure: iDMRG (chi = 16) of the pinned quantised row at kappa = 0.98, l = 0.95 l*: zigzag order parameter and
correlation length vs quantum strength tau (scan_tau_*.out). Usage: python plot_dmrg.py -> dmrg_scan.png"""
import re, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
pts = {}
for f, taus in [('scan_tau_coarse.out', None), ('scan_tau_fine.out', None), ('scan_tau_window.out', [0.0076, 0.0073, 0.0079])]:
    vals = [l for l in open(f) if l.startswith('tau')]
    for i, l in enumerate(vals):
        t = taus[i] if taus else float(re.search(r'tau ([\d.]+)', l).group(1))
        eta = abs(float(re.search(r'eta\s+(-?[\d.]+)', l).group(1))); xi = float(re.search(r'xi\s+([\d.]+)', l).group(1))
        pts[t] = (eta, xi)
t = np.array(sorted(pts)); eta = np.array([pts[x][0] for x in t]); xi = np.array([pts[x][1] for x in t])
plt.rcParams.update({'font.size': 9})
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
ax[0].plot(t*1e3, eta, color='#2a78d6', lw=2, marker='o', ms=5); ax[0].set_ylabel('zigzag order $|\\eta|$')
ax[1].plot(t*1e3, xi, color='#eb6834', lw=2, marker='o', ms=5); ax[1].set_ylabel('correlation length $\\xi$ (sites)')
for a in ax:
    a.set_xlabel('quantum strength $\\tau = k\\varepsilon^2$  ($\\times 10^{-3}$)'); a.grid(color='#e6e5e0', lw=0.6)
    a.spines[['top', 'right']].set_visible(False); a.axvline(7.9, color='#6b6a64', lw=0.8, ls=':')
fig.suptitle('Pinned quantised row, $\\kappa$ = 0.98, $\\ell = 0.95\\,\\ell^*$ (iDMRG, $\\chi$ = 16): quantum melting of the zigzag at $\\tau_c \\approx 0.0079$', fontsize=9)
fig.tight_layout(); fig.savefig('dmrg_scan.png', dpi=150); print('ok', list(zip(t, eta, xi)))
