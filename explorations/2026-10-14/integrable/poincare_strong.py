"""Discriminating test: N = 3 Poincare sections at LARGE amplitude for the vortex row (phi = 2/(1+Y^2)) vs the control
phi = 1/(1+Y^2)^2. Near-integrable systems develop chaotic layers far from equilibrium; exactly integrable ones do not.
Usage: python poincare_strong.py -> poincare_strong.png, poincare_strong.out"""
import numpy as np
from phi_family import rhs_phi
from poincare_phi import section, N, L
PHIS = {'2/(1+Y^2)  (vortex row)': (lambda Y: 2/(1 + Y**2), lambda Y: -4*Y/(1 + Y**2)**2),
        '1/(1+Y^2)^2  (control)': (lambda Y: 1/(1 + Y**2)**2, lambda Y: -4*Y/(1 + Y**2)**3)}
rng = np.random.default_rng(9); Y0 = -1.2
ics = [(Y0 + rng.uniform(-1.2, 1.2), Y0 + rng.uniform(-1.2, 1.2), Y0 + rng.uniform(-1.2, 1.2), L/3 + rng.uniform(-0.5, 0.5)) for _ in range(24)]
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4)); cols = plt.cm.viridis(np.linspace(0, 0.95, len(ics))); out = open('poincare_strong.out', 'w')
for a, (name, (phi, dphi)) in zip(ax, PHIS.items()):
    secs = section(rhs_phi(phi, dphi, L), ics, 6000.0)
    for c, sct in zip(cols, secs):
        if len(sct) > 20: a.plot(sct[:, 0], sct[:, 1], '.', ms=0.8, color=c)
    a.set_title(f'phi = {name}, large amplitude', fontsize=9); a.set_xlabel('$D_2$ at $D_1 = L/3$'); a.set_ylabel('$Y_1 - Y_2$')
    a.spines[['top', 'right']].set_visible(False)
    line = f'{name}: crossings per trajectory {[len(s) for s in secs]}'; print(line, flush=True); out.write(line + '\n')
fig.tight_layout(); fig.savefig('poincare_strong.png', dpi=140); print('saved')
